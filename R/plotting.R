#various functions for plotting phyloseq objects

library(ggplot2)

#Nice high contrast 20 colours
nice20 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',   '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',   '#000075', '#808080')
#larger set from rcolor brewer - a bit ugly but necessary
bigcolset<-c(RColorBrewer::brewer.pal(8,"Set2"), RColorBrewer::brewer.pal(12,"Set3"), RColorBrewer::brewer.pal(9,"Set1"), RColorBrewer::brewer.pal(12,"Paired"))
#a last resort!!!!!!
ridiculouslybigcolset <- c(rep(c(nice20,bigcolset),10))

#a nice theme for plotting purposes. can be added to the plots after they are made
theme_pp <- list(
  theme_bw(),
  theme(strip.background =element_rect(fill="white")),
  theme(axis.text.y = element_blank()),
  theme(axis.ticks = element_blank()),
  scale_x_discrete(expand = c(0,0)),
  scale_y_continuous(expand = c(0,0)),
  theme(panel.spacing.x = unit(0,"lines")))

#remove x labels from tgeme
theme_pp_blankx <- c(theme_pp,
                     theme(axis.text.x = element_blank()),
                     theme(axis.title.x = element_blank())
)

#' DADA2 output to physloseq object
#'
#' @param directory
#' @param metadata
#' @param ranks
#'
#' @return
#' @export
#'
#' @examples
ps_from_ampliseq <- function(directory,metadata=NULL,ranks){
  #asvs
  asv <- otu_table(as.matrix(read.table(paste0(directory,"/dada2/ASV_table.tsv"),header = TRUE,row.names = 1)),taxa_are_rows = TRUE)

  #taxonomy
  tax <- readr::read_tsv(paste0(directory,"/dada2/ASV_tax.tsv"),col_names = TRUE)
  asvnames <- tax$ASV_ID
  tax <- tax[,ranks]
  rownames(tax) <- asvnames
  tax <- tax_table(as.matrix(tax))

  if(!is.null(metadata)){
    #metadata - it's on you to format this properly.
    metadata <- sample_data(metadata)
  }

  ps <- phyloseq(asv,tax,metadata)

  return(ps)
}


#Strongly suspect this is bad practise
#' Return 20 pretty colours
#'
#' @return
#' @export
#'
#' @examples
get_nice_20 <- function(){
  return(c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',   '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',   '#000075', '#808080'))
}

#' Return a lot of ugly colours
#'
#' @return
#' @export
#'
#' @examples
get_big_colours <- function(){
  bigcolset<-c(RColorBrewer::brewer.pal(8,"Set2"), RColorBrewer::brewer.pal(12,"Set3"), RColorBrewer::brewer.pal(9,"Set1"), RColorBrewer::brewer.pal(12,"Paired"))
  #a last resort!!!!!!
  ridiculouslybigcolset <- c(rep(c(nice20,bigcolset),10))
  return(ridiculouslybigcolset)
  }

#to save time recalculating the agglomerated objects throughout the analysis.
#returns named list of phyloseq objects
#' Get Gloms
#'
#' @param ps
#' @param ranks
#'
#' @return
#' @export
#'
#' @examples
getgloms <- function(ps,ranks= c("Phylum" ,"Class","Order","Family","Genus")){
  ps@tax_table@.Data[is.na(ps@tax_table@.Data)] <- "Unassigned"
  gloms <- lapply(ranks, function(x) tax_glom(ps,taxrank=x,NArm=FALSE)) #agglomerate counts per taxon at each of the specied ranks
  names(gloms) <- ranks
  return(gloms)
}

#code for making abundance plots with phyloseq objects.
#by abundance specifies whether to order the taxa by decreasing ferquency in the plot.
#if no glom is supplied the function will calculate one instead
#' Plot Taxa Abundance
#'
#' @param ps
#' @param rank
#' @param xsep
#' @param wrap
#' @param n
#' @param colno
#' @param byabundance
#'
#' @return
#' @export
#'
#' @examples
plot_taxa_abundance <- function(ps,rank,xsep, wrap = NULL, n=20, colno = NULL,byabundance=TRUE,abs=FALSE,size=10){

  #set cols
  cols.n <- c(c(nice20, ridiculouslybigcolset)[1:n],"lightgrey")

  #check if phyloseq agglomerated at appropriate rank
  #and if not, perform agglomeration
  if(length(unique(as.data.frame(ps@tax_table@.Data)[,rank])) == length(unique(rownames((as.data.frame(ps@tax_table@.Data)))))){
    print("Agglomerating at specified rank...")
    ps@tax_table@.Data[is.na(ps@tax_table@.Data)] <- "Unassigned"
    glom <-tax_glom(ps, taxrank = rank, NArm = FALSE)
    print("Done.")
  }else{
    print("Using provided agglomerated object.")
  }

  #get the specified no of ASVs.
  topnotus <- names(sort(taxa_sums(glom),TRUE)[1:min(nrow(glom@tax_table),n)])

  #create a taxonomy table containing these ASVs with everything else set to "Other"
  taxtabn = cbind(tax_table(glom), taxon = "Other")
  taxtabn[topnotus, "taxon"] <- as(tax_table(glom)[topnotus, rank],
                                   "character")
  tax_table(glom) <- tax_table(taxtabn)
  melt <- psmelt(glom)

  #get names for reordering factors
  labels <- (unique(melt$taxon))
  nlabs <- length(labels)

  #move "Unassigned" to end
  if ("Unassigned" %in% labels){
    labels <- labels[!labels %in% "Unassigned"]
    labels[nlabs] = "Unassigned"
  }
  #if > n taxa, move "Other" to end
  if("Other" %in% melt$taxon){
    labels <- labels[!labels %in% "Other"]
    labels[nlabs] <- "Other"
  }
  for(i in seq(1:nlabs)){ #making other and unassigned have consistent colours
    if(labels[i] == "Unassigned"){
      cols.n[i] <- "grey"
    }
    else if(labels[i] == "Other"){
      cols.n[i] <- "lightblue"
    }
  }

  #labels by abundance w/ unassigned and otehr at end
  melt$taxon <- fct_relevel(melt$taxon,labels)

  #make the stacked bar chart
  i <- ggplot(melt, aes_string(x = xsep, y = "Abundance", fill = "taxon")) +
    geom_bar(stat = "identity", width = 1, position = position_fill()) +
    scale_fill_manual(values=cols.n, na.value = "grey")+
    theme(axis.title.x = element_blank())+
    labs(fill=rank)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "bottom")+
    theme(text = element_text(size = size))
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))

  if(abs){
    i <- i+ geom_bar(stat = "identity", width = 1, position = "stack")
  }

  if (is.character(wrap)){
    i <- i + facet_grid(as.formula(paste("~", wrap)), scales = "free_x", space = "free")
  }

  return(i)
}

#Bray curtis plots. PCOA is table from "ordinate(ps)".
#' Bray Plot
#'
#' @param ps
#' @param pcoa
#' @param colour
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
brayplot_pp <- function(ps, pcoa, colour, shape=NULL){
  p <- plot_ordination(ps, pcoa, color = colour, axes = c(1,2),shape = shape) +
    geom_point(size = 2) +
    labs(title = colour, color = colour)+
    scale_color_manual(values = nice20)


  return(p)
}

#alpha diverstiy plots. richtable is from estimate_richness(ps).
#' Alpha Diversity Plot
#'
#' @param richtable
#' @param xvar
#' @param yvar
#' @param wrapvar
#' @param colno
#'
#' @return
#' @export
#'
#' @examples
alphaplot_pp <- function(richtable, xvar, yvar, wrapvar = NULL, colno = NULL){
  p <- ggplot(richtable, aes_string(x = xvar, y = yvar, fill = xvar)) +
    geom_boxplot() +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none")

  if (is.character(wrapvar)){
    p <- p + facet_grid(as.formula(paste("~", wrapvar)),scales = "free_x")
  }
  return(p)
}

#for generating a table from kruskall significance tests. richtable is generated by estimate_richness(ps)
#' Do Kruskall Tests
#'
#' @param richtable
#' @param testvars
#' @param metrics
#'
#' @return
#' @export
#'
#' @examples
makekruskalltests <- function(richtable,testvars,metrics){
  #create df
  df <- data.frame()
  for (k in metrics) df[[k]] <- as.character()

  for (v in testvars){
    results <- vector()
    for(i in metrics){
      results[[i]] <- kruskal.test(formula(paste(i, "~ ", v)), data = richtable)$p.value
    }
    newrow <- nrow(df) + 1
    df[newrow,] <- results
    rownames(df)[newrow] <- v
  }
  return(df)
}

#plot expression scatter from DESEq2
deseq2_scatter <- function(sigtab,dds,x,wrap = NULL){
  plots <- list()
  for(t in rownames(sigtab)){
    title <- paste(sigtab[t,]$Family, sigtab[t,]$Genus,sep = "_")
    data <- plotCounts(dds = dds,gene = t,intgroup = colnames(dds@colData), main= title, returnData = TRUE)

    p <- ggplot(data, aes_string(x=x, y="count", color=x)) +
      scale_y_log10() +
      theme_classic() +
      geom_point(position=position_jitter(width=.1,height=0), size=5)+
      ggtitle(title) +
      theme(legend.position = "none") +
      theme(text = element_text(size = 10))

    if (is.character(wrap)){
      p <- p + facet_grid(as.formula(paste("~", wrap)), scales = "free_x", space="free")
    }
    plots[[t]] <- p
  }
  return(plots)
}

tax_to_box <- function(x,wrap,data){
  p <- ggplot(data, aes_string(x=x, y="count", fill=x)) +
    scale_y_log10() +
    theme_classic() +
    geom_boxplot()+
    geom_point(position=position_jitter(width=.1,height=0), size=5)+
    ggtitle(title) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 15))

  if (is.character(wrap)){
    p <- p + facet_grid(as.formula(paste("~", wrap)), scales = "free_x", space="free")
  }
  return(p)
}

#plot expression box from DESEq2
deseq2_box <- function(sigtab,dds,x,wrap = NULL){
  plots <- list()
  for(t in rownames(sigtab)){
    title <- paste(sigtab[t,]$Order,sigtab[t,]$Family, sigtab[t,]$Genus,sep = "_")
    data <- plotCounts(dds = dds,gene = t,intgroup = colnames(dds@colData), main= title, returnData = TRUE)

    p <- ggplot(data, aes_string(x=x, y="count", fill=x)) +
      scale_y_log10() +
      theme_classic() +
      geom_boxplot()+
      geom_point(position=position_jitter(width=.1,height=0), size=5)+
      ggtitle(title) +
      theme(legend.position = "none") +
      theme(text = element_text(size = 15))

    if (is.character(wrap)){
      p <- p + facet_grid(as.formula(paste("~", wrap)), scales = "free_x", space="free")
    }
    plots[[t]] <- p
  }
  return(plots)
}


sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

#function for plotting DA_bars from a sig table generated in DESeq2
da_bars <- function(sigtab, title,cols = NULL, x_limits = NULL,taxa_order = NULL,setcolor = FALSE,orderbyphylum = FALSE){
  bar_width =.8
  sigtabgen = subset(sigtab, !is.na(Genus) & Genus != "Unassigned" & Genus != "")

  # Specify custom order for Genus
  if (!is.null(taxa_order)) {
    sigtabgen <- sigtabgen[sigtabgen$Genus %in% taxa_order,]
    sigtabgen$Genus = factor(sigtabgen$Genus, (levels = rev(taxa_order)))

  } else if (orderbyphylum) {

    # Order the data by Phylum and then by log2FoldChange within each Phylum
    sigtabgen <- sigtabgen[rev(order(sigtabgen$Phylum, -sigtabgen$log2FoldChange)), ]

    # Create a custom factor order based on the sorted data
    genus_order <- sigtabgen$Genus

    # Reorder Genus factor levels based on the custom order
    sigtabgen$Genus <- factor(sigtabgen$Genus, levels = genus_order)

  } else {
    x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
    x = sort(x, TRUE)
    sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = rev(names(x)))
  }

  num_taxa = length(unique(sigtabgen$Genus))
  plot_height = max(2,num_taxa * 2.5)  # Minimum height of 3 units

  p <- ggplot(sigtabgen, aes(y = Genus, x = log2FoldChange, color = Phylum, fill = Phylum)) +
    theme_pubr() +
    geom_vline(xintercept = 0.0, color = "black", size = 1) +
    geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
    theme(axis.text.x = element_text(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 15)) +
    ggtitle(title) +
    guides(fill = guide_legend(override.aes = list(size = 5))) #+
  #coord_fixed(ratio = plot_height / num_taxa,expand = FALSE)

  if(! is.null(x_limits)){
    p <- p +
      xlim(x_limits)
  }

  if(!is.null(cols)){
    p <- p +
      scale_fill_manual(values = cols) +
      scale_color_manual(values = cols)
  }

  #this allows consistent colours across plots
  if(setcolor){
    p <- p +
      scale_fill_manual(values = c(Firmicutes = "orange",Actinobacteria="blue",Proteobacteria="purple",Bacteroidetes= "forestgreen")) +
      scale_color_manual(values = c(Firmicutes = "orange",Actinobacteria="blue",Proteobacteria="purple",Bacteroidetes= "forestgreen"))
  }

  return(p)
}

