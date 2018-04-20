##### This scripts consists of all the plotting functions #####

# Assigns colours
## check this for selecting colors http://colorbrewer2.org/

assignColours <- function(MSnSetObj, colourBy="SampleGroup"){
  if (!is(MSnSetObj, "MSnSet")){ stop("MSnSetObj has to be of class MSnSet..") }
  colourGroups <- as.character(pData(MSnSetObj)[,colourBy]) %>% sort() %>% unique()
  len <- length(colourGroups)
  if(len<3){
      sampleColours <- setNames(brewer.pal(3, "Dark2")[seq_len(len)], colourGroups)
  }else if(len<=8){
      sampleColours <- setNames(brewer.pal(len, "Dark2"), colourGroups)
  }else{
      coloursF <- brewer.pal(8, "Dark2") %>% colorRampPalette()
      sampleColours <- setNames(coloursF(len), colourGroups)
  }
  return(sampleColours)
}

# Intensity distribution plot
## intensities is a data frame containing columns for each sample
## sampleColours is a named vector that maps samples to colours
intensityPlot <- function(MSnSetObj, sampleColours=NULL, title="", colourBy="SampleGroup", 
                          transform=TRUE, xlab="log2(intensity)",
                          trFunc=log2xplus1){
  if (!is(MSnSetObj, "MSnSet")){ stop("MSnSetObj has to be of class MSnSet..") }
  if(!transform){ trFunc <- c }
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  exprs(MSnSetObj) %>% 
    as.data.frame() %>% 
    gather(SampleName, Intensity) %>%
    left_join(pData(MSnSetObj), "SampleName") %>% 
    mutate_at(vars(colourBy), funs(as.factor)) %>%
    mutate(Intensity=trFunc(Intensity)) %>% 
    filter(!is.na(Intensity)) %>%
    ggplot(aes_string(x="Intensity", group="SampleName", colour=colourBy)) + 
    stat_density(geom="line", position="identity") +
    scale_colour_manual(values=sampleColours, breaks=names(sampleColours)) +
    labs(x=xlab, title=title) +
    theme_bw() +
    theme(plot.title=element_text(size=14,hjust=0.5),
      axis.text.x=element_text(size=12),
      axis.title.x=element_text(size=14),
      axis.text.y=element_blank(),
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.title=element_blank(),
      legend.key=element_blank(),
      legend.key.width=unit(0.5, "cm"),
      legend.key.height=unit(0.5, "cm"),
      legend.text=element_text(size=11),
      legend.justification=c(1,1),
      legend.position=c(1,1),
      legend.background=element_rect(fill="transparent")
      )
}

# Peptide intensity plot
## intensities is a data frame containing peptide intensities with columns for each sample
## summarizedIntensities is a data frame containing summarized protein-level intensities
## protein is the protein for which intensities will be plotted
## samples is a list of samples to use in the plot

peptideIntensityPlot <- function(MSnSetObj, ProteinID, ProteinName, 
                                 combinedIntensities=NULL, selectedSequence=NULL, 
                                 selectedModifications=NULL){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!is.null(combinedIntensities) && !is(combinedIntensities,"MSnSet")){
      stop('combinedIntensities has to either NULL or object of class MSnSet..')
  }
  
  ## get peptide intensities for plotting
  intensities <- exprs(MSnSetObj) %>% 
    as.data.frame() %>% 
    rownames_to_column("PeptideID") %>%
    gather(SampleName, Intensity, -PeptideID) %>% 
    left_join(fData(MSnSetObj) %>% rownames_to_column("PeptideID"), "PeptideID") %>% 
    mutate(logIntensity=log2xplus1(Intensity)) %>% 
    filter(Accessions == ProteinID)
  
  ## get intensities for selected sequences if present
  seqIntensities <-  filter(intensities, Sequences %in% selectedSequence)
  if (!is.null(selectedModifications)){ 
    seqIntensities %<>% filter(Modifications %in% selectedModifications)
  }
  
  ## get intensity for protein level (summarised) data if present
  summIntensities <- data.frame(SampleName=vector(), logIntensity=vector(), Protein=vector())
  if (!is.null(combinedIntensities)){
    summIntensities <- exprs(combinedIntensities) %>% 
      as.data.frame() %>% 
      rownames_to_column("Protein") %>%
      gather(SampleName, Intensity, -Protein) %>% 
      left_join(fData(combinedIntensities), "Protein") %>% 
      mutate(logIntensity=log2xplus1(Intensity)) %>% 
      filter(Protein == ProteinID)
  }
  
  if (nrow(intensities) == 0) {
    warning("No peptides were found for ", ProteinID)
    return(NULL)
  }
  
  ggplot(intensities, aes(x=SampleName, y=logIntensity)) +
    geom_line(aes(group=PeptideID, colour=Sequences), size=0.6, alpha=0.5, linetype=2) +
    geom_point(aes(fill=Sequences), shape=21, size=2) +
    ## plot the selected sequences if present
    geom_line(data=seqIntensities, aes(group=PeptideID), colour="#6666FF", size=1.2, linetype=6)  +
    geom_point(data=seqIntensities, fill="#0000FF", shape=21, size=2.5) +
    ## plot the summarised intensities if present
    geom_line(data=summIntensities, aes(group=Protein), colour="#AAAAAA", size=1.2, linetype=6) +
    geom_point(data=summIntensities, fill="#888888", shape=21, size=2.5) +
    labs(y="log2(Intensity)", title=ProteinName) +
    theme_bw() +
    theme(plot.title=element_text(size=14, hjust=0.5),
          axis.text.x=element_text(size=11, angle=45, hjust=1),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.position="none")
}


# PCA plot
pcaPlot <- function(MSnSetObj, omitIgG=FALSE, sampleColours=NULL, transFunc=log2xplus1, transform=TRUE, 
                    colourBy="SampleGroup", title="", labelColumn="BioRep", labelsize=4, pointsize=4, x.nudge=4, x.PC=1){   
    
  if(!is(MSnSetObj,"MSnSet")){ stop('data has to be of class MSnSet..') }
  if(!colourBy%in%colnames(pData(MSnSetObj))){ stop("colourBy must a column names in the pData of the MSnSetObj")}
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  
  ## Remove IgG samples is requested
  if(omitIgG){ MSnSetObj <- MSnSetObj[,toupper(MSnSetObj$SampleGroup)!="IGG"] }
  if(!transform){ transFunc <- as.data.frame }
  intensities <- exprs(MSnSetObj) %>% as.data.frame() %>% na.omit() %>% transFunc()
  if (nrow(intensities) == 0){ return(NULL) }
  pca <- intensities %>% t() %>% prcomp()
  pcaVariance <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100)
  plotDat <- as.data.frame(pca$x) %>% 
    rownames_to_column("SampleName") %>% 
    left_join(pData(MSnSetObj), "SampleName") %>% 
    mutate_at(vars(colourBy), funs(as.factor))
  xPC <- paste0("PC", x.PC)
  yPC <- paste0("PC", x.PC+1)
  ggplot(plotDat, aes_string(x=xPC, y=yPC, fill=colourBy, label=labelColumn)) + 
    geom_point(pch=21, colour="black", size=pointsize) +
    { if(!is.null(labelColumn))
      geom_text(hjust=0, size=labelsize, nudge_x=x.nudge)
    } +
    scale_fill_manual(values=sampleColours, breaks=names(sampleColours)) +
    labs(x=paste0(xPC, ", ", pcaVariance[x.PC], "% variance"),
         y=paste0(yPC, ", ", pcaVariance[x.PC+1], "% variance"),
         fill=NULL, title=title) + 
    theme_bw() + 
    theme(text=element_text(size=14),
          plot.title=element_text(size=14,hjust=0.5),
          aspect.ratio=1)
}



# MA or volcano plot
maVolPlot <- function(diffstats, contrast, title="", controlGroup = NULL,
                      selectedGenes=NULL, fdrCutOff=0.05,
                      lfcCutOff=1, controlLfcCutOff=1, plotType="MA"){
  # For plotting we will assign the proteins to one of 7 groups:
  # A - selected (user specified in `selectedGenes`) - highlighted blue
  # B - significant - pass cutoffs - highligted red
  # C - non-significant - fail cutoffs - small and grey
  
  if(!plotType%in%c("MA", "Volcano")){ stop("plotType should be 'MA' or 'Volcano'..") }
  
  testSignficant <- function(dat){
    '%in%'(dat$adj.P.Val<=fdrCutOff, TRUE) & abs(dat$log2FC)>=lfcCutOff &
      abs(dat$controlLogFoldChange)>=controlLfcCutOff
  }
  
  daResTab <- suppressMessages(getContrastResults(diffstats=diffstats, contrast=contrast, 
                                                  controlGroup = controlGroup)) %>% 
    bind_rows(data.frame(controlLogFoldChange=vector())) %>% 
    mutate(controlLogFoldChange=ifelse(is.na(controlLogFoldChange), Inf, controlLogFoldChange)) %>% 
    mutate(group=ifelse(testSignficant(.), "Significant", "Non-significant")) %>%
    mutate(group=ifelse(GeneSymbol%in%selectedGenes, "Selected", group)) %>% 
    mutate(group=factor(group, levels=c("Selected", "Significant", "Non-significant"))) %>% 
    arrange(desc(group)) %>% 
    mutate(phredPval=-log10(adj.P.Val))
  
  if(plotType=="MA"){ 
    xFactor <- "AvgIntensity"; yFactor <- "log2FC" 
    xLab <- "average log2(Intensity)"; yLab <- "log2(Fold Change)"
  }
  if(plotType=="Volcano"){ 
    xFactor <- "log2FC"; yFactor <- "phredPval" 
    xLab <- "log2(Fold Change)"; yLab <- "-log10(Adjusted P value)"
  }
  
  xNudge <- diff(range(daResTab[,xFactor]))/100
  
  ggplot(daResTab, aes_string(x=xFactor, y=yFactor, 
                              colour="group", size="group", shape="group", alpha="group", fill="group")) +
    geom_hline(yintercept=0, color="gray50", size=0.5) +
    geom_point() +
      geom_text(data=subset(daResTab, group=="Selected"), aes(label=GeneSymbol), hjust=0, 
              size=3.5, nudge_x = xNudge) +
      scale_colour_manual(values=c("black", "black", "gray50"), drop=FALSE) +
      scale_size_manual(values=c(1.8, 1.5, 0.9), drop=FALSE) +
      scale_shape_manual(values=c(21, 21, 20), drop=FALSE) +
      scale_alpha_manual(values=c(1, 1, 0.6), drop=FALSE)  +
      scale_fill_manual(values=c("cyan", "red"), limits=levels(daResTab$group)[1:2], 
                        breaks=levels(droplevels(daResTab$group)), name="") +
      labs(x=xLab, y=yLab, title=title) + 
      theme_bw() +
      theme(
          text=element_text(size=16),
          plot.title=element_text(size=14,hjust=0.5)
      ) +
      guides(colour="none", size="none", shape="none", alpha="none",
             fill = guide_legend(override.aes = list(shape = 21)))
}

# Correlation plot
corrPlot <- function(MSnSetObj, addValues=TRUE, title=""){
  if(!is(MSnSetObj,"MSnSet")){ stop('data has to be of class MSnSet..') }
  if(!is.logical(addValues)){ stop('addValues has to be either TRUE or FALSE..') }
  
  col2Cols <- c("#FFFFFF", "#B90505")
  cor(exprs(MSnSetObj)) %>% 
    as.data.frame() %>% 
    rownames_to_column("X") %>% 
    gather("Y", "Cor", -1) %>% 
    mutate(addValues=addValues) %>% 
    mutate(CorTxt=ifelse(addValues==TRUE, round(Cor, 3), "")) %>%
    ggplot(aes(x=X, y=Y, fill=Cor)) +
    geom_tile(col="grey") +
    geom_text(aes(label=CorTxt)) +
      scale_fill_gradientn(colors=col2Cols, breaks=seq(0, 1, 0.2)) +
      labs(x=NULL, y=NULL, title=title) + 
      guides(fill = guide_colorbar(barheight = 10)) +
      theme(aspect.ratio=1, 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=13),
        axis.text.y = element_text(size=13),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank())
}

# Hierachical clustering plot
hierarchicalPlot <- function(MSnSetObj, sampleColours=NULL, colourBy="SampleGroup", horizontal=TRUE,
                             title=""){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!colourBy%in%colnames(pData(MSnSetObj))){ stop("colourBy must a column names in the pData of the MSnSetObj")}
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  
  dendro.dat <- t(log2xplus1(exprs(MSnSetObj))) %>% 
    dist(method = "euclidean") %>%
    hclust() %>% 
    as.dendrogram() %>% 
    dendro_data()
  labelDat <- dendro.dat$labels %>% 
    mutate(SampleName=as.character(label)) %>% 
    left_join(pData(MSnSetObj), "SampleName") %>% 
    mutate_at(vars(colourBy), funs(as.factor))
  axisBreaks <- pretty(dendro.dat$segments$yend)[-1] %>% head(-1)
  
  if(horizontal){ hj=0; ny=1; ang=0 }
  if(!horizontal){ hj=1; ny=-1; ang=90}
  hcPlot <- ggplot(dendro.dat$segment) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(data = labelDat, 
               aes_string(x = "x", y = "y", label = "SampleName", colour=colourBy), 
               hjust=hj, nudge_y=ny, angle=ang) +
    guides(colour=FALSE) +
    scale_fill_manual(values = sampleColours, breaks=names(sampleColours)) +
    labs(x=NULL, y="Distance", title=title)
  if(horizontal){
    hcPlot <- hcPlot +
      scale_y_reverse(expand = c(0.2, 0), breaks=axisBreaks) +
      coord_flip() +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank())
  }else{ 
    hcPlot <- hcPlot + 
      scale_y_continuous(expand = c(0.2, 0), breaks=axisBreaks) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank())
  }
   return(hcPlot)     
}


# Mean variance plot
plotMeanVar <- function(MSnSetObj, title=""){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  
  intensities <- log(exprs(MSnSetObj)+0.001)
  mvDat <- data.frame(
    Mean=rowMeans(intensities),
    Variance=apply(intensities, 1, var)
  )
  
  ssDat <- smooth.spline(x=mvDat$Mean, y=mvDat$Variance, spar=1) %$% 
    data.frame(x=x, y=y)
  
  ggplot(mvDat, aes(x=Mean, y=Variance)) + 
    geom_point(size=0.5, alpha=0.6, colour="darkblue") + 
    geom_line(data=ssDat, aes(x=x, y=y), colour="red", size=0.5) +
    theme_bw() +
    labs(title=title) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Protein coverage plot
coveragePlot <- function(MSnSetObj, ProteinID, ProteinName, fastaFile, myCol="brown"){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!is(ProteinID,"character")){ stop('ProteinID has to be of class character..') }
  if(!is(ProteinName,"character")){ stop('ProteinName has to be of class character..') }
  if(!is(fastaFile,"character")){ stop('fastaFile has to be of class character..') }
  if(!is(myCol,"character")){ stop('myCol has to be of class character..') }
  if(!"Sequences"%in%colnames(fData(MSnSetObj))){ 
    stop('MSnSetObj feature data must include a column of peptide sequences labelled "Sequences"..') 
    }
  
  ## read protein sequence from fastafile
  Protein_seq <- readAAStringSet(fastaFile)
  
  ## extract peptide sequence from MsnSet object and match peptide sequence with protein sequence and store co-ordinates
  getPosition <- function(peptideSeq, ProteinSeq=Protein_seq){
    vmatchPattern(peptideSeq, ProteinSeq) %>% 
          as.data.frame() %>% 
          select(start, end) %>% 
          return()
  }
  features <- fData(MSnSetObj) %>% 
      filter(Accessions==ProteinID) %>% 
      use_series(Sequence) %>% 
      gsub("^\\[.\\]\\.([A-Z]+)\\.\\[.\\]$", "\\1", .) %>% 
      map_dfr(getPosition)
  
  # get percent coverage 
  protWidth <- width(Protein_seq)
  coverage <- GRanges("feature",IRanges(features$start,features$end)) %>% reduce() %>% width() %>% sum()
  Perct <- round(coverage/protWidth*100,2)
  SubTitle <- paste0("Number of Unique Peptides: ", nrow(features), "\n% Coverage: ", Perct)
  
  # set the tick positions for the plot
  nTicks <- min(c(7, ceiling(protWidth/50) + 1))
  brkTicks <- round(seq(0,protWidth, length.out = nTicks), 0)
  
  features %>% 
    distinct() %>% 
    ggplot() +
    geom_rect(aes(xmin=start-1, xmax=end, ymin=0, ymax=10), fill=myCol) +
    geom_rect(xmin=0, xmax=protWidth, ymin=0, ymax=10, colour="black", fill=NA, size=0.15) +
    labs(title=ProteinName, subtitle=SubTitle) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_rect(fill = "white"), 
          panel.border = element_blank(),
          plot.subtitle = element_text(hjust = 0.5),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits=c(0,protWidth), breaks=brkTicks) +
    scale_y_continuous(limits=c(0,10), breaks=c(0,10), expand = c(0, 0))
}

# Intensity boxplots
intensityBoxplot <- function(MSnSetObj, title="", sampleColours=NULL, colourBy="SampleGroup"){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!colourBy%in%colnames(pData(MSnSetObj))){ stop("colourBy must a column names in the pData of the MSnSetObj")}
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  
  exprs(MSnSetObj) %>% 
    as.data.frame() %>% 
    gather("SampleName", "Intensity", everything()) %>% 
    mutate(logInt=log2(Intensity)) %>% 
    filter(is.finite(logInt)) %>% 
    left_join(pData(MSnSetObj), "SampleName") %>% 
    mutate_at(vars(colourBy), funs(as.factor)) %>% 
    ggplot() +
    geom_boxplot(aes_string(x="SampleName", y="logInt", fill=colourBy), alpha=0.6) +
      scale_fill_manual(values=sampleColours, breaks=names(sampleColours)) +
      labs(x="Sample", y="log2(Intensity)", title=title)  +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title=element_text(hjust=0.5)) +
      guides(fill=FALSE)
}

# Relative log intensity plot
rliPlot <- function(MSnSetObj, title="", sampleColours=NULL, colourBy="SampleGroup", omitIgG=TRUE){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!colourBy%in%colnames(pData(MSnSetObj))){ stop("colourBy must a column names in the pData of the MSnSetObj")}
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  
  # Remove IgG samples is requested
  if(omitIgG){ MSnSetObj <- MSnSetObj[,toupper(MSnSetObj$SampleGroup)!="IGG"] }
  exprs(MSnSetObj) %>% 
    as.data.frame() %>% 
    rownames_to_column("RowID") %>% 
    gather("SampleName", "Intensity", -RowID) %>% 
    mutate(logInt=log2(Intensity)) %>% 
    filter(is.finite(logInt)) %>% 
    group_by(RowID) %>% 
    mutate(medianLogInt=median(logInt)) %>% 
    ungroup() %>% 
    mutate(RLI=logInt-medianLogInt) %>% 
    left_join(pData(MSnSetObj), "SampleName") %>% 
    mutate_at(vars(colourBy), funs(as.factor)) %>%
    ggplot(aes_string(x="SampleName", y="RLI", fill=colourBy)) +
      geom_boxplot(alpha=0.6) +
      scale_fill_manual(values=sampleColours, breaks=names(sampleColours)) +
      labs(x="Sample", y="RLI", title=title)  +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title=element_text(hjust=0.5))  +
      guides(fill=FALSE)
}
