##### This scripts consists of all the plotting functions #####

## Assigns colours to samples in groups
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
      coloursF <- brewer.pal(8, "Dark2") %>% grDevices::colorRampPalette()
      sampleColours <- setNames(coloursF(len), colourGroups)
  }
  return(sampleColours)
}

# Intensity distribution plot

# intensities is a data frame containing columns for each sample
# sampleColours is a named vector that maps samples to colours
intensityPlot <- function(MSnSetObj, sampleColours=NULL, colourBy="SampleGroup", 
                          transform=TRUE, xlab="log2(intensity)",
                          trFunc=log2xplus1){
  if (!is(MSnSetObj, "MSnSet")){ stop("MSnSetObj has to be of class MSnSet..") }
  if(!transform){ trFunc <- c }
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  exprs(MSnSetObj) %>% 
    as.data.frame() %>% 
    gather(SampleName, Intensity) %>%
    left_join(pData(MSnSetObj), "SampleName") %>% 
    mutate(Intensity=trFunc(Intensity)) %>% 
    filter(!is.na(Intensity)) %>%
    ggplot(aes_string(x="Intensity", group="SampleName", colour=colourBy)) + 
    stat_density(geom="line", position="identity") +
    scale_colour_manual(values=sampleColours, breaks=names(sampleColours)) +
    labs(x=xlab) +
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
      ) %>% 
    return()
}

# peptide intensity plot

# intensities is a data frame containing peptide intensities with columns for each sample
# summarizedIntensities is a data frame containing summarized protein-level intensities
# protein is the protein for which intensities will be plotted
# samples is a list of samples to use in the plot

peptideIntensityPlot <- function(MSnSetObj, combinedIntensities=NULL, protein,
                                 selectedSequence=NULL, selectedModifications=NULL){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!is.null(combinedIntensities) && !is(combinedIntensities,"MSnSet")){
      stop('combinedIntensities has to either NULL or object of class MSnSet..')
  }
  
  intensities <- exprs(MSnSetObj) %>% 
    as.data.frame() %>% 
    rownames_to_column("PeptideID") %>%
    gather(SampleName, Intensity, -PeptideID) %>% 
    left_join(fData(MSnSetObj) %>% rownames_to_column("PeptideID"), "PeptideID") %>% 
    mutate(logIntensity=log2xplus1(Intensity)) %>% 
    filter(Accessions == protein)
  
  seqIntensities <-  filter(intensities, Sequences %in% selectedSequence)
  if (!is.null(selectedModifications)){ 
    seqIntensities %<>% filter(Modifications %in% selectedModifications)
  }
  
  summIntensities <- data.frame(SampleName=vector(), logIntensity=vector(), Protein=vector())
  if (!is.null(combinedIntensities)){
    summIntensities <- exprs(combinedIntensities) %>% 
      as.data.frame() %>% 
      rownames_to_column("Protein") %>%
      gather(SampleName, Intensity, -Protein) %>% 
      left_join(fData(combinedIntensities), "Protein") %>% 
      mutate(logIntensity=log2xplus1(Intensity)) %>% 
      filter(Protein == protein)
  }
  
  if (nrow(intensities) == 0) {
    warning("No peptides were found for ", protein)
    return(NULL)
  }
  
  intPepPlot <- ggplot(intensities) +
    geom_line(aes(x=SampleName, y=logIntensity, group=PeptideID, colour=Sequences),
              size=0.6, alpha=0.5, linetype=2) + 
    geom_point(aes(x=SampleName, y=logIntensity, fill=Sequences), shape=21, size=2) + 
    geom_line(data=seqIntensities, aes(x=SampleName, y=logIntensity, group=PeptideID),
              colour="#6666FF", size=1.2, linetype=6)  + 
    geom_point(data=seqIntensities, aes(x=SampleName, y=logIntensity),
               fill="#0000FF", shape=21, size=2.5) + 
    geom_line(data=summIntensities, aes(x=SampleName, y=logIntensity, group=Protein),
              colour="#AAAAAA", size=1.2, linetype=6) +
    geom_point(data=summIntensities, aes(x=SampleName, y=logIntensity),
               fill="#888888", shape=21, size=2.5) + 
    labs(y="log2(Intensity)") +
    theme_bw() +
    theme(plot.title=element_text(size=14,hjust=0.5),
          axis.text.x=element_text(size=11, angle=45, hjust=1),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.position="none") %>% 
    return()
}


# PCA plot
pcaPlot <- function(MSnSetObj, omitIgG=FALSE, sampleColours=NULL, transFunc=log2xplus1, transform=TRUE, 
                    colourBy="SampleGroup", labelColumn=NULL, labelsize=4, x.nudge=4, x.PC=1){   
    
  if(!is(MSnSetObj,"MSnSet")){ stop('data has to be of class MSnSet..') }
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  
  # Remove IgG samples is requested
  if(omitIgG){ MSnSetObj <- MSnSetObj[,toupper(MSnSetObj$SampleGroup)!="IGG"] }
  if(!transform){ transFunc <- as.data.frame }
  intensities <- exprs(MSnSetObj) %>% as.data.frame() %>% na.omit() %>% transFunc()
  if (nrow(intensities) == 0){ return(NULL) }
  pca <- intensities %>% t() %>% prcomp()
  pcaVariance <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100)
  plotDat <- as.data.frame(pca$x) %>% 
    rownames_to_column("SampleName") %>% 
    left_join(pData(MSnSetObj), "SampleName")
  xPC <- paste0("PC", x.PC)
  yPC <- paste0("PC", x.PC+1)
  pcaPt <- ggplot(plotDat) + 
    geom_point(aes_string(x=xPC, y=yPC, fill=colourBy), pch=21, colour="black", size=labelsize) +
    geom_text(aes_string(x=xPC, y=yPC, label=labelColumn), hjust=0, size=labelsize, nudge_x=x.nudge) +
    scale_fill_manual(values=sampleColours, breaks=names(sampleColours)) +
    labs(x=paste0(xPC, ", ", pcaVariance[x.PC], "% variance"),
         y=paste0(yPC, ", ", pcaVariance[x.PC+1], "% variance"),
         fill=NULL) + 
    theme_bw() + 
    theme(text=element_text(size=14),
          plot.title=element_text(size=14,hjust=0.5),
          aspect.ratio=1) 
  return(pcaPt)
}



# MA or volcano plot
maVolPlot <- function(diffstats, contrast, controlGroup = NULL, ann,
                      selectedGenes=NULL, fdrCutOff=0.05,
                      lfcCutOff=1.5, controlLfcCutOff=1, plotType="MA"){
  # For plotting we will assign the proteins to one of 7 groups:
  # A - selected (user specified in `selectedGenes`) & significant
  # B - selected (user specified in `selectedGenes`) & non-significant
  # C - significant, specific
  # D - significant, non-specific
  # E - non-significant, large log2FC, specific
  # F - non-significant, large log2FC, non-specific
  # G - specific
  # H - non-specific
  altLevels <- function(x, newLev){ levels(x) <- newLev; return(x) }
  
  daResTab <- suppressMessages(getContrastResults(diffstats=diffstats, contrast=contrast, 
                                                  controlGroup = controlGroup, ann = ann)) %>% 
    bind_rows(data.frame(controlLogFoldChange=vector())) %>% 
    mutate(controlLogFoldChange=ifelse(is.na(controlLogFoldChange), Inf, controlLogFoldChange)) %>% 
    mutate(group=ifelse(controlLogFoldChange>=controlLfcCutOff, 7, 8)) %>% 
    mutate(group=ifelse(group==8 & '%in%'(adj.P.Val<=fdrCutOff, T), 4, group)) %>%
    mutate(group=ifelse(group==7 & '%in%'(adj.P.Val<=fdrCutOff, T), 3, group)) %>%
    mutate(group=ifelse(group==8 & abs(log2FC)>=lfcCutOff, 6, group)) %>%
    mutate(group=ifelse(group==7 & abs(log2FC)>=lfcCutOff, 5, group)) %>%
    mutate(group=ifelse(GeneSymbol%in%selectedGenes, 2, group)) %>% 
    mutate(group=ifelse(group==2 & '%in%'(adj.P.Val<=fdrCutOff, T), 1, group)) %>% 
    arrange(desc(group)) %>% 
    mutate(group=factor(LETTERS[group], levels=LETTERS[1:8])) %>% 
    mutate(group=altLevels(group, c("User Selected & Significant", "User Selected", 
                                    "Specific & Significant", "Non-Specific & Significant", 
                                    LETTERS[5:8]))) %>% 
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
      geom_text(data=subset(daResTab, group%in%levels(group)[1:2]), aes(label=GeneSymbol), hjust=0, 
              size=3.5, nudge_x = xNudge) +
      scale_colour_manual(values=rep(c("black", "gray50"), each=4), drop=FALSE) +
      scale_size_manual(values=rep(c(1.8, 1.5, 1.2, 0.9), each=2), drop=FALSE) +
      scale_shape_manual(values=rep(c(21, 20), each=4), drop=FALSE) +
      scale_alpha_manual(values=rep(c(1, 1, 0.8, 0.6), each=2), drop=FALSE)  +
      scale_fill_manual(values=c("cyan", "purple", "red", "orange"), limits=levels(daResTab$group)[1:4], 
                        breaks=levels(droplevels(daResTab$group)), name="") +
      labs(x=xLab, y=yLab) + 
      theme_bw() +
      theme(
          text=element_text(size=16),
          plot.title=element_text(size=14,hjust=0.5)
      ) +
      guides(colour="none", size="none", shape="none", alpha="none",
             fill = guide_legend(override.aes = list(shape = 21))) %>% 
    return()
}

# MA plot wrapper
maPlot <- function(diffstats, contrast, controlGroup = NULL, ann,
                   selectedGenes=NULL, fdrCutOff=0.05,
                   lfcCutOff=1.5, controlLfcCutOff=1){
    
    maVolPlot(diffstats=diffstats, contrast=contrast, controlGroup = controlGroup, 
              ann = ann,  selectedGenes=selectedGenes, fdrCutOff=fdrCutOff, 
              lfcCutOff=lfcCutOff, controlLfcCutOff=controlLfcCutOff, plotType="MA") %>% 
        return()
}

# Volcano plot wrapper
volPlot <- function(diffstats, contrast, controlGroup = NULL, ann, 
                    fdrCutOff=0.05, lfcCutOff=1.5, controlLfcCutOff=1, selectedGenes=NULL){
    maVolPlot(diffstats=diffstats, contrast=contrast, controlGroup = controlGroup, 
              ann = ann,  selectedGenes=selectedGenes, fdrCutOff=fdrCutOff, 
              lfcCutOff=lfcCutOff, controlLfcCutOff=controlLfcCutOff, plotType="Volcano") %>% 
        return()
}

# correlation plot
corrPlot <- function(MSnSetObj, method="shade", title="Correlation plot"){
  if(!is(MSnSetObj,"MSnSet")){ stop('data has to be of class MSnSet..') }
  
  col2Cols <- c("#FFFFFF", "#DCF0F3", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")
  corPlotObj <- cor(exprs(MSnSetObj)) %>% 
    as.data.frame() %>% 
    rownames_to_column("X") %>% 
    gather("Y", "Cor", -1) %>%
    ggplot(aes(x=X, y=Y, fill=Cor)) +
    geom_tile()
  
  corPlotObj +
    geom_tile(col="grey") +
    scale_fill_gradientn(colors=col2Cols, breaks=seq(0, 1, 0.2)) +
    labs(x=NULL, y=NULL) + 
    guides(fill = guide_colorbar(barheight = 10)) +
    theme(aspect.ratio=1, 
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=13),
      axis.text.y = element_text(size=13),
      panel.background = element_blank()) %>% 
    return()
}

# hierachical clustering plot
hierarchicalPlot <- function(MSnSetObj, sampleColours=NULL, colourBy="SampleGroup"){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  
  dendro.dat <- t(log2xplus1(exprs(MSnSetObj))) %>% 
    dist(method = "euclidean") %>%
    hclust() %>% 
    as.dendrogram() %>% 
    dendro_data()
  labelDat <- dendro.dat$labels %>% 
    rename(SampleName=label) %>% 
    left_join(pData(MSnSetObj), "SampleName")
  axisBreaks <- pretty(dendro.dat$segments$yend)[-1] %>% head(-1)
  
  ggplot(dendro.dat$segment) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_label(data = labelDat, 
               aes_string(x = "x", y = "y", label = "SampleName", fill=colourBy), 
               vjust = 0.5, hjust=0, nudge_y=1, alpha=0.4) + 
    scale_y_reverse(expand = c(0.2, 0), breaks=axisBreaks) +
    guides(fill=F) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank()) +
    scale_fill_manual(values = sampleColours, breaks=names(sampleColours)) +
    coord_flip() + 
    labs(x=NULL, y="Distance")
}


# plot mean var
plotMeanVar <- function(MSnSetObj){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  
  intensities <- log(exprs(MSnSetObj)+0.001)
  mvDat <- data.frame(
    Mean=rowMeans(intensities),
    Variance=apply(intensities, 1, var)
  )
  
  ssDat <- smooth.spline(x=mvDat$Mean, y=mvDat$Variance, spar=1) %$% 
    data.frame(x=x, y=y)
  
  ggplot(mean.var_res, aes(x=Mean, y=Variance)) + 
    geom_point(size=0.5, alpha=0.6, colour="darkblue") + 
    geom_line(data=ssdat, aes(x=x, y=y), colour="red", size=0.5) +
    theme_bw()
}

### Protein coverage plot
coveragePlot <- function(MSnSetObj, ProteinID, fastaFile, myCol="brown"){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!is(ProteinID,"character")){ stop('ProteinID has to be of class character..') }
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
          dplyr::select(start, end) %>% 
          return()
  }
  features <- fData(MSnSetObj) %>% 
      filter(Accessions==ProteinID) %>% 
      use_series(Sequence) %>% 
      gsub("^\\[.\\]\\.([A-Z]+)\\.\\[.\\]$", "\\1", .) %>% 
      purrr::map_dfr(getPosition)
  
  # get percent coverage 
  protWidth <- width(Protein_seq)
  coverage <- GRanges("feature",IRanges(features$start,features$end)) %>% IRanges::reduce() %>% width() %>% sum()
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
    labs(title=ProteinID, subtitle=SubTitle) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_rect(fill = "white"), 
          panel.border = element_blank(),
          plot.subtitle = element_text(hjust = 0.5),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits=c(0,protWidth), breaks=brkTicks) +
    scale_y_continuous(limits=c(0,10), breaks=c(0,10), expand = c(0, 0)) %>% 
    return()
}

# Intensity boxplots
intensityBoxplot <- function(MSnSetObj, sampleColours=NULL, colourBy="SampleGroup"){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  
  exprs(MSnSetObj) %>% 
    as.data.frame() %>% 
    gather("SampleName", "Intensity", everything()) %>% 
    mutate(logInt=log2(Intensity)) %>% 
    filter(is.finite(logInt)) %>% 
    left_join(pData(MSnSetObj), "SampleName") %>% 
    ggplot() +
    geom_boxplot(aes_string(x="SampleName", y="logInt", fill=colourBy), alpha=0.6) +
      scale_fill_manual(values=sampleColours, breaks=names(sampleColours)) +
      labs(x="Sample", y="log2(Intensity)")  +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      guides(fill=F) %>% 
    return()
}

# relative log intensity plot
rliPlot <- function(MSnSetObj, sampleColours=NULL, colourBy="SampleGroup"){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(is.null(sampleColours)){ sampleColours <- assignColours(MSnSetObj, colourBy=colourBy) }
  
  exprs(MSnSetObj) %>% 
    as.data.frame() %>% 
    rownames_to_column("RowID") %>% 
    select(-starts_with("IgG")) %>% 
    gather("SampleName", "Intensity", -RowID) %>% 
    mutate(logInt=log2(Intensity)) %>% 
    filter(is.finite(logInt)) %>% 
    group_by(RowID) %>% 
    mutate(medianLogInt=median(logInt)) %>% 
    ungroup() %>% 
    mutate(RLI=logInt-medianLogInt) %>% 
    left_join(pData(MSnSetObj), "SampleName") %>% 
    ggplot() +
      geom_boxplot(aes_string(x="SampleName", y="RLI", fill=colourBy), alpha=0.6) +
      scale_fill_manual(values=sampleColours, breaks=names(sampleColours)) +
      labs(x="Sample", y="RLI")  +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
      guides(fill=F) %>% 
    return()
}
