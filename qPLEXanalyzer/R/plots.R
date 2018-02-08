##### This scripts consists of all the plotting functions #####

## Assigns colours to samples in groups
## check this for selecting colors http://colorbrewer2.org/

assignColours <- function(data, groupColourSchemes)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  phenodata <- pData(data)
  groups <- levels(phenodata$Label)
  sampleColours <- character(0)
  # TODO add check for number of groups exceeding the number of colour schemes
  for (i in 1:length(groups))
  {
    samples <- phenodata %>% filter(Label == groups[i]) %>% select(Experiment) %>% unlist(use.names = FALSE)
    # TODO add check for number of samples exceeding the maximum number of colours for the selected colour scheme
    colours <- c(brewer.pal(length(samples) + 2, groupColourSchemes[i]))[2:(length(samples) + 1)]
    names(colours) <- samples
    sampleColours <- c(sampleColours, colours)
  }
  sampleColours
}


# Intensity distribution plot

# intensities is a data frame containing columns for each sample
# sampleColours is a named vector that maps samples to colours
intensityPlot <- function(data, sampleColours, log2Transform=TRUE, title = "", 
                                      xlab = "log2(intensity)", minIntensity = NA, maxIntensity = NA)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  samples <- as.character(data$Experiment)  
  colours <- sampleColours[samples] %>% as.character
  intensities <- as.data.frame(exprs(data))
  if(log2Transform)
  {
    log2xplus1 <- function(x) { log2(x + 1) }
    intensities <- log2xplus1(intensities) 
  }
  intensities <- intensities %>%
    select(one_of(samples)) %>%
    gather(sample, intensity) %>%
    filter(complete.cases(.))

  intensities$sample <- factor(intensities$sample, levels = samples)

  plot <- ggplot(intensities, aes(x = intensity, colour = sample))
  plot <- plot + stat_density(geom = "line", position = "identity")
  plot <- plot + scale_colour_manual(values = colours, drop = FALSE)
  plot <- plot + ggtitle(title)
  plot <- plot + xlab(xlab)
  plot <- plot + theme_bw()
  plot <- plot + theme(
    plot.title = element_text(size = 14,hjust=0.5),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 11),
    legend.justification = c(1,1),
    legend.position = c(1,1),
    legend.background = element_rect(fill = "transparent")
  )
  if (!is.na(minIntensity) && !is.na(maxIntensity))
    plot <- plot + xlim(minIntensity, maxIntensity)

  return(plot)
}


# peptide intensity plot

# intensities is a data frame containing peptide intensities with columns for each sample
# summarizedIntensities is a data frame containing summarized protein-level intensities
# protein is the protein for which intensities will be plotted
# samples is a list of samples to use in the plot

peptideIntensityPlot <- function(peptideIntensities, combinedIntensities=NULL, protein,
                                 title = "", ylab = "log2(intensity)",
                                 minIntensity = NA, maxIntensity = NA,
                                 selectedSequence = NULL, selectedModifications = NULL)
{
  if(!is(peptideIntensities,"MSnSet"))
    stop('peptideIntensities has to be of class MSnSet..')
  if(!is.null(combinedIntensities) && !is(combinedIntensities,"MSnSet"))
    stop('combinedIntensities has to either NULL or object of class MSnSet..')
  log2xplus1 <- function(x) { log2(x + 1) }
  intensities <- as.data.frame(exprs(peptideIntensities))
  intensities <- log2xplus1(intensities) 
  intensities <- cbind(fData(peptideIntensities),intensities)
  samples <- as.character(pData(peptideIntensities)$Experiment)
  
  
  intensities <- intensities %>%
    filter(Master.Protein.Accessions == protein)


  if (nrow(intensities) == 0) {
    intensities[1,2] <- ""
    intensities[1,1] <- "No peptides"
  }

  intensities <- intensities %>%
    mutate(id = 1:nrow(intensities)) %>%
    gather(sample, intensity, one_of(samples))

  intensities$sample <- factor(intensities$sample, levels = samples)

  plot <- ggplot()
  plot <- plot + ggtitle(title)
  plot <- plot + ylab(ylab)
  plot <- plot + geom_line(data = intensities, aes(x = sample, y = intensity, group = id, colour = Annotated.Sequence), size = 0.6)
  plot <- plot + geom_point(data = intensities, aes(x = sample, y = intensity, colour = Annotated.Sequence), size = 2)
  if (!is.na(minIntensity) && !is.na(maxIntensity))
  {
    plot <- plot + coord_cartesian(ylim = c(minIntensity, maxIntensity))
    plot <- plot + scale_y_continuous(breaks = seq(minIntensity, maxIntensity, 2))
  }
  plot <- plot + theme_bw()
  plot <- plot + theme(
    plot.title = element_text(size = 14,hjust=0.5),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

  if (!is.null(combinedIntensities))
  {
    summarizedIntensities <- as.data.frame(exprs(combinedIntensities))
    summarizedIntensities <- log2xplus1(summarizedIntensities) 
    summarizedIntensities <- cbind(fData(combinedIntensities),summarizedIntensities)
    samples <- as.character(pData(combinedIntensities)$Experiment)
  
    summarizedIntensities <- summarizedIntensities %>%
      filter(Protein == protein)
    
    if (sum(summarizedIntensities$Count) > 0)
    {
      
      summarizedIntensities <- summarizedIntensities %>%
        select(one_of(samples)) %>%
        gather(sample, intensity)

      summarizedIntensities$sample <- factor(summarizedIntensities$sample, levels = samples)

      plot <- plot + geom_point(data = summarizedIntensities, aes(x = sample, y = intensity), colour = "grey50", size = 2.5)
      plot <- plot + geom_line(data = summarizedIntensities, aes(x = sample, y = intensity, group = 1), colour = "grey50", size = 1.5)
    }
  }

  if (!is.null(selectedSequence))
  {
    selectedIntensities <- intensities %>% filter(Annotated.Sequence %in% selectedSequence)
    if (!is.null(selectedModifications))
      selectedIntensities <- selectedIntensities %>% filter(Modifications %in% selectedModifications)
    if (nrow(selectedIntensities) > 0)
    {
      plot <- plot + geom_point(data = selectedIntensities, aes(x = sample, y = intensity), colour = "blue", size = 2.5)
      plot <- plot + geom_line(data = selectedIntensities, aes(x = sample, y = intensity, group = id), colour = "blue", size = 1.5)
    }
  }

  return (plot)
}


# PCA plot

# groupColours is a named vector that maps groups to colours

pcaPlot <- function(data, groupColours, title = "", labels = NULL, legend = TRUE, labelsize=2.5, logTransform = TRUE)
{
  samples <- as.character(data$Experiment)  
  intensities <- as.data.frame(exprs(data))
  groups <- data$Label
  colours <- groupColours %>% as.character
  intensities <- intensities %>%
    select(one_of(samples)) %>%
    filter(complete.cases(.))
  if (nrow(intensities) == 0) return(ggplot() + theme_bw())
  if (logTransform)
  {
    intensities <- log(intensities) 
  }
  pca <- intensities %>% t %>% prcomp
  pcaVariance <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100)
  pca <- as.data.frame(pca$x)
  pca$Group <- groups
  plot <- ggplot(pca)
  plot <- plot + ggtitle(title)
  
  if (!is.null(labels))
    plot <- plot + geom_text(aes(x = PC1, y = PC2, label = paste("  ", labels)), hjust = 0, size = labelsize)
  plot <- plot + geom_point(aes(x = PC1, y = PC2, colour = Group), size = labelsize)
  plot <- plot + scale_colour_manual(values = colours)
  plot <- plot + scale_x_continuous(expand = c(0.2, 0))
  plot <- plot + scale_y_continuous(expand = c(0.2, 0))
  plot <- plot + xlab(paste("PC1, ", pcaVariance[1], "% variance", sep = ""))
  plot <- plot + ylab(paste("PC2, ", pcaVariance[2], "% variance", sep = ""))
  plot <- plot + theme_bw()
  plot <- plot + theme(
    text = element_text(size = 14),
    plot.title = element_text(size = 14,hjust=0.5))
  if (legend)
    plot <- plot + theme(
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = "right"
    )
  else
    plot <- plot + theme(legend.position = "none")
  return(plot)
}


# MA plot
maPlot <- function(differentialExpressionResults, selectedGenes = NULL,
                   xlab = "average log2(intensity)",
                   ylab = "log2FC",
                   significanceLevel = 0.05, title = "",
                   minLogFoldChangeForLabelling = 1.5,
                   controlLogFoldChangeThreshold = -Inf,
                   pointSize = 1.5)
{
#  cat("MA plot\n")

  if (!"controlLogFoldChange" %in% colnames(differentialExpressionResults))
    differentialExpressionResults$controlLogFoldChange <- Inf

#  differentialExpressionResults

  differentialExpressionResults$group <-
    ifelse(differentialExpressionResults$GeneSymbol %in% selectedGenes, 0,
    ifelse(!is.na(differentialExpressionResults$adj.P.Val) & differentialExpressionResults$adj.P.Val <= significanceLevel & differentialExpressionResults$controlLogFoldChange >= controlLogFoldChangeThreshold, 1,
    ifelse(!is.na(differentialExpressionResults$adj.P.Val) & differentialExpressionResults$adj.P.Val <= significanceLevel, 2,
    ifelse(abs(differentialExpressionResults$log2FC) >= minLogFoldChangeForLabelling & differentialExpressionResults$controlLogFoldChange >= controlLogFoldChangeThreshold, 3,
    ifelse(abs(differentialExpressionResults$log2FC) >= minLogFoldChangeForLabelling, 4,
    ifelse(differentialExpressionResults$controlLogFoldChange >= controlLogFoldChangeThreshold, 5, 6))))))

  # 0 - selected
  # 1 - significant, specific
  # 2 - significant, non-specific
  # 3 - non-significant, large log2FC, specific
  # 4 - non-significant, large log2FC, non-specific
  # 5 - specific
  # 6 - non-specific

  differentialExpressionResults <- differentialExpressionResults %>% arrange(desc(group))

  differentialExpressionResults$group <- factor(differentialExpressionResults$group, 0:6)

  differentialExpressionResults <- differentialExpressionResults %>%
    mutate(label = paste(" ", GeneSymbol))

  plot <- ggplot(differentialExpressionResults, aes(x = AvgIntensity, y = log2FC, colour = group, size = group, shape = group, alpha = controlLogFoldChange))
  plot <- plot + geom_hline(yintercept = 0, color="gray50", size = 0.5)
  plot <- plot + geom_point()
  plot <- plot + geom_text(data = subset(differentialExpressionResults, group != 0 & (adj.P.Val <= significanceLevel | abs(log2FC) >= minLogFoldChangeForLabelling)), aes(label = label), hjust = 0, size = 3.5)
  plot <- plot + geom_point(data = subset(differentialExpressionResults, group == 0), alpha = 1.0)
  plot <- plot + geom_text(data = subset(differentialExpressionResults, group == 0), aes(label = label), hjust = 0, size = 5.5, alpha = 1.0)
  plot <- plot + scale_colour_manual(values = c("blue", "deeppink3", "deeppink2", "gray50", "gray50", "gray50", "gray50"), drop = FALSE)
  plot <- plot + scale_size_manual(values = c(1.2 * pointSize, pointSize, pointSize, 0.8 * pointSize, 0.8 * pointSize, 0.6 * pointSize, 0.6 * pointSize), drop = FALSE)
  plot <- plot + scale_shape_manual(values = c(16, 16, 1, 16, 1, 16, 1), drop = FALSE)
  plot <- plot + scale_alpha(range = c(0.2, 1.0))
  plot <- plot + xlab(xlab)
  plot <- plot + ylab(ylab)
  plot <- plot + ggtitle(title)
  plot <- plot + theme_bw()
  plot <- plot + theme(
    text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 14,hjust=0.5)
  )
  return (plot)
}


# correlation plot

corrPlot <- function(data, method="shade", title="Correlation plot")
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  cor_mat_bn <- cor(exprs(data))
  corrplot(cor_mat_bn,method=method, title=title, mar=c(5.1,4.1,4.1,2.1))
}

hierarchicalPlot <- function(data,label_color,branchlength=20,title)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  clust_mat <- t(log2(exprs(data)))
  dist_mat <- dist(clust_mat, method = "euclidean")
  hclust.obj <- hclust(dist_mat)
  ColorDendrogram(hc=hclust.obj, y=label_color,labels=hclust.obj$labels,
                  branchlength=branchlength, main=title, xlab="Samples")
}


plotMeanVar <- function(data,title,color,ylim=NULL,xlim=NULL)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  if(!is.character(title))
    stop('title has to be group name of class character ..')
  if(!is.character(color))
    stop('color has to be of class character ..')
  if(!is.null(ylim) && !is.numeric(ylim))
    stop('ylim has to be either NULL or of class numeric ..')
  if(!is.null(xlim) && !is.numeric(xlim))
    stop('xlim has to be either NULL or of class numeric ..')
  
  intensities <- log(as.data.frame(exprs(data)+0.001))
  mean.var <- function(x)
  {
    mean.variance <- cbind(mean=apply(x,1,mean),var=apply(x,1,var))
    return(as.data.frame(mean.variance))
  }
  
  mean.var_res <- mean.var(intensities)
  smoothingSpline = smooth.spline(x=mean.var_res$mean, y=mean.var_res$var, spar=1)
  if(is.null(ylim) && is.null(xlim))
  {
    smoothScatter(x=mean.var_res$mean, y=mean.var_res$var,xlab="Mean",ylab="Variance",main=title)
  }
  else 
    smoothScatter(x=mean.var_res$mean, y=mean.var_res$var,ylim=ylim,xlim=xlim,xlab="Mean",ylab="Variance",main=title)
  lines(smoothingSpline,col=color)
}

### Protein coverage plot

coveragePlot <- function(data,ProteinID,name,fastaFile, col="brown")
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  if(!is(ProteinID,"character"))
    stop('ProteinID has to be of class character..')
  if(!is(name,"character"))
    stop('name has to be of class character..')
  if(!is(fastaFile,"character"))
    stop('fastaFile has to be of class character..')
  if(!is(col,"character"))
    stop('col has to be of class character..')
  name <- as.character(name)
  # extract peptide sequence from MsnSet object
  ind <- which(fData(data)$Master.Protein.Accessions==ProteinID)
  ref_seq <- as.character(fData(data)$Annotated.Sequence[ind])
  peptide_seq <- sapply(ref_seq,function(x) unlist(strsplit(x,split="[.]"))[2])
  
  ## read protein sequence from fastafile
  Protein_seq <- readAAStringSet(fastaFile)
  
  ## match peptide sequence with protein sequence and store co-ordinates
  sind <- array()
  eind <- array()
  for (i in 1:length(peptide_seq))
  {
    sind[i] <- unlist(startIndex(vmatchPattern(peptide_seq[i],Protein_seq)))
    eind[i] <- unlist(endIndex(vmatchPattern(peptide_seq[i],Protein_seq)))
  }
  sind <- sind-1
  # list of names...
  names <- list(
    description = paste("Number of Unique Peptides",":",length(unique(peptide_seq))),
    name = name)
  
  ## create the feature data frame with all the information to plot 
  begin <- sind
  end <- eind
  col <- rep(col, length(sind))
  features <- data.frame(begin, end, col)
  features$col <- as.character(features$col)  
  
  # sort features in order of where they begin
  features <- features[order(features$begin),]
  features <- unique(features)
  gr <- GRanges("feature",IRanges(features$begin,features$end))
  gr <- union(gr,gr)
  Perct <- round((sum(end(gr)-start(gr)))/width(Protein_seq)*100,2)
  # add additonal peptide co-ordinate to encompass entire protein sequence for plotting purpose
  total <- data.frame(1,width(Protein_seq),"NA")
  colnames(total) <- colnames(features)
  features <- rbind(features,total)
  
  ## draw the diagram
  screen.width <- width(Protein_seq)
  screen.height <- 30  # this is a bit arbitary
  if(width(Protein_seq) <= 50)
    lengthOut <- 3 else if (width(Protein_seq) > 50 & width(Protein_seq) <= 100)
      lengthOut <- 4 else if(width(Protein_seq) > 100 & width(Protein_seq) <= 200)
        lengthOut <- 5 else 
          lengthOut <- 7
  
  
  xticks <- round(seq(1,width(Protein_seq),length.out = lengthOut))
  toadd <- round(xticks[2:(length(xticks)-1)], digits = -1) 
  final_ticks <- c(head(xticks,n=1),toadd,tail(xticks,n=1))
  #par(mar=c(10,4.1,1,2.1))
  plot(0,xlim=c(0,screen.width),ylim=c(0,100),type="n",axes = F,ann = F)
  axis(side = 1, at = final_ticks)
  
  # make the rectangles in a loop
  for (i in 1:length(features$begin) ) {
    rect(xleft   = features$begin[i],
         ytop    = 0,
         ybottom = 10,
         xright  = features$end[i],
         col = features$col[i],lwd=1,border=NA)
  }
  rect(xleft   = features$begin[i],
       ytop    = 0,
       ybottom = 10,
       xright  = features$end[i],
       col = features$col[i],lwd=1)
  
  # add text to the top of the illustration with the recommended name
  text(max(features$end)/2, screen.height+10, names$name, cex=1.5)
  
  # add information about number of peptides
  text(max(features$end)/2, screen.height-2, names$description, cex=1)
  
  # add information about percentage coverage of protein
  text(max(features$end)/2, screen.height-12 , paste("% Coverage:", Perct), cex=0.8)
}
