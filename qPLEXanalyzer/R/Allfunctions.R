######### This script contains all the functions required for data processing and analysis ###############


convertToMSnset <- function(data,metadata,indExpData,Sequences,Accessions,rmMissing=TRUE)
{
  if (!is.data.frame(data))
    stop("data has to be of class dataframe")
  if (!is.data.frame(metadata))
    stop("metadata has to be of class dataframe")
  if(!is.numeric(indExpData))
    stop('indExpData has to be of class numeric ..')
  if(!is.numeric(Sequences))
    stop('Sequences has to be of class numeric ..')
  if(!is.numeric(Accessions))
    stop('Accessions has to be of class numeric ..')
  columns <- c("SampleName","SampleGroup","BioRep","TechRep")
  if(length(which((columns %in% colnames(metadata))==FALSE)) > 0)
    stop('metadata must contain "SampleName","SampleGroup", "BioRep" and "TechRep" columns ..')
  if (!is.logical(rmMissing))
    stop("rmMissing has to be of class logical")
  colnames(data)[Sequences] <- "Sequences"
  colnames(data)[Accessions] <- "Accessions"
  if(rmMissing)
    data %<>% filter_at(vars(indExpData), all_vars(!is.na(.)))
  obj <- readMSnSet2(data,ecol=indExpData)
  rownames(metadata) <- as.character(metadata$SampleName)
  check <- all(rownames(metadata) %in% sampleNames(obj))
  if(check)
  {
    ind <- match(sampleNames(obj), rownames(metadata))
    metadata <- metadata[ind,]
  }
  else
    stop("Problem consistency between column names \n            in expression data and SampleName in metadata")
  pData(obj) <- metadata
  featureNames(obj) <- paste0("peptide_",featureNames(obj))
  return(obj)
}  

  

##### Summarization function ##########

# functions for summarizing intensity data

# Summarizes multiple peptide measurements for a protein.
# Assumes that there are columns for each of the samples specified and uses
# Protein column for grouping peptide-level measurements.
# Filters any rows with missing values.
# Typical summarization functions are sum, mean and median.
# For successful running of this function the annotation file must have four column "Protein","Gene","Description" and "GeneSymbol"
# In addition the PD output must have columns "Sequences" and "Accessions"

summarizeIntensities <- function(data, summarizationFunction, annotation)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  if(!is.data.frame(annotation))
    stop('annotation has to be of class data frame..')
  columns <- c("Sequences","Accessions")
  if(length(which((columns %in% colnames(fData(data)))==FALSE)) > 0)
    stop('This MSnSet is not a peptide dataset ..')
  Proteins <- as.character(fData(data)$Accessions)
  features <- fData(data)
  features <- as.data.frame(features[,c("Sequences","Accessions")], stringsAsFactors=FALSE)
  features <- unique(features)
  features$Sequences <- as.character(features$Sequences)
  features$Accessions <- as.character(features$Accessions)
  counts <- features %>% count(Accessions) %>% rename(Protein=Accessions,Count = n)
  intensities <- cbind.data.frame(exprs(data),Protein=Proteins)
  summIntensities <- intensities %>%
    group_by(Protein) %>%
    summarize_all(funs(summarizationFunction))
  summIntensities$Protein <- as.character(summIntensities$Protein)
  summarizedProteinIntensities <- left_join(counts, summIntensities, by ="Protein")
  summarizedProteinIntensities <- right_join(annotation, summarizedProteinIntensities, by = "Protein")
  obj <- readMSnSet2(summarizedProteinIntensities,ecol=c(6:ncol(summarizedProteinIntensities)))
  pData(obj) <- pData(data)
  featureNames(obj) <- fData(obj)$Protein
  sampleNames(obj) <- pData(obj)$SampleName
  return(obj)
}


############### Normalization functions ################

# functions for normalizing intensity data

# Performs quantile normalization on the intensities within columns

normalizeQuantiles <- function(data)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  exprs(data) <- normalize.quantiles(exprs(data))
  return(data)
}

# Performs scaling normalization on the intensities within columns (mean, median or sum)
normalizeScaling <- function(data, func, Protein = NULL)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  intensities <- as.data.frame(exprs(data))
  intensitiesForScaling <- intensities
  
  if (!is.null(Protein))
  {
    featuredata <- fData(data)
    ### use protein identifier here
    ind <- which(featuredata$Master.Protein.Accessions == Protein)
    if(length(ind)==0)
      stop('Protein not found')
    intensitiesForScaling <- intensities[ind,] 
  }
    
  scaledIntensities <- intensitiesForScaling %>%
    summarize_all(funs(func)) %>%
    mutate_all(funs(log)) %>%
    as.numeric

  scalingFactors <- exp(scaledIntensities - mean(scaledIntensities))
  normalizedIntensities <- t(t(intensities) / scalingFactors)
  exprs(data) <- normalizedIntensities
  return(data)
}


# Performs scaling normalization on the intensities within group (median or mean)

groupScaling <- function(data,func,Grp="SampleGroup")
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  if(!is.character(Grp))
    stop('Grp has to be of class character..')
  intensities <- as.data.frame(exprs(data))
  allgrps <- split(data,Grp)
  scalingFactors <- numeric()
  for (i in 1:length(allgrps))
  {
    intensitiesForScaling <- as.data.frame(exprs(allgrps[[i]]))
    scaledIntensities <- intensitiesForScaling %>%
      summarize_all(funs(func)) %>%
      mutate_all(funs(log)) %>%
      as.numeric
    grpsFactors <- exp(scaledIntensities - mean(scaledIntensities))
    names(grpsFactors) <- pData(allgrps[[i]])$SampleName
    scalingFactors <- c(scalingFactors,grpsFactors)
  }
  ind <- match(data$SampleName,names(scalingFactors))
  scalingFactors <- scalingFactors[ind]
  normalizedIntensities <- t(t(intensities) / scalingFactors)
  exprs(data) <- normalizedIntensities
  return(data)
}

#### Row scaling based on mean or median of row
rowScaling <- function(data,func)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  intensities <- exprs(data)
  rwm <- apply(intensities,1,func)
  res <- intensities/rwm
  exprs(data) <- log2(res+0.0001)
  return(data)
}


#### Function to regress expression values based on single protein ####
regressIntensity <- function(data,controlInd=NULL,ProteinId)
{
  if(!is(data,"MSnSet"))
    stop('data has to be of class MSnSet..')
  if(!is.null(controlInd))
    data <- data[,-controlInd]
  if(!is.character(ProteinId))
    stop('ProteinId has to be of class character')
  ind <- which(fData(data)$Protein==ProteinId)
  if(length(ind)==0)
    stop('ProteinId is not found or this is not summarized protein intensities dataset...')
  prot <- exprs(data)[ind,]
  dep <- exprs(data)
  indep <- exprs(data)
  for (i in 1:ncol(indep))
  {
    indep[,i] <- prot[i]
  }
  combdata <- cbind(dep,indep)
  Original_Correlation <- apply(dep[-ind,],1,function(x) cor(x,dep[ind,]))
  par(mfrow=c(1,2))
  hist(Original_Correlation,main = "Corr Raw data")
  residuals <- apply(combdata, 1, function (x) resid(lm(x[1:ncol(dep)]~x[(ncol(dep)+1):ncol(combdata)])))
  residuals <- t(residuals)
  exprs(data) <- residuals
  pData(data)$SampleGroup <- factor(pData(data)$SampleGroup)
  reg_dep <- exprs(data)
  Transformed_Correlation <- apply(reg_dep[-ind,],1,function(x) cor(x,dep[ind,]))
  hist(Transformed_Correlation,main = "Corr Regressed data")
  return(data)
}

############## Differential Expression #############


# Fits a linear model to the intensity data using limma.
# The PhenoData table must contain SampleName and SampleGroup columns.
# The intensities table must contain column headings for each sample in PhenoData
# Fits a linear model to the intensity data using limma.
# The PhenoData table must contain SampleName and SampleGroup columns.
# The intensities table must contain column headings for each sample in PhenoData
computeDiffStats <- function(data, batchEffect = NULL, applyLog2Transform = TRUE, contrasts, 
                             trend = TRUE, robust = TRUE)
{
  cat("Fitting linear model\n")
  samples <- as.character(pData(data)$SampleName)
  intensities <- as.data.frame(exprs(data))
  if (applyLog2Transform)
  {
    log2xplus1 <- function(x) { log2(x + 1) }
    intensities <- log2xplus1(intensities) 
  }
  batchEffect <- unique(c("SampleGroup", batchEffect))
  model <- as.formula(paste(c("~ 0", batchEffect), collapse = " + "))
  design <- model.matrix(model, data = pData(data))
  colnames(design) <- colnames(design) %>%
    sub(pattern = "^SampleGroup", replacement = "") %>%
    gsub(pattern = " ", replacement = "_")
  if(length(which(is.na(pData(data)$TechRep)==FALSE))>0)
  {
    dupcor <- duplicateCorrelation(intensities,design,block=pData(data)$TechRep)
    fit <- lmFit(intensities, design = design, weights = NULL, correlation=dupcor$consensus)
  }
  else
    fit <- lmFit(intensities, design = design, weights = NULL)
  
  cat("Fitting contrasts\n")
  
  contrasts <- contrasts %>%
    gsub(pattern = " ", replacement = "_") %>%
    sub(pattern = "_-_", replacement = " - ") %>%
    sub(pattern = "_vs_", replacement = " - ")
  
  contrasts <- makeContrasts(contrasts = contrasts, levels = fit$design)
  contrastsfit <- contrasts.fit(fit, contrasts)
  
  cat("Computing empirical Bayes statistics for differential expression")
  # if trend = TRUE, fit an intensity-dependent trend to prior variances 
  # (not assuming constant variance across intensity range)
  fittedContrasts <- eBayes(contrastsfit, trend = trend, robust = robust)
  return(diffstats <- list(data=data,fittedLM=fit,fittedContrasts=fittedContrasts))
}


getContrastResults <- function(diffstats, contrast, controlGroup = NULL, applyLog2Transform = FALSE, 
                               writeFile= FALSE)
{
  cat("Obtaining results for contrast", contrast, "\n")
  
  contrast <- contrast %>%
    gsub(pattern = " ", replacement = "_") %>%
    sub(pattern = "_-_", replacement = " - ") %>%
    sub(pattern = "_vs_", replacement = " - ")
  
  data <- diffstats$data
  fittedContrasts <- diffstats$fittedContrasts
  fittedLinearModel <- diffstats$fittedLM
  results <- topTable(fittedContrasts, coef = contrast, number = Inf, sort.by = "none")
  contrastGroups <- contrast %>% strsplit(" - ") %>% unlist
  fittedIntensities <- as.data.frame(fittedLinearModel$coefficients)
  contrastIntensities <- select(fittedIntensities, one_of(contrastGroups))
  
  if (!is.null(controlGroup))
  {
    controlIntensity <- fittedIntensities[, controlGroup]
    results$controlLogFoldChange = apply(contrastIntensities - controlIntensity, 1, max)
  }
  
  intensities <- as.data.frame(exprs(data))
  if (applyLog2Transform)
  {
    log2xplus1 <- function(x) { log2(x + 1) }
    intensities <- log2xplus1(intensities) 
  }
  results <- cbind(fData(data),intensities,results)
  samples <- as.character(data$SampleName)
  results <- results %>%
    arrange(desc(B))
  results <- results %>%
    mutate_at(funs(round(., digits = 2)), .vars=c("logFC", "t", "B",samples)) %>%
    mutate_at(funs(signif(., digits = 2)), .vars=c("P.Value", "adj.P.Val"))
  colnames(results)[match(c("AveExpr","logFC","Count"),colnames(results))] <- 
    c("AvgIntensity","log2FC","Unique_Peptides")
  
  if(writeFile == TRUE)
    write.table(results, paste0(names(contrast),".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  return(results)
}