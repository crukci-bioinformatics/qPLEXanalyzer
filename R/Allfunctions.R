######### This script contains all the functions required for data processing and analysis ###############

# Log2 with addition of 1e-10 count to deal with zeros
log2xplus1 <- function(x) { log2(x + 1e-10) }

convertToMSnset <- function(ExpObj,metadata,indExpData,Sequences=NULL,Accessions,
                            type="peptide",rmMissing=TRUE){
  if (!is.data.frame(ExpObj)){ stop("ExpObj has to be of class dataframe") }
  if (!is.data.frame(metadata)){ stop("metadata has to be of class dataframe") }
  if(!is.numeric(indExpData)){ stop('indExpData has to be of class numeric ..') }
  if(!is.numeric(Sequences) & !is.null(Sequences)){ stop('Sequences has to be either numeric or NULL..') }
  if(!is.numeric(Accessions)){ stop('Accessions has to be of class numeric ..') }
  if(!is.character(type)){ stop('type has to be of class character ..') }
  if(!(type%in%c("peptide","protein"))) stop('type has to be either peptide or protein')
  columns <- c("SampleName","SampleGroup","BioRep","TechRep")
  if(!all(columns%in%colnames(metadata))) stop('metadata must contain"', columns, '" columns ..')
  if (!is.logical(rmMissing)){ stop("rmMissing has to be of class logical") }
  
  colnames(ExpObj)[Accessions] <- "Accessions"
  if(type=="peptide") 
  {
    if(!is.numeric(Sequences)){ stop('Sequences has to be numeric when type is peptide..') }
    colnames(ExpObj)[Sequences] <- "Sequences"
  }
  if(rmMissing){ ExpObj %<>% filter_at(vars(indExpData), all_vars(!is.na(.))) }
  obj <- readMSnSet2(ExpObj,ecol=indExpData)
  rownames(metadata) <- as.character(metadata$SampleName)
  check <- all(rownames(metadata) %in% sampleNames(obj))
  if(check){
    ind <- match(sampleNames(obj), rownames(metadata))
    metadata <- metadata[ind,]
  }else{
    stop("Column names in expression data do not match SampleNames in metadata")
  }
  pData(obj) <- metadata
  if(type=="protein")
    featureNames(obj) <- fData(obj)$Accessions
  else
    featureNames(obj) <- paste0("peptide_",featureNames(obj))
  return(obj)
}   

  

##### Summarization function ##########

# Functions for summarizing intensity data

# Summarizes multiple peptide measurements for a protein.
# Assumes that there are columns for each of the samples specified and uses
# Accessions column for grouping peptide-level measurements.
# Filters any rows with missing values.
# Typical summarization functions are sum, mean and median.
# For successful running of this function the annotation file must have four column "Accessions","Gene","Description" and "GeneSymbol"
# In addition the MSnSetObj must have column "Sequences" denoting its a peptide dataset

summarizeIntensities <- function(MSnSetObj, summarizationFunction, annotation){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!is.data.frame(annotation)){ stop('annotation has to be of class data frame..') }
  if(!"Sequences" %in% colnames(fData(MSnSetObj))){ stop('This MSnSet is not a peptide dataset ..') }
  
  counts <- as.data.frame(fData(MSnSetObj)) %>%
    select(Accessions, Sequences) %>% 
    distinct %>% 
    count(Accessions) %>%
    mutate(Accessions=as.character(Accessions)) %>% 
    rename(Count = n)
  
  summarizedIntensities <- as.data.frame(exprs(MSnSetObj)) %>%
    mutate(Accessions=as.character(fData(MSnSetObj)$Accessions)) %>%
    group_by(Accessions) %>%
    summarize_all(funs(summarizationFunction)) %>% 
    left_join(counts, by="Accessions") %>% 
    left_join(annotation, by = "Accessions") %>% 
    select(Accessions, colnames(annotation), Count, everything())
  
  expInd <- seq(grep("Count", colnames(summarizedIntensities)) + 1,
                ncol(summarizedIntensities))
  
  obj <- readMSnSet2(summarizedIntensities, ecol=expInd)
  pData(obj) <- pData(MSnSetObj)
  featureNames(obj) <- fData(obj)$Accessions
  sampleNames(obj) <- pData(obj)$SampleName
  return(obj)
}



############### Normalization functions ################

# functions for normalizing intensity data

# Performs quantile normalization on the intensities within columns

normalizeQuantiles <- function(MSnSetObj){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  exprs(MSnSetObj) <- normalize.quantiles(exprs(MSnSetObj))
  return(MSnSetObj)
}

# Performs scaling normalization on the intensities within columns (mean, median or sum)
normalizeScaling <- function(MSnSetObj, scalingFunction, ProteinId = NULL){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  intensities <- as.data.frame(exprs(MSnSetObj))
  intensitiesForScaling <- intensities
  
  if (!is.null(ProteinId)){
    featuredata <- fData(MSnSetObj)
    ### use protein identifier here
    ind <- which(featuredata$Accessions == ProteinId)
    if(!length(ind))
      stop('Protein not found')
    intensitiesForScaling <- intensities[ind,] 
  }
    
  scaledIntensities <- intensitiesForScaling %>%
    summarize_all(funs(scalingFunction)) %>%
    mutate_all(funs(log)) %>%
    as.numeric

  scalingFactors <- exp(scaledIntensities - mean(scaledIntensities))
  normalizedIntensities <- t(t(intensities) / scalingFactors)
  exprs(MSnSetObj) <- normalizedIntensities
  return(MSnSetObj)
}


# Performs scaling normalization on the intensities within group (median or mean)

groupScaling <- function(MSnSetObj, scalingFunction=median, groupingColumn="SampleGroup"){
    if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
    if(!is.character(groupingColumn)){ stop('groupingColumn has to be of class character..') }
    if(!is.function(scalingFunction)){ stop('scalingFunction should be a summary function, e.g. mean or median') }

    exprs(MSnSetObj) <- as.data.frame(exprs(MSnSetObj)) %>%
        rownames_to_column("PeptideID")  %>%
        gather("SampleName", "RawIntensity", -PeptideID) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        rename_at(vars(groupingColumn), ~"Grouping_column") %>%
        group_by(SampleName) %>%
        mutate(scaledIntensity=scalingFunction(RawIntensity) %>% log) %>%
        group_by(Grouping_column) %>%
        mutate(meanscaledIntensity=mean(scaledIntensity)) %>%
        ungroup() %>%
        mutate(scalingFactors=exp(scaledIntensity-meanscaledIntensity)) %>%
        mutate(normalizedIntensities=RawIntensity/scalingFactors) %>%
        select(PeptideID, SampleName, normalizedIntensities) %>%
        spread(SampleName, normalizedIntensities) %>%
        arrange(factor(PeptideID, levels=rownames(MSnSetObj))) %>%
        as.data.frame() %>%
        column_to_rownames("PeptideID") %>%
        select(colnames(MSnSetObj)) %>%
        as.matrix()
    return(MSnSetObj)
}

#### Row scaling based on mean or median of row
rowScaling <- function(MSnSetObj, scalingFunction){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  intensities <- exprs(MSnSetObj)
  rwm <- apply(intensities,1,scalingFunction)
  res <- intensities/rwm
  exprs(MSnSetObj) <- log2(res+0.0001)
  return(MSnSetObj)
}


#### Function to regress expression values based on single protein ####
regressIntensity <- function(MSnSetObj,controlInd=NULL,ProteinId){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!is.null(controlInd)){ MSnSetObj <- MSnSetObj[,-controlInd] }
  if(!is.character(ProteinId)){ stop('ProteinId has to be of class character') }
  ind <- which(fData(MSnSetObj)$Accessions==ProteinId)
  if(!length(ind)){ stop('ProteinId is not found or this is not summarized protein intensities dataset...') }
  prot <- exprs(MSnSetObj)[ind,]
  dep <- exprs(MSnSetObj)
  indep <- exprs(MSnSetObj)
  for (i in seq_len(ncol(indep))){
    indep[,i] <- prot[i]
  }
  combdata <- cbind(dep,indep)
  Original_Correlation <- apply(dep[-ind,],1,function(x) cor(x,dep[ind,]))
  par(mfrow=c(1,2))
  hist(Original_Correlation,main = "Corr Raw data")
  residuals <- apply(combdata, 1, function (x) resid(lm(x[seq_len(ncol(dep))]~x[seq(ncol(dep)+1, ncol(combdata))])))
  residuals <- t(residuals)
  exprs(MSnSetObj) <- residuals
  pData(MSnSetObj)$SampleGroup <- factor(pData(MSnSetObj)$SampleGroup)
  reg_dep <- exprs(MSnSetObj)
  Transformed_Correlation <- apply(reg_dep[-ind,],1,function(x) cor(x,dep[ind,]))
  hist(Transformed_Correlation,main = "Corr Regressed data")
  return(MSnSetObj)
}

############## Differential Expression #############


# Fits a linear model to the intensity data using limma.
# The PhenoData table must contain SampleName and SampleGroup columns.
# The intensities table must contain column headings for each sample in PhenoData
# Fits a linear model to the intensity data using limma.
# The PhenoData table must contain SampleName and SampleGroup columns.
# The intensities table must contain column headings for each sample in PhenoData
computeDiffStats <- function(MSnSetObj, batchEffect = NULL, transform = TRUE, contrasts, 
                             trend = TRUE, robust = TRUE){
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  if(!is.logical(transform)){ stop('transform has to either TRUE or FALSE..') }
  if(!is.logical(trend)){ stop('trend has to be of either TRUE or FALSE..') }
  if(!is.logical(robust)){ stop('robust has to be of either TRUE or FALSE..') }
  if(!is(MSnSetObj,"MSnSet")){ stop('MSnSetObj has to be of class MSnSet..') }
  message("Fitting linear model\n")
  intensities <- as.data.frame(exprs(MSnSetObj))
  if (transform){
    intensities <- log2xplus1(intensities) 
  }
  batchEffect <- unique(c("SampleGroup", batchEffect))
  model <- as.formula(paste(c("~ 0", batchEffect), collapse = " + "))
  design <- model.matrix(model, data = pData(MSnSetObj))
  colnames(design) <- colnames(design) %>%
    sub(pattern = "^SampleGroup", replacement = "") %>%
    gsub(pattern = " ", replacement = "_")
  if(length(which(is.na(pData(MSnSetObj)$TechRep)==FALSE))>0){
    dupcor <- duplicateCorrelation(intensities,design,block=pData(MSnSetObj)$TechRep)
    fit <- lmFit(intensities, design = design, weights = NULL, correlation=dupcor$consensus)
  }else{
    fit <- lmFit(intensities, design = design, weights = NULL)
  }
  
  message("Fitting contrasts\n")
  
  contrasts <- contrasts %>%
    gsub(pattern = " ", replacement = "_") %>%
    sub(pattern = "_-_", replacement = " - ") %>%
    sub(pattern = "_vs_", replacement = " - ")
  
  contrasts <- makeContrasts(contrasts = contrasts, levels = fit$design)
  contrastsfit <- contrasts.fit(fit, contrasts)
  
  message("Computing empirical Bayes statistics for differential expression\n")
  fittedContrasts <- eBayes(contrastsfit, trend = trend, robust = robust)
  return(diffstats <- list(MSnSetObj=MSnSetObj,fittedLM=fit,fittedContrasts=fittedContrasts))
}


getContrastResults <- function(diffstats, contrast, controlGroup = NULL, transform = TRUE, 
                               writeFile= FALSE){
  if(!is.logical(transform)){ stop('transform has to either TRUE or FALSE..') }
  if(!is.logical(writeFile)){ stop('writeFile has to either TRUE or FALSE..') }
  message("Obtaining results for contrast", contrast, "\n")
  
  contrast <- contrast %>%
    gsub(pattern = " ", replacement = "_") %>%
    sub(pattern = "_-_", replacement = " - ") %>%
    sub(pattern = "_vs_", replacement = " - ")
  
  MSnSetObj <- diffstats$MSnSetObj
  fittedContrasts <- diffstats$fittedContrasts
  fittedLinearModel <- diffstats$fittedLM
  results <- topTable(fittedContrasts, coef = contrast, number = Inf, sort.by = "none")
  contrastGroups <- contrast %>% strsplit(" - ") %>% unlist
  fittedIntensities <- as.data.frame(fittedLinearModel$coefficients)
  contrastIntensities <- select(fittedIntensities, one_of(contrastGroups))
  
  if (!is.null(controlGroup)){
    controlIntensity <- fittedIntensities[, controlGroup]
    results$controlLogFoldChange = apply(contrastIntensities - controlIntensity, 1, max)
  }
  
  intensities <- as.data.frame(exprs(MSnSetObj))
  if (transform){
    intensities <- log2xplus1(intensities) 
  }
  SamplesCol <- as.character(MSnSetObj$SampleName)
  results <- cbind(fData(MSnSetObj),intensities,results)
  results <- results %>%
    arrange(desc(B))
  results <- results %>%
    mutate_at(funs(round(., digits = 2)), .vars=c("logFC", "t", "B",SamplesCol)) %>%
    mutate_at(funs(signif(., digits = 2)), .vars=c("P.Value", "adj.P.Val"))
  colnames(results)[match(c("Count","AveExpr","logFC"),colnames(results))] <- c("Unique_peptides",
                                                                                "AvgIntensity","log2FC")
  if(writeFile == TRUE)
    write.table(results, paste0(names(contrast),".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  return(results)
}