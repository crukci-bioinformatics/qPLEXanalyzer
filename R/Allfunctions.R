##############
# This script contains all functions required for data processing and analysis 
###############

# Log2 with addition of 1e-10 count to deal with zeros
log2xplus1 <- function(x) {
    log2(x + 1e-10)
}

convertToMSnset <- function(ExpObj, metadata, indExpData, Sequences=NULL, 
                            Accessions, type="peptide", rmMissing=TRUE) {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_convertToMSnset, .args)
    
    colnames(ExpObj)[Accessions] <- "Accessions"
    colnames(ExpObj)[Sequences] <- "Sequences"
    
    if (rmMissing) {
        ExpObj %<>% filter_at(vars(indExpData), all_vars(!is.na(.)))
    }
    obj <- readMSnSet2(ExpObj, ecol = indExpData)
    rownames(metadata) <- as.character(metadata$SampleName)
    ind <- match(sampleNames(obj), rownames(metadata))
    metadata <- metadata[ind, ]
    pData(obj) <- metadata
    if (type == "protein") {
        featureNames(obj) <- fData(obj)$Accessions
    } else {
        featureNames(obj) <- paste0("peptide_", featureNames(obj))
    }
    return(obj)
}



##### Summarization function ##########

# Functions for summarizing intensity data

# Summarizes multiple peptide measurements for a protein.
# Assumes that there are columns for each of the samples specified and uses
# Accessions column for grouping peptide-level measurements.
# Filters any rows with missing values.
# Typical summarization functions are sum, mean and median.
# For successful running of this function the annotation file must have at least
# two columns: "Accessions" and "GeneSymbol"
# The MSnSetObj must have column "Sequences" denoting its a peptide dataset

summarizeIntensities <- function(MSnSetObj, summarizationFunction, annotation) {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_summarizeIntensities, .args)
    
    counts <- as.data.frame(fData(MSnSetObj)) %>%
        select(Accessions, Sequences) %>%
        distinct() %>%
        count(Accessions) %>%
        mutate(Accessions = as.character(Accessions)) %>%
        rename(Count = n)
    
    summarizedIntensities <- as.data.frame(exprs(MSnSetObj)) %>%
        mutate(Accessions = as.character(fData(MSnSetObj)$Accessions)) %>%
        group_by(Accessions) %>%
        summarize_all(funs(summarizationFunction)) %>%
        left_join(counts, by = "Accessions") %>%
        left_join(annotation, by = "Accessions") %>%
        select(Accessions, colnames(annotation), Count, everything())
    
    expInd <- seq(
        grep("Count", colnames(summarizedIntensities)) + 1,
        ncol(summarizedIntensities)
    )
    
    obj <- readMSnSet2(summarizedIntensities, ecol = expInd)
    pData(obj) <- pData(MSnSetObj)
    featureNames(obj) <- fData(obj)$Accessions
    sampleNames(obj) <- pData(obj)$SampleName
    return(obj)
}



############### Normalization functions ################

# functions for normalizing intensity data

# Performs quantile normalization on the intensities within columns

normalizeQuantiles <- function(MSnSetObj) {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_normalizeQuantiles, .args)
    
    exprs(MSnSetObj) <- normalize.quantiles(exprs(MSnSetObj))
    return(MSnSetObj)
}

# Performs scaling normalization on the intensities within columns
normalizeScaling <- function(MSnSetObj, scalingFunction, ProteinId = NULL) {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_normalizeScaling, .args)
    
    intensities <- as.data.frame(exprs(MSnSetObj))
    intensitiesForScaling <- intensities
    
    if (!is.null(ProteinId)) {
        featuredata <- fData(MSnSetObj)
        ### use protein identifier here
        ind <- which(featuredata$Accessions == ProteinId)
        if (!length(ind)) {
            stop("Protein not found")
        }
        intensitiesForScaling <- intensities[ind, ]
    }
    
    scaledIntensities <- intensitiesForScaling %>%
        summarize_all(funs(scalingFunction)) %>%
        mutate_all(funs(log)) %>%
        as.numeric()
    
    scalingFactors <- exp(scaledIntensities - mean(scaledIntensities))
    normalizedIntensities <- t(t(intensities) / scalingFactors)
    exprs(MSnSetObj) <- normalizedIntensities
    return(MSnSetObj)
}


# Performs scaling normalization on the intensities within group 

groupScaling <- function(MSnSetObj, scalingFunction=median, 
                         groupingColumn="SampleGroup") {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_groupScaling, .args)
    
    exprs(MSnSetObj) <- as.data.frame(exprs(MSnSetObj)) %>%
        rownames_to_column("PeptideID") %>%
        gather("SampleName", "RawIntensity", -PeptideID) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        rename_at(vars(groupingColumn), ~ "Grouping_column") %>%
        group_by(SampleName) %>%
        mutate(scaledIntensity = scalingFunction(RawIntensity) %>% log()) %>%
        group_by(Grouping_column) %>%
        mutate(meanscaledIntensity = mean(scaledIntensity)) %>%
        ungroup() %>%
        mutate(scalingFactors = exp(scaledIntensity - meanscaledIntensity)) %>%
        mutate(normalizedIntensities = RawIntensity / scalingFactors) %>%
        select(PeptideID, SampleName, normalizedIntensities) %>%
        spread(SampleName, normalizedIntensities) %>%
        arrange(factor(PeptideID, levels = rownames(MSnSetObj))) %>%
        as.data.frame() %>%
        column_to_rownames("PeptideID") %>%
        select(colnames(MSnSetObj)) %>%
        as.matrix()
    return(MSnSetObj)
}

#### Row scaling based on mean or median of row
rowScaling <- function(MSnSetObj, scalingFunction) {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_rowScaling, .args)
    
    intensities <- exprs(MSnSetObj)
    rwm <- apply(intensities, 1, scalingFunction)
    res <- intensities / rwm
    exprs(MSnSetObj) <- log2(res + 0.0001)
    return(MSnSetObj)
}


#### Function to regress expression values based on single protein ####
regressIntensity <- function(MSnSetObj, controlInd=NULL, ProteinId) {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_regressIntensity, .args)
    
    ind <- which(fData(MSnSetObj)$Accessions == ProteinId)
    prot <- exprs(MSnSetObj)[ind, ]
    dep <- exprs(MSnSetObj)
    indep <- exprs(MSnSetObj)
    for (i in seq_len(ncol(indep))) {
        indep[, i] <- prot[i]
    }
    combdata <- cbind(dep, indep)
    originalCorrelation <- apply(dep[-ind, ], 1, function(x) cor(x, dep[ind, ]))
    par(mfrow = c(1, 2))
    hist(originalCorrelation, main = "Corr Raw data")
    calculateResiduals <- function(x){
        resid(lm(x[seq_len(ncol(dep))] ~ x[seq(ncol(dep) + 1, ncol(combdata))]))
    }
    residuals <- apply(combdata, 1, calculateResiduals)
    residuals <- t(residuals)
    exprs(MSnSetObj) <- residuals
    pData(MSnSetObj)$SampleGroup <- factor(pData(MSnSetObj)$SampleGroup)
    reg_dep <- exprs(MSnSetObj)
    transformedCorrelation <- apply(reg_dep[-ind, ], 1, 
                                     function(x) cor(x, dep[ind, ]))
    hist(transformedCorrelation, main = "Corr Regressed data")
    return(MSnSetObj)
}

############## Differential Expression #############


# Fits a linear model to the intensity data using limma.
# The PhenoData table must contain SampleName and SampleGroup columns.
# The intensities table must contain a column for each sample in PhenoData
computeDiffStats <- function(MSnSetObj, batchEffect = NULL, transform = TRUE, 
                             contrasts, trend = TRUE, robust = TRUE) {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_computeDiffStats, .args)
    
    message("Fitting linear model")
    intensities <- as.data.frame(exprs(MSnSetObj))
    if (transform) {
        intensities <- log2xplus1(intensities)
    }
    batchEffect <- unique(c("SampleGroup", batchEffect))
    model <- as.formula(paste(c("~ 0", batchEffect), collapse = " + "))
    design <- model.matrix(model, data = pData(MSnSetObj))
    colnames(design) <- colnames(design) %>%
        sub(pattern = "^SampleGroup", replacement = "") %>%
        gsub(pattern = " ", replacement = "_")
    if (length(which(is.na(pData(MSnSetObj)$TechRep) == FALSE)) > 0) {
        dupcor <- duplicateCorrelation(intensities, 
                                       design = design, 
                                       block = pData(MSnSetObj)$TechRep)
        fit <- lmFit(intensities, 
                     design = design, 
                     weights = NULL, 
                     correlation = dupcor$consensus)
    } else {
        fit <- lmFit(intensities, design = design, weights = NULL)
    }
    
    message("Fitting contrasts\n")
    
    contrasts <- contrasts %>%
        gsub(pattern = " ", replacement = "_") %>%
        sub(pattern = "_-_", replacement = " - ") %>%
        sub(pattern = "_vs_", replacement = " - ")
    
    contrasts <- makeContrasts(contrasts = contrasts, levels = fit$design)
    contrastsfit <- contrasts.fit(fit, contrasts)
    
    message("Computing empirical Bayes statistics for differential expression")
    fittedContrasts <- eBayes(contrastsfit, trend = trend, robust = robust)
    return(diffstats <- list(MSnSetObj = MSnSetObj, 
                             fittedLM = fit, 
                             fittedContrasts = fittedContrasts))
}


getContrastResults <- function(diffstats, contrast, controlGroup = NULL, 
                               transform = TRUE, writeFile= FALSE) {
    .args <- as.list(match.call()[-1])
    do.call(checkArg_getContrastResults, .args)
    
    message("Obtaining results for contrast", contrast, "\n")
    contrast <- contrast %>%
        gsub(pattern = " ", replacement = "_") %>%
        sub(pattern = "_-_", replacement = " - ") %>%
        sub(pattern = "_vs_", replacement = " - ")
    
    MSnSetObj <- diffstats$MSnSetObj
    fittedContrasts <- diffstats$fittedContrasts
    fittedLinearModel <- diffstats$fittedLM
    results <- topTable(fittedContrasts, 
                        coef = contrast, 
                        number = Inf, 
                        sort.by = "none")
    contrastGroups <- contrast %>% strsplit(" - ") %>% unlist()
    fittedIntensities <- as.data.frame(fittedLinearModel$coefficients)
    contrastIntensities <- select(fittedIntensities, one_of(contrastGroups))
    
    if (!is.null(controlGroup)) {
        controlIntensity <- fittedIntensities[, controlGroup]
        results$controlLogFoldChange <- 
            apply(contrastIntensities - controlIntensity, 1, max)
    }
    
    intensities <- as.data.frame(exprs(MSnSetObj))
    if (transform) {
        intensities <- log2xplus1(intensities)
    }
    SamplesCol <- as.character(MSnSetObj$SampleName)
    results <- cbind(fData(MSnSetObj), intensities, results)
    results <- results %>%
        arrange(desc(B))
    results <- results %>%
        mutate_at(funs(round(., digits = 2)), 
                  .vars = c("logFC", "t", "B", SamplesCol)) %>%
        mutate_at(funs(signif(., digits = 2)), 
                  .vars = c("P.Value", "adj.P.Val")) %>% 
        rename_all(function(x){str_replace(x, 
                                           "^Count$", 
                                           "Unique_peptides")}) %>%
        rename(AvgIntensity=AveExpr, log2FC=logFC)
    
    if (writeFile == TRUE) {
        write.table(results, paste0(names(contrast), ".txt"), 
                    quote = FALSE, 
                    sep = "\t", 
                    row.names = FALSE)
    }
    return(results)
}
