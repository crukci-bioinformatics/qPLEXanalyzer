# Fits a linear model to the intensity data using limma.
# The PhenoData table must contain SampleName and SampleGroup columns.
# The intensities table must contain a column for each sample in PhenoData
# Log2 with addition of 1e-10 count to deal with zeros
log2xplus1 <- function(x) {
    log2(x + 1e-10)
}

computeDiffStats <- function(MSnSetObj, batchEffect = NULL, transform = TRUE, 
                             contrasts, trend = TRUE, robust = TRUE) {
    checkArg_computeDiffStats(MSnSetObj, batchEffect, transform, contrasts,
                              trend, robust)
    
    message("Fitting linear model")
    intensities <- as.data.frame(exprs(MSnSetObj))
    if (transform) {
        intensities <- log2xplus1(intensities)
    }
    batchEffect <- unique(c("SampleGroup", batchEffect))
    model <- as.formula(paste(c("~0", batchEffect), collapse = " + "))
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

