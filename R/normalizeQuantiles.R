# Performs quantile normalization on the intensities within columns

normalizeQuantilesnormalizeQuantiles <- function(MSnSetObj) {
    checkArg_normalizeQuantiles(MSnSetObj)
    
    exprs(MSnSetObj) <- normalize.quantiles(exprs(MSnSetObj))
    return(MSnSetObj)
}

