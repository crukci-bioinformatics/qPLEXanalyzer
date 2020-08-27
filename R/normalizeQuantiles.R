# Performs quantile normalization on the intensities within columns

normalizeQuantiles <- function(MSnSetObj) {
    checkArg_normalizeQuantiles(MSnSetObj)
    
    exprs(MSnSetObj) <- normalize.quantiles(exprs(MSnSetObj))
    return(MSnSetObj)
}

