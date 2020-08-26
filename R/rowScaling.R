# Row scaling based on mean or median of row
rowScaling <- function(MSnSetObj, scalingFunction) {
    checkArg_rowScaling(MSnSetObj, scalingFunction)
    
    intensities <- exprs(MSnSetObj)
    rwm <- apply(intensities, 1, scalingFunction)
    res <- intensities / rwm
    exprs(MSnSetObj) <- log2(res + 0.0001)
    return(MSnSetObj)
}

