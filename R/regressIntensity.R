#### Function to regress expression values based on single protein ####

regressIntensity <- function(MSnSetObj, ProteinId, controlInd=NULL, plot=TRUE) {
    checkArg_regressIntensity(MSnSetObj, controlInd, ProteinId)
    
    ind <- which(fData(MSnSetObj)$Accessions == ProteinId)
    prot <- exprs(MSnSetObj)[ind, ]
    dep <- exprs(MSnSetObj)
    indep <- exprs(MSnSetObj)
    for (i in seq_len(ncol(indep))) {
        indep[, i] <- prot[i]
    }
    combdata <- cbind(dep, indep)
    originalCorrelation <- apply(dep[-ind, ], 1, function(x) cor(x, dep[ind, ]))
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
    if(plot){
        tryCatch({
            par(mfrow = c(1, 2))
            hist(originalCorrelation, main = "Corr Raw data")
            hist(transformedCorrelation, main = "Corr Regressed data")
            }, 
            error = function(err) {
                message("QC histograms not plotted:")
                message(err)
                }
            )
    }
    return(MSnSetObj)
}

