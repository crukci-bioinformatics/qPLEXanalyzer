# Argument check function
checkArg_regressIntensity <- function(MSnSetObj, controlInd, ProteinId){
    assert_that(is_MSnSet(MSnSetObj), is_ProteinSet(MSnSetObj))
    assert_that(is_validControlColumn(controlInd, MSnSetObj))
    assert_that(is.string(ProteinId))
    assert_that(is_validProteinId(ProteinId, MSnSetObj))
}

#' Regression based analysis
#' 
#' Performs linear regression on protein intensities based on selected protein
#' (qPLEX-RIME bait)
#'
#' This function performs regression based analysis upon protein intensities
#' based on a selected protein. In qPLEX RIME this method could be used to
#' regress out the effect of target protein on other interactors. This function
#' corrects this dependency of many proteins on the target protein levels by
#' linear regression. It sets the target protein levels as the independent
#' variable (x) and each of the other proteins as the dependent variable (y).
#' The resulting residuals of the linear regressions y=ax+b are the protein
#' levels corrected for target protein dependency.
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param controlInd numeric; index of IgG within MSnSet
#' @param ProteinId character; Uniprot protein ID
#' @param plot character; Whether or not to plot the QC histograms
#' @return An object of class \code{MSnSet} (see \code{\link{MSnSet-class}}).
#' This consists of corrected protein levels. In addition, the function can
#' also plot histograms of correlation of target protein with all other
#' proteins before and after this correction.
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' MSnset_P <- summarizeIntensities(MSnSet_data, sum, human_anno)
#' IgG_ind <- which(pData(MSnset_P)$SampleGroup == "IgG")
#' MSnset_reg <- regressIntensity(MSnset_P, 
#'                                controlInd=IgG_ind, 
#'                                ProteinId="P03372")
#' 
#' @importFrom Biobase exprs exprs<- fData pData pData<-
#' @importFrom graphics hist par
#' @importFrom stats cor lm resid
#'
#' @export regressIntensity
regressIntensity <- function(MSnSetObj, ProteinId, controlInd=NULL, plot=TRUE) {
    checkArg_regressIntensity(MSnSetObj, controlInd, ProteinId)
    
    if (!is.null(controlInd)) { MSnSetObj <- MSnSetObj[, -controlInd] }
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

