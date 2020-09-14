# Performs quantile normalization on the intensities within columns



#' Quantile normalization
#' 
#' Performs quantile normalization on the intensities within columns
#' 
#' The peptide intensities are roughly replaced by the order statics on their
#' abundance.  This normalization technique has the effect of making the
#' distributions of intensities from the different samples identical in terms
#' of their statistical properties. It is the strongest normalization method
#' and should be used carefully as it erases most of the difference between the
#' samples.
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @return An object of class \code{MSnSet} (see \code{\link{MSnSet-class}}) 
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' MSnset_norm <- normalizeQuantiles(MSnSet_data)
#' 
#' @export normalizeQuantiles
normalizeQuantiles <- function(MSnSetObj) {
    checkArg_normalizeQuantiles(MSnSetObj)
    
    exprs(MSnSetObj) <- normalize.quantiles(exprs(MSnSetObj))
    return(MSnSetObj)
}

