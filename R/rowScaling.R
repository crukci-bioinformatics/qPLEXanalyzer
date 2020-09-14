# Row scaling based on mean or median of row


#' Normalization by scaling peptide/protein intensity across all samples
#' 
#' Divide each peptide/protein by the row mean/median and transform to log2
#' 
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param scalingFunction function; median or mean
#' @return An object of class \code{MSnSet} (see \code{\link{MSnSet-class}}).
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' MSnset_norm <- rowScaling(MSnSet_data, scalingFunction=median)
#' 
#' @export rowScaling
rowScaling <- function(MSnSetObj, scalingFunction) {
    checkArg_rowScaling(MSnSetObj, scalingFunction)
    
    intensities <- exprs(MSnSetObj)
    rwm <- apply(intensities, 1, scalingFunction)
    res <- intensities / rwm
    exprs(MSnSetObj) <- log2(res + 0.0001)
    return(MSnSetObj)
}

