# Argument check function
#' Mean variance plot
#' 
#' Computes and plots variance v mean intensity for peptides in MSnset
#'
#' The intensities are transformed with the natural logarithm after adding
#' 0.001 (\code{log(x + 0.001)}), and the mean and variance of each feature are
#' computed across all samples. A smoothing spline fitted to variance against
#' mean is shown in red. The intensities must not contain missing values.
#'
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param title character: title for the plot
#' @return An object created by \code{ggplot}
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' plotMeanVar(MSnSet_data, title="Mean_Variance")
#' 
#' @import ggplot2
#' @importFrom Biobase exprs
#' @importFrom stats smooth.spline var
#'
#' @export plotMeanVar
plotMeanVar <- function(MSnSetObj, title="") {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    
    intensities <- log(exprs(MSnSetObj) + 0.001)
    mvDat <- data.frame(
        Mean = rowMeans(intensities),
        Variance = apply(intensities, 1, var)
    )
    
    ss <- smooth.spline(x = mvDat$Mean, y = mvDat$Variance, spar = 1)
    ssDat <- data.frame(x = ss$x, y = ss$y)
    
    ggplot(mvDat, aes(x = Mean, y = Variance)) +
        geom_point(size = 0.5, alpha = 0.6, colour = "darkblue") +
        geom_line(data = ssDat, aes(x = x, y = y),
                  colour = "red", linewidth = 0.5) +
        theme_bw() +
        labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5))
}

