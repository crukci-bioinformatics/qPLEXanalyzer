#' Correlation plot
#' 
#' Computes and display correlation plot for samples within MSnSet
#' 
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param addValues logical: adds correlation values to the plot
#' @param title character; title of the plot
#' @param low_cor_colour colour; colour for lowest correlation in scale
#' @param high_cor_colour colour; colour for highest correlation in scale
#' @param low_cor_limit numeric; lower limit for correlation in colour scale
#' @param high_cor_limit numeric; upper limit for correlation in colour scale
#' @param textsize integer: set the size of correlation values text
#' @return An object created by \code{ggplot}
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
#'     metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'     indExpData=c(7:16), 
#'     Sequences=2, 
#'     Accessions=6)
#' corrPlot(MSnSet_data, addValues=TRUE, title="Correlation plot")
#' 
#' # change colours
#' corrPlot(MSnSet_data, addValues=TRUE, title="Correlation plot", 
#'     low_cor_colour="yellow", high_cor_colour="pink")
#' 
#' @import ggplot2
#' @importFrom Biobase exprs
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @importFrom scales squish
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#'
#' @export corrPlot
corrPlot <- function(MSnSetObj, addValues=TRUE, title="", 
                     low_cor_colour="#FFFFFF", high_cor_colour="#B90505",
                     low_cor_limit=0, high_cor_limit=1,
                     textsize=3) {
    checkArg_corrPlot(MSnSetObj, addValues, title, 
                      low_cor_colour, high_cor_colour,
                      low_cor_limit, high_cor_limit,
                      textsize)
    
    plotDat <- cor(exprs(MSnSetObj)) %>%
        as.data.frame() %>%
        rownames_to_column("X") %>%
        pivot_longer(names_to = "Y", values_to = "Cor", -X) %>%
        mutate(AddValues = addValues) %>%
        mutate(CorTxt = ifelse(AddValues, round(Cor, 3), ""))
    
    col2Cols <- c(low_cor_colour, high_cor_colour)
    
    fillScaleLimits <- c(low_cor_limit, high_cor_limit)
    if(max(plotDat$Cor) > high_cor_limit | min(plotDat$Cor) < low_cor_limit){
        warning(str_c("Some correlation values fall outside the provided scale ",
                      "limits.\n",
                      "Out of bounds values will be coloured according to the ",
                      "nearest limit."))
    }
    
    ggplot(plotDat, aes(x = X, y = Y, fill = Cor)) +
        geom_tile(col = "grey") +
        geom_text(aes(label = CorTxt), size=textsize) +
        scale_fill_gradientn(colors = col2Cols, 
                             breaks = seq(0, 1, 0.2),
                             limits = fillScaleLimits, 
                             oob = squish) +
        labs(x = NULL, y = NULL, title = title) +
        guides(fill = guide_colorbar(barheight = 10)) +
        theme(
            aspect.ratio = 1,
            axis.text.x = element_text(angle = 90, 
                                       hjust = 1, 
                                       vjust = 0.5, 
                                       size = 13),
            axis.text.y = element_text(size = 13),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank()
        )
}

