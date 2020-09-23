#' Get differential statistics results
#' 
#' Get differential statistics results for given contrasts.
#' 
#' 
#' @param diffstats list; output of computeDiffStats function
#' @param contrast character; contrast of interest for which to retrieve
#' differential statistics results
#' @param controlGroup character; control group such as IgG
#' @param transform logical; apply log2 transformation to the raw intensities
#' @param writeFile logical; whether to write the results into a text file
#' @return A \code{\link{data.frame}} object and text file containing the
#' result of the differential statistics.
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' MSnset_norm <- groupScaling(MSnSet_data, scalingFunction=median)
#' MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
#' contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
#' diffstats <- computeDiffStats(MSnset_Pnorm, contrasts=contrasts)
#' diffexp <- getContrastResults(diffstats=diffstats, contrast=contrasts)
#' 
#' @import limma
#' @importFrom Biobase exprs fData
#' @importFrom dplyr across arrange bind_cols desc mutate rename rename_with
#' select
#' @importFrom magrittr %>%
#' @importFrom readr write_tsv
#' @importFrom stringr str_replace str_replace_all
#' @importFrom tidyselect one_of
#'
#' @export getContrastResults
getContrastResults <- function(diffstats, contrast, controlGroup = NULL, 
                               transform = TRUE, writeFile= FALSE) {
    checkArg_getContrastResults(diffstats, contrast, controlGroup, transform,
                                writeFile)
    
    message("Obtaining results for contrast", contrast, "\n")
    contrast <- contrast %>%
        str_replace_all(pattern = " ", replacement = "_") %>%
        str_replace(pattern = "_-_", replacement = " - ") %>%
        str_replace(pattern = "_vs_", replacement = " - ")
    
    MSnSetObj <- diffstats$MSnSetObj
    fittedContrasts <- diffstats$fittedContrasts
    fittedLinearModel <- diffstats$fittedLM
    results <- topTable(fittedContrasts, 
                        coef = contrast, 
                        number = Inf, 
                        sort.by = "none",
                        confint=TRUE)
    contrastGroups <- contrast %>% strsplit(" - ") %>% unlist()
    fittedIntensities <- as.data.frame(fittedLinearModel$coefficients)
    contrastIntensities <- select(fittedIntensities, one_of(contrastGroups))
    
    if (!is.null(controlGroup)) {
        controlIntensity <- fittedIntensities[, controlGroup]
        results$controlLogFoldChange <- 
            apply(contrastIntensities - controlIntensity, 1, max)
    }
    
    intensities <- as.data.frame(exprs(MSnSetObj))
    if (transform) {
        intensities <- log2xplus1(intensities)
    }
    SamplesCol <- as.character(MSnSetObj$SampleName)
    results <- bind_cols(fData(MSnSetObj), intensities, results) %>%
        arrange(desc(B)) %>%
        mutate(across(c("logFC", "t", "B", SamplesCol), round, digits = 2)) %>%
        mutate(across(c("P.Value", "adj.P.Val"), signif, digits = 2)) %>% 
        # Count column may not be present
        rename_with(str_replace, 
                    pattern = "^Count$", 
                    replacement = "Unique_peptides") %>%
        rename(AvgIntensity=AveExpr, log2FC=logFC)
    
    if (writeFile == TRUE) {
        write_tsv(results, str_c(names(contrast), ".txt"))
    }
    return(results)
}
