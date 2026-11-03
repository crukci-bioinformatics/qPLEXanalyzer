# Argument check function
checkArg_getContrastResults <- function(diffstats, 
                                        contrast, 
                                        controlGroup, 
                                        transform, 
                                        writeFile){
    assert_that(is_validDiffstats(diffstats))
    assert_that(is_validContrast(contrast, diffstats))
    assert_that(is_validControlGroup(controlGroup, diffstats))
    assert_that(is.flag(transform))
    assert_that(is.flag(writeFile))
}

#' Get differential statistics results
#' 
#' Get differential statistics results for given contrasts.
#'
#' The results for the contrast are obtained with
#' \code{\link[limma]{topTable}} (with confidence intervals) and combined with
#' the feature data and the sample intensities. Rows are sorted by decreasing
#' B statistic.
#'
#' If \code{controlGroup} is given, a \code{controlLogFoldChange} column is
#' added. It is computed from the group means of the fitted linear model: for
#' each group involved in the contrast the log2 fold change over the control
#' group is calculated, and the maximum of these is reported.
#'
#' @param diffstats list; output of computeDiffStats function
#' @param contrast character; name of contrast of interest for which to retrieve
#' differential statistics results
#' @param controlGroup character; control group such as IgG
#' @param transform logical; apply log2 transformation (\code{log2(x + 1e-10)})
#' to the sample intensity columns included in the output. This does not
#' affect the statistics, which are computed by \code{computeDiffStats}.
#' @param writeFile logical; whether to also write the results to a
#' tab-separated text file in the working directory, named after the contrast
#' (with " - " replaced by "_vs_") with the extension ".txt"
#' @return A \code{\link{data.frame}} with one row per feature, containing the
#' feature data columns (a "Count" column is renamed "Unique_peptides"; for
#' data summarised with \code{\link{summarizeIntensities}} it is the number of
#' distinct peptide sequences), one column of intensities per sample, and the
#' \code{topTable} results: \code{log2FC} (renamed from logFC),
#' \code{CI.L}, \code{CI.R}, \code{AvgIntensity} (renamed from AveExpr),
#' \code{t}, \code{P.Value}, \code{adj.P.Val} and \code{B}, plus
#' \code{controlLogFoldChange} if \code{controlGroup} is given. log2FC, t, B
#' and the sample intensities are rounded to 2 decimal places, and P.Value and
#' adj.P.Val to 2 significant digits.
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
#' diffexp <- getContrastResults(diffstats=diffstats, contrast="tam.24h_vs_vehicle")
#' 
#' @import limma
#' @importFrom Biobase exprs fData
#' @importFrom dplyr across arrange bind_cols desc mutate rename rename_with
#' select
#' @importFrom readr write_tsv
#' @importFrom stringr str_replace str_replace_all
#' @importFrom tidyselect all_of
#'
#' @export getContrastResults
getContrastResults <- function(diffstats, contrast, controlGroup = NULL, 
                               transform = TRUE, writeFile= FALSE) {
    checkArg_getContrastResults(diffstats, contrast, controlGroup, transform,
                                writeFile)
    
    message("Obtaining results for contrast  ", contrast, "\n")
    
    MSnSetObj <- diffstats$MSnSetObj
    fittedContrasts <- diffstats$fittedContrasts
    fittedLinearModel <- diffstats$fittedLM
    results <- topTable(fittedContrasts, 
                        coef = contrast, 
                        number = Inf, 
                        sort.by = "none",
                        confint=TRUE)
    # Get the groups involved in the contrast
    contrastCoeff <- fittedContrasts$contrasts[, contrast]
    contrastGroups <- names(contrastCoeff)[contrastCoeff != 0]

    # Calculate control log fold change if control group is provided
    if (!is.null(controlGroup)) {
        fittedIntensities <- as.data.frame(fittedLinearModel$coefficients)
        contrastIntensities <- select(fittedIntensities, all_of(contrastGroups))
        controlIntensity <- fittedIntensities[, controlGroup]
        results$controlLogFoldChange <- 
            apply(contrastIntensities - controlIntensity, 1, max)
    }
    
    intensities <- as.data.frame(exprs(MSnSetObj))
    if (transform) {
        intensities <- log2xplus1(intensities)
    }
    SamplesCol <- as.character(MSnSetObj$SampleName)
    results <- bind_cols(fData(MSnSetObj), intensities, results) |>
        arrange(desc(B)) |>
        mutate(across(c("logFC", "t", "B", all_of(SamplesCol)),
                      function(x) round(x, digits = 2))) |>
        mutate(across(c("P.Value", "adj.P.Val"),
                      function(x) signif(x, digits = 2))) |>
        # Count column may not be present
        rename_with(str_replace, 
                    pattern = "^Count$", 
                    replacement = "Unique_peptides") |>
        rename(AvgIntensity=AveExpr, log2FC=logFC)
    

    if (writeFile) {
        outFile <- contrast |> 
            str_replace(pattern = " - ", replacement = "_vs_") |> 
            str_c(".txt")
        write_tsv(results, outFile)
    }
    return(results)
}
