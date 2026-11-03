# Argument check function
checkArg_computeDiffStats <- function(MSnSetObj, 
                                      batchEffect, 
                                      transform, 
                                      contrasts, 
                                      trend, 
                                      robust){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validBatchEffect(batchEffect, MSnSetObj))
    assert_that(is.flag(transform))
    assert_that(is.character(contrasts))
    assert_that(is.flag(trend))
    assert_that(is.flag(robust))
}

# Log2 with addition of 1e-10 count to deal with zeros
log2xplus1 <- function(x) {
    log2(x + 1e-10)
}

# Rewrite user contrasts so that they refer to the design matrix column names.
# `levelNames` are the design column names before spaces are replaced by "_".
# - Group names containing spaces are replaced by their "_" version (longest
#   names first, so a name contained in a longer one is not split).
# - A contrast that is exactly "<level>_vs_<level>" is read as
#   "<level> - <level>", unless the whole string is itself a level name.
# Everything else in the contrast is left untouched, so contrasts with
# brackets, several terms or group names containing "_vs_" are kept as given.
formatContrasts <- function(contrasts, levelNames) {
    levelsFixed <- gsub(" ", "_", levelNames)
    spaced <- levelNames[grepl(" ", levelNames)]
    spaced <- spaced[order(nchar(spaced), decreasing = TRUE)]
    for (lev in spaced) {
        contrasts <- gsub(lev, gsub(" ", "_", lev), contrasts, fixed = TRUE)
    }
    for (i in seq_along(contrasts)) {
        parts <- trimws(strsplit(contrasts[i], "_vs_", fixed = TRUE)[[1]])
        if (!contrasts[i] %in% levelsFixed &&
            length(parts) == 2 && all(parts %in% levelsFixed)) {
            contrasts[i] <- paste(parts, collapse = " - ")
        }
    }
    contrasts
}

#' Compute differential statistics
#' 
#' Compute differential statistics on the given contrasts, based on
#' \code{\link{limma}} functions.
#' 
#' A statistical analysis for the identification of differentially regulated or
#' bound proteins is carried out using limma based analysis. It uses linear
#' models to assess differential expression in the context of multifactor
#' designed experiments.  Firstly, a linear model is fitted for each protein
#' where the model includes a variable for each sample group ("SampleGroup")
#' plus any additional variables given in \code{batchEffect}. Then, log2 fold
#' changes between comparisons are estimated. Multiple testing correction of
#' p-values are applied using the Benjamini-Hochberg method to control the
#' false discovery rate (FDR).
#'
#' In order to correct for batch effects (e.g. MS run), additional variable(s)
#' can be given in \code{batchEffect}. Each must correspond to a column name in
#' pData(MSnSetObj). "SampleGroup" is always included in the model; the
#' variables in \code{batchEffect} are added to it.
#'
#' Repeated measurements of the same biological unit, such as technical
#' replicates of a sample or several experiments performed on the same cell
#' line, are accounted for through the "TechRep" column of pData(MSnSetObj).
#' If this column contains any non-NA values, the correlation between
#' measurements of the same unit is estimated with
#' \code{\link[limma]{duplicateCorrelation}} using "TechRep" as the blocking
#' variable, and the linear model is then fitted with
#' \code{\link[limma]{lmFit}} using the same block and the consensus
#' correlation. Samples sharing a "TechRep" value are treated as correlated
#' measurements of one biological unit, so all measurements of the same unit
#' must have the same value, and different units must have different values.
#' A unit may be measured in more than one group (for example the same subject
#' measured under two conditions, as in the "Multi-level experiments" section
#' of the limma User's Guide); what must be avoided is reusing a replicate or
#' cell line number for different units. For example, with two cell lines per
#' group, four experiments per cell line and each cell line belonging to one
#' group only, label the samples by cell line (e.g. "Control_1", "Control_2",
#' "POLG_1", "POLG_2") rather than by a cell line number that repeats across
#' groups (e.g. 1, 2). If "TechRep" is used, every sample must have a value: a
#' sample without repeated measurements should be given its own unique value.
#' If all "TechRep" values are NA, the linear model is fitted without blocking.
#' A variable nested within "SampleGroup", such as cell line in the example
#' above, should not be given in \code{batchEffect}, because its effects
#' cannot be estimated separately from the groups.
#'
#' @param MSnSetObj MSnSet; An object of class MSnSet
#' @param batchEffect character; vector of variable(s) to correct for batch
#' effect, Default : "SampleGroup". All variables in the model, including
#' "SampleGroup", are treated as categorical: numeric or logical columns
#' (e.g. a run number stored as an integer) are converted to factors.
#' @param transform logical; apply log2 transformation to the raw intensities
#' before fitting. The transformation is \code{log2(x + 1e-10)}, so zero
#' intensities become about -33; set to FALSE if the intensities are already
#' on a log scale.
#' @param contrasts character; named character vector of contrasts for
#' differential statistics, written in terms of the SampleGroup names
#' (e.g. "tam.24h - vehicle"). Group names containing spaces can be written
#' as they are or with the spaces replaced by "_". A contrast written as
#' "groupA_vs_groupB", where both are existing group names, is read as
#' "groupA - groupB". Any other expression accepted by
#' \code{\link[limma]{makeContrasts}} can be used, e.g.
#' "(A - B) - (C - D)".
#' @param trend logical; passed to \code{\link[limma]{eBayes}}: should an
#' intensity-dependent trend be allowed for the prior variance? If FALSE the
#' prior variance is constant.
#' @param robust logical; passed to \code{\link[limma]{eBayes}}: should the
#' estimation of the prior degrees of freedom and prior variance be
#' robustified against outlier sample variances?
#' @return A list object containing three components: MSnSetObj of class
#' \code{MSnSet} (see \code{\link{MSnSet-class}}) object), fittedLM (fitted
#' linear model) and fittedContrasts. This object should be input into
#' getContrastResults function to get differential results.  See
#' \code{\link{eBayes}} function of \code{\link{limma}} for more details on
#' differential statistics.
#'
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
#' contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle", 
#'                tam.6h_vs_vehicle = "tam.6h - vehicle")
#' diffstats <- computeDiffStats(MSnSetObj=MSnset_Pnorm, contrasts=contrasts)
#'
#' @import limma
#' @importFrom Biobase exprs pData
#' @importFrom stats as.formula model.matrix
#' 
#' @export computeDiffStats
computeDiffStats <- function(MSnSetObj, batchEffect = NULL, transform = TRUE, 
                             contrasts, trend = TRUE, robust = TRUE) {
    checkArg_computeDiffStats(MSnSetObj, batchEffect, transform, contrasts,
                              trend, robust)
    
    message("Fitting linear model")
    intensities <- as.data.frame(exprs(MSnSetObj))
    if (transform) {
        intensities <- log2xplus1(intensities)
    }
    batchEffect <- unique(c("SampleGroup", batchEffect))
    model <- as.formula(paste(c("~0", batchEffect), collapse = " + "))
    # model variables are categorical: numeric or logical columns (e.g. BioRep
    # or Run stored as integers) would otherwise be fitted as a single slope
    modelData <- pData(MSnSetObj)
    numericVars <- batchEffect[vapply(modelData[batchEffect], function(x) {
        is.numeric(x) || is.logical(x)
    }, logical(1))]
    if (length(numericVars) > 0) {
        message("Treating numeric column(s) ",
                paste(numericVars, collapse = ", "),
                " as categorical (factor) in the linear model")
        modelData[numericVars] <- lapply(modelData[numericVars], factor)
    }
    design <- model.matrix(model, data = modelData)
    levelNames <- sub("^SampleGroup", "", colnames(design))
    colnames(design) <- gsub(" ", "_", levelNames)
    if (length(which(is.na(pData(MSnSetObj)$TechRep) == FALSE)) > 0) {
        techRepBlock <- pData(MSnSetObj)$TechRep
        dupcor <- duplicateCorrelation(intensities,
                                       design = design,
                                       block = techRepBlock)
        fit <- lmFit(as.matrix(intensities),
                     design = design,
                     weights = NULL,
                     block = techRepBlock,
                     correlation = dupcor$consensus)
    } else {
        fit <- lmFit(as.matrix(intensities), design = design, weights = NULL)
    }
    
    message("Fitting contrasts\n")
    
    contrasts <- formatContrasts(contrasts, levelNames)

    contrasts <- makeContrasts(contrasts = contrasts, levels = fit$design)
    contrastsfit <- contrasts.fit(fit, contrasts)
    
    message("Computing empirical Bayes statistics for differential expression")
    fittedContrasts <- eBayes(contrastsfit, trend = trend, robust = robust)
    return(diffstats <- list(MSnSetObj = MSnSetObj, 
                             fittedLM = fit, 
                             fittedContrasts = fittedContrasts))
}

