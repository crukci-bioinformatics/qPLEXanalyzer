# Fits a linear model to the intensity data using limma.
# The PhenoData table must contain SampleName and SampleGroup columns.
# The intensities table must contain a column for each sample in PhenoData
# Log2 with addition of 1e-10 count to deal with zeros
log2xplus1 <- function(x) {
    log2(x + 1e-10)
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
#' where the model includes variables for each group and MS run. Then, log2
#' fold changes between comparisions are estimated. Multiple testing correction
#' of p-values are applied using the Benjamini-Hochberg method to control the
#' false discovery rate (FDR).
#' 
#' In order to correct for batch effect, variable(s) can be defined. It should
#' corresponds to a column name in pData(MSnSetObj). The default variable is
#' "SampleGroup" that distinguish between two groups.  If more variables are
#' defined they are added to default.
#' 
#' @param MSnSetObj MSnSet; An object of class MSnSet
#' @param batchEffect character; vector of variable(s) to correct for batch
#' effect, Default : "SampleGroup"
#' @param transform logical; apply log2 transformation to the raw intensitites
#' @param contrasts character; named character vector of contrasts for
#' differential statistics
#' @param trend logical; TRUE or FALSE
#' @param robust logical; TRUE or FALSE
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
#' metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#' indExpData=c(7:16), Sequences=2, Accessions=6)
#' MSnset_norm <- groupScaling(MSnSet_data, scalingFunction=median)
#' MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
#' contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle", tam.6h_vs_vehicle = "tam.6h - vehicle")
#' diffstats <- computeDiffStats(MSnSetObj=MSnset_Pnorm, contrasts=contrasts)
#'
#' @import limma
#' @importFrom Biobase exprs pData
#' @importFrom magrittr %>%
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
    design <- model.matrix(model, data = pData(MSnSetObj))
    colnames(design) <- colnames(design) %>%
        sub(pattern = "^SampleGroup", replacement = "") %>%
        gsub(pattern = " ", replacement = "_")
    if (length(which(is.na(pData(MSnSetObj)$TechRep) == FALSE)) > 0) {
        dupcor <- duplicateCorrelation(intensities, 
                                       design = design, 
                                       block = pData(MSnSetObj)$TechRep)
        fit <- lmFit(intensities, 
                     design = design, 
                     weights = NULL, 
                     correlation = dupcor$consensus)
    } else {
        fit <- lmFit(intensities, design = design, weights = NULL)
    }
    
    message("Fitting contrasts\n")
    
    contrasts <- contrasts %>%
        gsub(pattern = " ", replacement = "_") %>%
        sub(pattern = "_-_", replacement = " - ") %>%
        sub(pattern = "_vs_", replacement = " - ")
    
    contrasts <- makeContrasts(contrasts = contrasts, levels = fit$design)
    contrastsfit <- contrasts.fit(fit, contrasts)
    
    message("Computing empirical Bayes statistics for differential expression")
    fittedContrasts <- eBayes(contrastsfit, trend = trend, robust = robust)
    return(diffstats <- list(MSnSetObj = MSnSetObj, 
                             fittedLM = fit, 
                             fittedContrasts = fittedContrasts))
}

