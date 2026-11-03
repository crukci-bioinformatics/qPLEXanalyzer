# Argument check function
checkArg_coefVar <- function(MSnSetObj){
  assert_that(is_MSnSet(MSnSetObj))
}

#' Calculating the coefficient of variation by utilizing expression data within individual sample groups
#'
#' Calculating the coefficient of variation by utilizing peptide/protein expression data within individual sample groups.
#' 
#' In this approach, we calculate distributions of the coefficient of variation (CV)
#' for the dataset. The CVs are determined based on peptides or proteins intensities
#' within each sample group, and the results are visualized through boxplots and
#' cumulative fraction plots for each sample group.
#'
#' The samples are grouped by the "SampleGroup" column of pData(MSnSetObj). For
#' each feature and group the CV is computed as 100 * standard deviation / mean
#' of the intensities as given (no transformation is applied, so raw intensities
#' should be used). The cumulative fraction plot shows, for each group, the
#' proportion of features with a CV less than or equal to the value on the
#' x-axis. CVs that cannot be computed (e.g. features with missing values, or
#' with a mean of zero) are excluded from the cumulative fraction plot, with a
#' message giving their number.
#'
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @return A \code{list} of two objects of class \code{ggplot}: \code{boxpl},
#' boxplots of the CVs of each sample group, and \code{cumulpl}, the cumulative
#' fraction plot.
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' # the IgG samples contain zero intensities, which give a mean of zero (and an
#' # undefined CV) for some peptides; a constant is added to avoid this
#' exprs(MSnSet_data) <- exprs(MSnSet_data)+1.1
#' res <- coefVar(MSnSet_data)
#' @import ggplot2
#' @importFrom Biobase exprs exprs<- pData
#' @importFrom dplyr arrange filter group_by mutate n row_number ungroup
#' @importFrom tidyr pivot_longer
#'
#' @export coefVar

coefVar <- function(MSnSetObj) {
  
  ## internal function to compute coefficient of variation
  computeCV <- function(x)
  {
    intens <- exprs(x)
    av <- rowMeans(intens)
    sd <- apply(intens,1,sd)
    cv <- 100 * sd / av
    return(cv)
  }
  
  ### split the object by different samplegroup
  allgrps <- split(MSnSetObj,"SampleGroup")
  
  ### compute CVs for each samplegroup
  res_CV <- lapply(allgrps,computeCV)
  res_CV_comb <- bind_cols(res_CV)
  
  ### reshape data for visualisation
  CV_comb <- pivot_longer(res_CV_comb, everything(),names_to = "SampleName",
                          values_to = "CV")

  ### cumulative fraction: proportion of features with CV <= x in each group;
  ### CVs that are NA/NaN (missing values or zero mean) cannot be ranked
  nNonFinite <- sum(!is.finite(CV_comb$CV))
  if (nNonFinite > 0) {
    message(nNonFinite, " CV value(s) are missing or not finite (e.g. missing ",
            "intensities or zero mean) and are excluded from the cumulative ",
            "fraction plot.")
  }
  CV_cumul <- CV_comb |>
    filter(is.finite(CV)) |>
    arrange(SampleName, CV) |>
    group_by(SampleName) |>
    mutate(fraction = row_number() / n()) |>
    ungroup()

  ### boxplot for each samplegroup
  boxpl <- ggplot(CV_comb,mapping = aes(x=SampleName,y=CV)) +
    geom_boxplot() + 
    theme_bw() + 
    labs(title = "Boxplot of Coefficient of Variation",
         x = "Sample Name",
         y = "Coefficient of Variation (%)") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 10)) 
    
  
  ### cumulative fraction plot for each samplegroup
  cumulpl <- ggplot(CV_cumul,mapping = aes(x=CV,y=fraction,color=SampleName)) +
    geom_line() + 
    theme_bw() +
    labs(title = "Cumulative Fraction Plot",
         x = "Coefficient of Variation (%)",
         y = "Cumulative Fraction") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 10))
  
  ## return both plots in the form of list
  return(list(boxpl=boxpl,cumulpl=cumulpl))
}
