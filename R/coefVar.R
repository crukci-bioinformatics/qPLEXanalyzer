# Argument check function
checkArg_coefVar <- function(MSnSetObj){
  assert_that(is_MSnSet(MSnSetObj))
}

#' Calculating the coefficient of variation by utilizing expression data within individual sample groups.
#' 
#' Calculating the coefficient of variation by utilizing peptide/protein expression data within individual sample groups.
#' 
#' In this approach, we calculate distributions of the coefficient of variation (CV) 
#' for the dataset. The CVs are determined based on peptides or proteins intensities 
#' within each sample group, and the results are visualized through boxplots and 
#' cumulative fraction plots for each sample group. 
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @return An object of class \code{list} consisting of object of class \code{ggplot}) 
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' res <- coefVar(MSnSet_data)                             
#' @import ggplot2
#' @importFrom Biobase exprs exprs<- pData
#' @importFrom dplyr arrange group_by mutate ungroup
#' @importFrom magrittr %>%
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
                          values_to = "CV") %>% 
    group_by(SampleName) %>% 
    arrange(CV) %>% 
    mutate(cum_sum=cumsum(CV),fraction=cum_sum/sum(CV)) %>% 
    ungroup()
  
  ### boxplot for each samplegroup
  boxpl <- ggplot(CV_comb,mapping = aes(x=SampleName,y=CV)) +
    geom_boxplot() + 
    theme_bw() + 
    labs(title = "Boxplot of Coefficient of Variation",
         x = "Sample Name",
         y = "Coefficient of Variation (%)")
  
  ### cumulative fraction plot for each samplegroup
  cumulpl <- ggplot(CV_comb,mapping = aes(x=CV,y=fraction,color=SampleName)) +
    geom_line() + 
    theme_bw() +
    labs(title = "Cumulative Fraction Plot",
         x = "Coefficient of Variation (%)",
         y = "Cumulative Fraction")
  
  ## return both plots in the form of list
  return(list(boxpl=boxpl,cumulpl=cumulpl))
}
