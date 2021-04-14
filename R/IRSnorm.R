# Argument check function
checkArg_IRSnorm <- function(MSnSetObj, IRSname, groupingColumn){
  assert_that(is_MSnSet(MSnSetObj))
  assert_that(is_validIRSname(IRSname,MSnSetObj))
  assert_that(is.string(groupingColumn))
  assert_that(is_validMetadataColumn(groupingColumn, MSnSetObj))
}

#' Batch Correction by Internal Reference Scale
#' 
#' Performs batch correction on multiple runs using Internal Reference Scale
#' 
#' Internal Reference Scale (IRS) is pooled sample made up of aliquots of protein from all samples. The IRS is then run and measured
#' in each TMT experiment. The normalization procedure makes different measurements of the IRS exactly the same and puts all of the 
#' reporter ions on the same "intensity scale".  
#' The argument 'IRSname' is used to define the name of the Reference group within the SampleGroup column. 
#' The argument "groupingColumn" takes one of the column of pData(MSnSetObj) to define separate batches to correct. 
#' The default variable name is "Plex".   
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param IRSname character; name of Reference group within SampleGroup column
#' @param groupingColumn character; the pData(MSnSetObj) column name used to define batches; default="Plex"
#' @return An object of class \code{MSnSet} (see \code{\link{MSnSet-class}}) 
#' @examples
#' 
#' data(human_anno)
#' data(ER_ARID1A_KO_MCF7)
#' MSnset_SET1 <- convertToMSnset(ER_ARID1A_KO_MCF7$intensities_Set1,metadata=ER_ARID1A_KO_MCF7$metadata_Set1,
#'                 indExpData=c(7:15),Sequences=2,Accessions=6)
#' MSnset_SET2 <- convertToMSnset(ER_ARID1A_KO_MCF7$intensities_Set2,metadata=ER_ARID1A_KO_MCF7$metadata_Set2,
#'                indExpData=c(7:15),Sequences=2,Accessions=6)
#' MSnset_SET3 <- convertToMSnset(ER_ARID1A_KO_MCF7$intensities_Set3,metadata=ER_ARID1A_KO_MCF7$metadata_Set3,
#'                 indExpData=c(7:15),Sequences=2,Accessions=6)
#' MSnset_SET4 <- convertToMSnset(ER_ARID1A_KO_MCF7$intensities_Set4,metadata=ER_ARID1A_KO_MCF7$metadata_Set4,
#'                 indExpData=c(7:14),Sequences=2,Accessions=6)
#' MSnset_SET5 <- convertToMSnset(ER_ARID1A_KO_MCF7$intensities_Set5,metadata=ER_ARID1A_KO_MCF7$metadata_Set5,
#'                 indExpData=c(7:15),Sequences=2,Accessions=6)
#' MSnset_SET1_norm <- normalizeScaling(MSnset_SET1, median)
#' MSnset_SET2_norm <- normalizeScaling(MSnset_SET2, median)
#' MSnset_SET3_norm <- normalizeScaling(MSnset_SET3, median)
#' MSnset_SET4_norm <- normalizeScaling(MSnset_SET4, median)
#' MSnset_SET5_norm <- normalizeScaling(MSnset_SET5, median)
#' MSnset_SET1_Pnorm <- summarizeIntensities(MSnset_SET1_norm, sum, human_anno)
#' MSnset_SET2_Pnorm <- summarizeIntensities(MSnset_SET2_norm, sum, human_anno)
#' MSnset_SET3_Pnorm <- summarizeIntensities(MSnset_SET3_norm, sum, human_anno)
#' MSnset_SET4_Pnorm <- summarizeIntensities(MSnset_SET4_norm, sum, human_anno)
#' MSnset_SET5_Pnorm <- summarizeIntensities(MSnset_SET5_norm, sum, human_anno)
#' MSnset_SET1_Pnorm <- updateSampleNames(updateFvarLabels(MSnset_SET1_Pnorm))
#' MSnset_SET2_Pnorm <- updateSampleNames(updateFvarLabels(MSnset_SET2_Pnorm))
#' MSnset_SET3_Pnorm <- updateSampleNames(updateFvarLabels(MSnset_SET3_Pnorm))
#' MSnset_SET4_Pnorm <- updateSampleNames(updateFvarLabels(MSnset_SET4_Pnorm))
#' MSnset_SET5_Pnorm <- updateSampleNames(updateFvarLabels(MSnset_SET5_Pnorm))
#' MSnset_comb <- MSnbase::combine(MSnset_SET1_Pnorm,MSnset_SET2_Pnorm,MSnset_SET3_Pnorm,MSnset_SET4_Pnorm,MSnset_SET5_Pnorm)
#' tokeep <- which(complete.cases(fData(MSnset_comb))==TRUE)
#' MSnset_comb <- MSnset_comb[tokeep,]
#' sampleNames(MSnset_comb) <- pData(MSnset_comb)$SampleName
#' fData(MSnset_comb) <- fData(MSnset_comb)[,c(2,3,6)]
#' colnames(fData(MSnset_comb)) <- c("Sequences","Modifications","Accessions")
#' MSnset_comb_corr <- IRSnorm(MSnset_comb,IRSname="Ref",groupingColumn="Run")

#' @export IRSnorm


IRSnorm <- function(MSnSetObj, IRSname="RefPool", 
                         groupingColumn="Plex") {
  checkArg_IRSnorm(MSnSetObj, IRSname, groupingColumn)
  
  Ref_Set <- MSnSetObj[,which(pData(MSnSetObj)$SampleGroup==IRSname)]
  plex_id <- unique(pData(MSnSetObj)[,groupingColumn])
  allgrps <- split(Ref_Set,groupingColumn)
  grpMean <-  BiocGenerics::lapply(allgrps, function(x) {apply(exprs(x),1,mean)})
  irs_data <- matrix(unlist(grpMean), ncol = length(grpMean), byrow = FALSE)
  irs_data_geomean <- apply(irs_data, 1, function(x) exp(mean(log(x))))
  irs_factors <- irs_data_geomean/irs_data

  for (i in 1:length(plex_id))
  {
    exprs(MSnSetObj)[,which(pData(MSnSetObj)[,groupingColumn]==i)] <- 
      exprs(MSnSetObj)[,which(pData(MSnSetObj)[,groupingColumn]==i)] * irs_factors[,i]
  }
  return(MSnSetObj)
}
