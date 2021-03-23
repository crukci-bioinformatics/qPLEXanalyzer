# Argument check function
checkArg_IRSnorm <- function(MSnSetObj, IRSname, groupingColumn){
  assert_that(is_MSnSet(MSnSetObj))
  assert_that(is_validIRSname(IRSname,MSnSetObj))
  assert_that(is_validGroupingColumn(groupingColumn, MSnSetObj))
}

# Perform IRS normalization

IRSnorm <- function(MSnSetObj, IRSname="RefPool", 
                         groupingColumn="Plex") {
  checkArg_IRSnorm(MSnSetObj, IRSname, groupingColumn)
  
  Ref_Set <- MSnSetObj[,which(pData(MSnSetObj)$SampleGroup==IRSname)]
  plex_id <- unique(pData(MSnSetObj)[,groupingColumn])
  allgrps <- split(Ref_Set,groupingColumn)
  grpMean <- lapply(allgrps, function(x) {apply(exprs(x),1,mean)})
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
