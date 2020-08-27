mergePeptides <- function(MSnSetObj, summarizationFunction, annotation, PosMasterProt=NULL) {
  checkArg_mergePeptides(MSnSetObj, summarizationFunction, annotation, PosMasterProt=NULL)
  
  if(!is.null(PosMasterProt)){
    colnames(fData(MSnSetObj))[PosMasterProt] <- "PosInMasterProtein"
    counts <- fData(MSnSetObj) %>%
      select(Accessions, Sequences, PosInMasterProtein) %>%
      mutate(phosseqid = str_c(Sequences, "_", Accessions)) %>%
      count(phosseqid, PosInMasterProtein, name = "Count")
  }
  else
  {
    counts <- fData(MSnSetObj) %>%
      select(Accessions, Sequences) %>%
      mutate(phosseqid = str_c(Sequences, "_", Accessions)) %>%
      count(phosseqid, name = "Count")
  }
  
  summarizedIntensities <- as.data.frame(exprs(MSnSetObj)) %>%
    mutate(phosseqid = str_c(fData(MSnSetObj)$Sequences,
                             "_",
                             fData(MSnSetObj)$Accessions)) %>%
    group_by(phosseqid) %>%
    summarize(across(everything(), summarizationFunction)) %>%
    left_join(counts, by = "phosseqid") 
  
  summarizedIntensities$Accessions <- unlist(lapply(strsplit(summarizedIntensities$phosseqid,split="_"),
                                                    function(x) x[2]))
  
  summarizedIntensities <- summarizedIntensities %>%
    left_join(annotation, by = "Accessions") %>%
    select(Accessions, colnames(annotation), Count, contains("PosInMasterProtein"), everything())
  
  if(!is.null(PosMasterProt)){
    expInd <- seq(
      grep("Count", colnames(summarizedIntensities)) + 3,
      ncol(summarizedIntensities)
    )
  }
  else
  {
    expInd <- seq(
      grep("Count", colnames(summarizedIntensities)) + 2,
      ncol(summarizedIntensities)
    )
  }
  
  obj <- readMSnSet2(summarizedIntensities, ecol = expInd)
  pData(obj) <- pData(MSnSetObj)
  featureNames(obj) <- fData(obj)$phosseqid
  sampleNames(obj) <- pData(obj)$SampleName
  return(obj)
}

