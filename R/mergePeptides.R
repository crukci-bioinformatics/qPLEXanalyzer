mergePeptides <- function(MSnSetObj, summarizationFunction, annotation, PosMasterProt=NULL) {
  checkArg_mergePeptides(MSnSetObj, summarizationFunction, annotation, PosMasterProt)
  
  concatUnique <- function(x){ unique(x) %>% str_c(collapse=";") }
  
  summarizedIntensities <- fData(MSnSetObj) %>%
    select(Accessions, Sequences, PosInMasterProtein=PosMasterProt) %>%
    mutate(across(everything(), as.character)) %>% 
    mutate(phosseqid = str_c(Sequences, "_", Accessions)) %>%  
    bind_cols(as.data.frame(exprs(MSnSetObj))) %>% 
    group_by(phosseqid, Accessions) %>%
    select(-Sequences) %>% 
    summarize(across(where(is.character), concatUnique),
              across(where(is.numeric), summarizationFunction), 
              Count=n())   %>% 
    ungroup() %>% 
    left_join(annotation, by = "Accessions") %>%
    select(Accessions, 
           colnames(annotation), 
           Count, 
           contains("PosInMasterProtein"), 
           everything())

  obj <- readMSnSet2(summarizedIntensities, ecol = sampleNames(MSnSetObj))
  pData(obj) <- pData(MSnSetObj)
  featureNames(obj) <- fData(obj)$phosseqid
  sampleNames(obj) <- pData(obj)$SampleName
  return(obj)
}

