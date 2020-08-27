mergePeptides <- function(MSnSetObj, summarizationFunction, annotation, keepCols=NULL) {
  checkArg_mergePeptides(MSnSetObj, summarizationFunction, annotation, keepCols)
  
  concatUnique <- function(x){ unique(x) %>% str_c(collapse=";") }
  
  summarizedIntensities <- fData(MSnSetObj) %>%
    select(Accessions, Sequences, all_of(keepCols)) %>%
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
           everything())

  obj <- readMSnSet2(summarizedIntensities, ecol = sampleNames(MSnSetObj))
  pData(obj) <- pData(MSnSetObj)
  featureNames(obj) <- fData(obj)$phosseqid
  return(obj)
}

