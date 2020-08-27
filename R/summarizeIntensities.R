##### Summarization function ##########

# Functions for summarizing intensity data

# Summarizes multiple peptide measurements for a protein.
# Assumes that there are columns for each of the samples specified and uses
# Accessions column for grouping peptide-level measurements.
# Filters any rows with missing values.
# Typical summarization functions are sum, mean and median.
# For successful running of this function the annotation file must have at least
# two columns: "Accessions" and "GeneSymbol"
# The MSnSetObj must have column "Sequences" denoting its a peptide dataset

summarizeIntensities <- function(MSnSetObj, summarizationFunction, annotation) {
    checkArg_summarizeIntensities(MSnSetObj, summarizationFunction, annotation)
    
   summarizedIntensities <- fData(MSnSetObj) %>%
       select(Accessions, Sequences) %>%
       mutate(across(c(Accessions, Sequences), as.character)) %>% 
       bind_cols(as.data.frame(exprs(MSnSetObj))) %>% 
       group_by(Accessions) %>%
       summarize(across(where(is.numeric), summarizationFunction), 
                 Count=n_distinct(Sequences)) %>% 
       left_join(annotation, by = "Accessions") %>%
       select(Accessions, colnames(annotation), Count, everything())
    
    obj <- readMSnSet2(summarizedIntensities, ecol = sampleNames(MSnSetObj))
    pData(obj) <- pData(MSnSetObj)
    featureNames(obj) <- fData(obj)$Accessions
    sampleNames(obj) <- pData(obj)$SampleName
    return(obj)
}

