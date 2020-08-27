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
    
    counts <- as.data.frame(fData(MSnSetObj)) %>%
        select(Accessions, Sequences) %>%
        distinct() %>%
        count(Accessions) %>%
        mutate(Accessions = as.character(Accessions)) %>%
        rename(Count = n)
    
    summarizedIntensities <- as.data.frame(exprs(MSnSetObj)) %>%
        mutate(Accessions = as.character(fData(MSnSetObj)$Accessions)) %>%
        group_by(Accessions) %>%
        summarize(across(everything(), summarizationFunction)) %>%
        left_join(counts, by = "Accessions") %>%
        left_join(annotation, by = "Accessions") %>%
        select(Accessions, colnames(annotation), Count, everything())
    
    expInd <- seq(
        grep("Count", colnames(summarizedIntensities)) + 1,
        ncol(summarizedIntensities)
    )
    
    obj <- readMSnSet2(summarizedIntensities, ecol = expInd)
    pData(obj) <- pData(MSnSetObj)
    featureNames(obj) <- fData(obj)$Accessions
    sampleNames(obj) <- pData(obj)$SampleName
    return(obj)
}

