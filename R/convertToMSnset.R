convertToMSnset <- function(ExpObj, metadata, indExpData, Sequences=NULL, 
                            Accessions, type="peptide", rmMissing=TRUE) {
    checkArg_convertToMSnset(ExpObj, metadata, indExpData, Sequences, 
                             Accessions, type, rmMissing)
    
    colnames(ExpObj)[Accessions] <- "Accessions"
    colnames(ExpObj)[Sequences] <- "Sequences"
    
    if (rmMissing) {
        ExpObj %<>% filter(across(indExpData, ~!is.na(.x)))
    }
    obj <- readMSnSet2(ExpObj, ecol = indExpData)

    # if the meta data is a tibble, we can't set rownames
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- as.character(metadata$SampleName)
    ind <- match(sampleNames(obj), rownames(metadata))
    metadata <- metadata[ind, ]
    pData(obj) <- metadata
    if (type == "protein") {
        featureNames(obj) <- fData(obj)$Accessions
    } else {
        featureNames(obj) <- paste0("peptide_", featureNames(obj))
    }
    return(obj)
}

