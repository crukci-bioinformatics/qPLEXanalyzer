convertToMSnset <- function(ExpObj, metadata, indExpData, Sequences=NULL, 
                            Accessions, type="peptide", rmMissing=TRUE) {
    checkArg_convertToMSnset(ExpObj, metadata, indExpData, Sequences, 
                             Accessions, type, rmMissing)
    
    colnames(ExpObj)[Accessions] <- "Accessions"
    colnames(ExpObj)[Sequences] <- "Sequences"
    
    if (rmMissing) {
        ExpObj <- ExpObj %>% 
            drop_na(indExpData)
    }
    obj <- readMSnSet2(ExpObj, ecol = indExpData)

    pData(obj) <- metadata %>% 
        mutate(rowname = SampleName) %>% 
        column_to_rownames() %>% 
        `[`(sampleNames(obj), TRUE)
            
    if (type == "protein") {
        featureNames(obj) <- fData(obj)$Accessions
    } else {
        featureNames(obj) <- str_c("peptide_", featureNames(obj))
    }
    return(obj)
}

