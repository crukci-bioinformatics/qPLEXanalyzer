getContrastResults <- function(diffstats, contrast, controlGroup = NULL, 
                               transform = TRUE, writeFile= FALSE) {
    checkArg_getContrastResults(diffstats, contrast, controlGroup, transform,
                                writeFile)
    
    message("Obtaining results for contrast", contrast, "\n")
    contrast <- contrast %>%
        gsub(pattern = " ", replacement = "_") %>%
        sub(pattern = "_-_", replacement = " - ") %>%
        sub(pattern = "_vs_", replacement = " - ")
    
    MSnSetObj <- diffstats$MSnSetObj
    fittedContrasts <- diffstats$fittedContrasts
    fittedLinearModel <- diffstats$fittedLM
    results <- topTable(fittedContrasts, 
                        coef = contrast, 
                        number = Inf, 
                        sort.by = "none",confint=TRUE)
    contrastGroups <- contrast %>% strsplit(" - ") %>% unlist()
    fittedIntensities <- as.data.frame(fittedLinearModel$coefficients)
    contrastIntensities <- select(fittedIntensities, one_of(contrastGroups))
    
    if (!is.null(controlGroup)) {
        controlIntensity <- fittedIntensities[, controlGroup]
        results$controlLogFoldChange <- 
            apply(contrastIntensities - controlIntensity, 1, max)
    }
    
    intensities <- as.data.frame(exprs(MSnSetObj))
    if (transform) {
        intensities <- log2xplus1(intensities)
    }
    SamplesCol <- as.character(MSnSetObj$SampleName)
    results <- cbind(fData(MSnSetObj), intensities, results)
    results <- results %>%
        arrange(desc(B))
    results <- results %>%
        mutate_at(funs(round(., digits = 2)), 
                  .vars = c("logFC", "t", "B", SamplesCol)) %>%
        mutate_at(funs(signif(., digits = 2)), 
                  .vars = c("P.Value", "adj.P.Val")) %>% 
        rename_all(function(x){str_replace(x, 
                                           "^Count$", 
                                           "Unique_peptides")}) %>%
        rename(AvgIntensity=AveExpr, log2FC=logFC)
    
    if (writeFile == TRUE) {
        write.table(results, paste0(names(contrast), ".txt"), 
                    quote = FALSE, 
                    sep = "\t", 
                    row.names = FALSE)
    }
    return(results)
}
