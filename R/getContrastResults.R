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
        mutate(across(c("logFC", "t", "B", SamplesCol), round, digits = 2)) %>%
        mutate(across(c("P.Value", "adj.P.Val"), signif, digits = 2)) %>% 
        rename_with(str_replace, 
                    pattern = "^Count$", 
                    replacement = "Unique_peptides") %>%
        rename(AvgIntensity=AveExpr, log2FC=logFC)
    
    if (writeFile == TRUE) {
        write.table(results, paste0(names(contrast), ".txt"), 
                    quote = FALSE, 
                    sep = "\t", 
                    row.names = FALSE)
    }
    return(results)
}
