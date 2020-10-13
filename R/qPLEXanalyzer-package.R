# Argument check function


#' exp2_Xlink dataset
#' 
#' An ER qPLEX-RIME experiment was performed to compare two different methods
#' of crosslinking.  MCF7 cells were double crosslinked with DSG/formaldehyde
#' (double) or with formaldehyde alone (single). Four biological replicates
#' were obtained for each condition along with two IgG pooled samples from each
#' replicate.
#' 
#' 
#' @name exp2_Xlink
#' @docType data
#' @format An object of class \code{\link{list}} related to peptides
#' quantification. It consists of qPLEX-RIME data of 10 samples divided into
#' three conditions (FA, DSG.FA and IgG).
#' @return An object of class \code{\link{list}} related to peptides
#' quantification.
#' @keywords data datasets
NULL





#' exp3_OHT_ESR1 dataset
#' 
#' Three ER qPLEX-RIME (10plex) experiments were performed to investigate the
#' dynamics of the ER complex assembly upon 4-hydrotamoxifen (OHT) treatment at
#' 2h, 6h and 24h or at 24h post-treatment with the drug-vehicle alone
#' (ethanol). Two biological replicates of each condition were included in each
#' experiment to finally consider a total of six replicates per time point.
#' Additionally, MCF7 cells were treated with OHT or ethanol and cross-linked
#' at 24h post-treatment in each experiment to be used for mock IgG pull-downs
#' and to enable discrimination of non-specific binding in the same experiment.
#' This is a timecourse experiment to study the effect of tamoxifen in ER
#' interactome using qPLEX-RIME method.
#' 
#' 
#' @name exp3_OHT_ESR1
#' @docType data
#' @format An object of class \code{\link{list}} related to peptides
#' quantification. It consists of qPLEX-RIME data from three experimental runs.
#' Each run contains 10 samples divided into five conditions (IgG, vehicle,
#' tam.2h, tam.6h and tam.24h).
#' @return An object of class \code{\link{list}} related to peptides
#' quantification.
#' @keywords data datasets
NULL





#' human_anno dataset
#' 
#' Uniprot Human protein annotation table.
#' 
#' 
#' @name human_anno
#' @docType data
#' @format An object of class \code{\link{data.frame}} consisting of uniprot
#' human protein annotation.
#' @keywords data datasets
NULL





#' Tools for qPLEX-RIME data analysis
#' 
#' Tools for quantitiative proteomics data analysis generated from qPLEX-RIME
#' method
#' The package offers the following functionalities
#' Data processing, normalization & analysis: 
#' \itemize{ 
#' \item \code{convertToMSnset}: Converts quantitative data to a MSnSet 
#' \item \code{summarizeIntensities}: Summarizes multiple peptide measurements 
#' for a protein 
#' \item \code{normalizeQuantiles}: Performs quantile normalization on the
#' peptides/proteins intensities 
#' \item \code{normalizeScaling}: Performs scaling normalization on the
#' peptides/proteins intensities (mean, median or sum) 
#' \item \code{groupScaling}: Performs scaling normalization on the 
#' peptides/proteins intensities within group (median or mean) 
#' \item \code{rowScaling}: Normalization by scaling peptide/protein intensity
#' across all samples 
#' \item \code{regressIntensity}: Performs linear regression on protein 
#' intensities based on selected protein 
#' \item \code{computeDiffStats}: Compute differential statistics for the given
#' contrasts 
#' \item \code{getContrastResults}: Get differential statistics results for 
#' given contrast 
#' }
#' Visualization: 
#' \itemize{ 
#' \item \code{assignColours}: Assigns colours to samples in groups 
#' \item \code{corrPlot}: Correlation plot of all the samples 
#' \item \code{coveragePlot}: Computes and display protein sequence coverage of
#' \item \code{hierarchicalPlot}: Hierarchical clustering plot of all the 
#' samples 
#' \item \code{intensityBoxplot}: Intensity distribution boxplot of all the 
#' samples
#' \item \code{intensityPlot}: Intensity distribution plot of all the samples
#' \item \code{maVolPlot}: MA or Volcano plot of differential analysis results 
#' \item \code{pcaPlot}: PCA plot of all the samples 
#' \item \code{peptideIntensityPlot}: Peptide intensity distribution plot of 
#' specific protein 
#' \item \code{plotMeanVar}: Computes and plots mean-variance for samples in 
#' MSnSet 
#' \item \code{rliPlot}: Relative intensity plot of all the samples 
#' selected protein in proteomics experiment 
#' }
#' 
#' @name qPLEXanalyzer-package
#' @aliases qPLEXanalyzer-package qPLEXanalyzer
#' @docType package
#' @author Matthew Eldridge, Kamal Kishore, Ashley Sawle (Maintainer)
#' 
#' \email{ads2202cu@@gmail.com}
#' @keywords package
NULL



