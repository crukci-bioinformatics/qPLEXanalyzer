\name{qPLEXanalyzer-package}
\alias{qPLEXanalyzer-package}
\alias{qPLEXanalyzer}
\docType{package}
\title{Tools for qPLEX-RIME data analysis}
\description{Tools for quantitiative proteomics data analysis generated from qPLEX-RIME method}
\details{
  \tabular{ll}{
    Package: \tab qPLEXanalyzer\cr
    Type: \tab Package\cr
    Version: \tab 1.0.2\cr
    Date: \tab 2018-03-07\cr
    License: \tab GPL\cr
    Depends: \tab methods\cr
  }
  
  The package offers the following functionalities
  
  Data processing, normalization & analysis:
  \itemize{
    \item convertToMSnset: Converts quantitative data to a MSnSet
    \item summarizeIntensities: Summarizes multiple peptide measurements for a protein
    \item normalizeQuantiles: Performs quantile normalization on the peptides/proteins intensities
    \item normalizeScaling: Performs scaling normalization on the peptides/proteins intensities (mean, median or sum)
    \item groupScaling: Performs scaling normalization on the peptides/proteins intensities within group (median or mean)
    \item rowScaling: Normalization by scaling peptide/protein intensity across all samples
    \item regressIntensity: Performs linear regression on protein intensities based on selected protein
    \item computeDiffStats: Compute differential statistics for the given contrasts
    \item getContrastResults: Get differential statistics results for given contrast
  }
  
  Visualization:
  \itemize{
    \item assignColours: Assigns colours to samples in groups
    \item intensityPlot: Intensity distribution plot of all the samples
    \item intensityBoxplot: Intensity distribution boxplot of all the samples
    \item peptideIntensityPlot: Peptide intensity distribution plot of specific protein
    \item pcaPlot: PCA plot of all the samples
    \item maVolPlot: MA or Volcano plot of differential analysis results
    \item corrPlot: Correlation plot of all the samples
    \item rliPlot: Relative intensity plot of all the samples
    \item hierarchicalPlot: Hierarchical clustering plot of all the samples
    \item plotMeanVar: Computes and plots mean-variance for samples in MSnSet
    \item coveragePlot: Computes and display protein sequence coverage of selected protein in proteomics experiment
  }
}
\author{
   Matthew Eldridge, Kamal Kishore (Maintainer), Ashley Sawle
  
  \email{kamal.fartiyal84@gmail.com}
}

\keyword{ package }
