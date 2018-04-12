\name{maVolPlot}
\alias{maVolPlot}
\title{MA or Volcano Plot}
\description{
  MA or Volcano plot of differential statistics results
}
\usage{maVolPlot(diffstats, contrast, title="", controlGroup = NULL,
                      selectedGenes=NULL, fdrCutOff=0.05,
                      lfcCutOff=1, controlLfcCutOff=1, plotType="MA")}
\arguments{
   \item{diffstats}{list; Output of computeDiffStats function}
   \item{contrast}{Character; Contrast of interest to plot differential statistics results}
   \item{title}{Character: Title for the plot}
   \item{controlGroup}{Character; control group such as IgG}
   \item{selectedGenes}{Character: a vector defining genes to plot}
   \item{fdrCutOff}{numeric: cutoff for adj.P.Val}
   \item{lfcCutOff}{numeric: cutoff for log2FC}
   \item{controlLfcCutOff}{numeric: plot genes above controlLogFoldChange cutoff}
   \item{plotType}{Character: "MA" or "Volcano"}
}
\value{
A "MA" or "Volcano" plot to visualize differential statistics results.
}
\examples{
data(human_anno)
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),Sequences=2,Accessions=6)
MSnset_norm <- groupScaling(MSnSet_data,func=median)
MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
diffstats <- computeDiffStats(MSnset_Pnorm, contrasts=contrasts)
maVolPlot(diffstats, contrast = contrasts, plotType="MA")
maVolPlot(diffstats, contrast = contrasts, plotType="Volcano")
}
