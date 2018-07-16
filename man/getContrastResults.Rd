\name{getContrastResults}
\alias{getContrastResults}
\title{Get differential statistics results}
\description{
  Get differential statistics results for given contrasts.
}
\usage{getContrastResults(diffstats, contrast, controlGroup = NULL, transform = TRUE, writeFile= FALSE)}
\arguments{
  \item{diffstats}{list; output of computeDiffStats function}
  \item{contrast}{character; contrast of interest for which to retrieve differential statistics results}
  \item{controlGroup}{character; control group such as IgG}
  \item{transform}{logical; apply log2 transformation to the raw intensities}                      
  \item{writeFile}{logical; whether to write the results into a text file}
}
\value{A \code{\link{data.frame}} object and text file containing the result of the differential statistics.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
MSnset_norm <- groupScaling(MSnSet_data, scalingFunction=median)
MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
diffstats <- computeDiffStats(MSnset_Pnorm, contrasts=contrasts)
diffexp <- getContrastResults(diffstats=diffstats, contrast=contrasts)
}
