\name{getContrastResults}
\alias{getContrastResults}
\title{Get differential statistics results}
\description{
  Get differential statistics results for given contrasts.
}
\usage{getContrastResults(diffstats, contrast, controlGroup = NULL, ann, applyLog2Transform = TRUE, writeFile= FALSE)}
\arguments{
  \item{diffstats}{list; Output of computeDiffStats function}
  \item{contrast}{character; Contrast of interest to see differential statistics results}
  \item{controlGroup}{character; control group such as IgG}
  \item{ann}{dataframe; A dataframe of protein annotation of four columns: 
                        "Protein", "Gene", "Description" and "GeneSymbol"}
  \item{applyLog2Transform}{Logical; apply log2 transformation to the raw intensitites}                      
  \item{writeFile}{logical; whether to write the results into a text file}
}
\value{A dataframe object and text file containing the result of the differential statistics.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
MSnset_norm <- groupScaling(MSnset_data,func=median)
MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
diffstats <- computeDiffStats(data=MSnset_Pnorm, contrasts=contrasts)
diffexp <- getContrastResults(diffstats=diffstats, contrast=contrasts, controlGroup = NULL, ann= human_anno, writeFile= FALSE)
}
