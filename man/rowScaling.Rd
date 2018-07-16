\name{rowScaling}
\alias{rowScaling}
\title{Normalization by scaling peptide/protein intensity across all samples}
\description{
  Divide each peptide/protein by the row mean/median and transform to log2
}
\usage{rowScaling(MSnSetObj, scalingFunction)}
\arguments{
  \item{MSnSetObj}{MSnSet; an object of class MSnSet}
  \item{scalingFunction}{function; median or mean}
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
MSnset_norm <- rowScaling(MSnSet_data, scalingFunction=median)
}
