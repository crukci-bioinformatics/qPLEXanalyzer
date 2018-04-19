\name{normalizeQuantiles}
\alias{normalizeQuantiles}
\title{Quantile normalization}
\description{
  Performs quantile normalization on the intensities within columns
}
\usage{normalizeQuantiles(MSnSetObj)}
\arguments{
  \item{MSnSetObj}{MSnSet; an object of class MSnSet}
}
\details{
  The peptide intensities are roughly replaced by the order statics on their abundance. 
  This normalization technique has the effect of making the distributions of intensities
  from the different samples identical in terms of their statistical properties. It is 
  the strongest normalization method and should be used carefully as it erases most of the
  difference between the samples.
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),Sequences=2,Accessions=6)
MSnset_norm <- normalizeQuantiles(MSnSet_data)
}

