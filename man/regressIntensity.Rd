\name{regressIntensity}
\alias{regressIntensity}
\title{Regression based analysis}
\description{
  Performs linear regression on protein intensities based on selected protein (qPLEX-RIME bait)
}
\usage{regressIntensity(MSnSetObj, ProteinId, controlInd=NULL, plot=TRUE)}
\arguments{
  \item{MSnSetObj}{MSnSet; an object of class MSnSet}
  \item{controlInd}{numeric; index of IgG within MSnSet}
  \item{ProteinId}{character; Uniprot protein ID}
  \item{plot}{character; Whether or not to plot the QC histograms}
}
\details{
This function performs regression based analysis upon protein intensities based on a selected protein. In qPLEX RIME this method could be used to regress out the effect of target protein on other interactors. This function corrects this dependency of many proteins on the target protein levels by linear regression. It sets the target protein levels as the independent variable (x) and each of the other proteins as the dependent variable (y). The resulting residuals of the linear regressions y=ax+b are the protein levels corrected for target protein dependency.
}
\value{An instance of class "MSnSet". This consists of corrected protein levels. In addition the function can also plot histograms of correlation of target protein with all other proteins before and after this correction.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
MSnset_P <- summarizeIntensities(MSnSet_data, sum, human_anno)
IgG_ind <- which(pData(MSnset_P)$SampleGroup == "IgG")
MSnset_reg <- regressIntensity(MSnset_P, controlInd=IgG_ind, ProteinId="P03372")
}
