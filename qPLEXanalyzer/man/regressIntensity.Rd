\name{regressIntensity}
\alias{regressIntensity}
\title{Regression based analysis}
\description{
  Performs linear regression on protein intensities based on selected protein (qPLEX-RIME bait)
}
\usage{regressIntensity(data,controlInd=NULL,ProteinId)}
\arguments{
  \item{data}{MSnSet; An object of class MSnSet}
  \item{controlInd}{numeric; Index of IgG within MSnSet}
  \item{ProteinId}{character; Uniprot protein ID}
}
\details{
This function performs regression based analysis upon protein intensities based on selected protein. In qPLEX RIME this method could be used to regress out affect of target protein on other interactors. This function corrects this dependency of many proteins on the target protein levels by linear regression. It sets target protein levels as the independent variable (x) and all other proteins as dependent variable (y). The resulting residuals of the linear regression y=ax+b are the protein levels corrected for target protein dependency.
}
\value{An instance of class "MSnSet". This consists of corrected protein levels. In addtion the function also plots histogram of correlation of target protein with all other proteins before and after this correction.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
MSnset_P <- summarizeIntensities(MSnset_data, sum, human_anno)
IgG_ind <- which(pData(MSnset_P)$Label == "IgG")
MSnset_reg <- regressIntensity(MSnset_P,controlInd=IgG_ind,ProteinId="P03372")
}
