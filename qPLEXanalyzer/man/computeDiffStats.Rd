\name{computeDiffStats}
\alias{computeDiffStats}
\title{Compute differential statistics}
\description{
  Compute differential statistics on the given contrasts, based on \code{\link{limma}} functions.
}
\usage{computeDiffStats(data, batchEffect = NULL, applyLog2Transform = TRUE, contrasts, trend = TRUE, robust = TRUE)}
\arguments{
  \item{data}{MSnSet; An object of class MSnSet}
  \item{batchEffect}{character; vector of variable(s) to correct for batch effect, Default : "Label"}
   \item{applyLog2Transform}{Logical; apply log2 transformation to the raw intensitites} 
  \item{contrasts}{character; comparisions for differential statistics}
  \item{trend}{logical; TRUE or FALSE}
  \item{robust}{logical; TRUE or FALSE}
}
\details{
  A statistical analysis for the identification of differentially regulated or bound proteins is 
  carried out using limma based analysis. It uses linear models to assess differential expression 
  in the context of multifactor designed experiments.  Firstly, a linear model is fitted for each 
  protein where the model includes variables for each group and MS run. Then, log2 fold changes 
  between comparisions are estimated. Multiple testing correction of p-values are applied using 
  the Benjamini-Hochberg method to control the false discovery rate (FDR).

  In order to correct for batch effect, variable(s) can be defined. It should corresponds to a 
  column name in pData(data). The default variable is "Label" that distinguish between two groups. 
  If more variables are defined they are added to default.
}
\value{A list object containing three components: data (\R{MSnSet} object), fittedLM (fitted linear model) 
and fittedContrasts. This object should be input into getContrastResults function to get differential results. 
See \code{\link{eBayes}} function of \code{\link{limma}} for more details on diffential statistics.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
MSnset_norm <- groupScaling(MSnset_data,func=median)
MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
diffstats <- computeDiffStats(data=MSnset_Pnorm, contrasts=contrasts)
}
