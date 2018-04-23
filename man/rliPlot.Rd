\name{rliPlot}
\alias{rliPlot}
\title{RLI plot}
\description{
  Relative log intensity (RLI) plots of the samples within MSnset
}
\usage{rliPlot(MSnSetObj, title="", sampleColours=NULL, colourBy="SampleGroup", omitIgG=TRUE)}
\arguments{
  \item{MSnSetObj}{MSnSet; an object of class MSnSet}
  \item{title}{character: title for the plot}
  \item{sampleColours}{character: a named vector of colours for samples}
  \item{colourBy}{character: column name to use for colouring samples from pData(MSnSetObj)}
  \item{omitIgG}{logical: whether to remove IgG from the RLI plot}
}
\value{An instance of class "MSnSet".}
\details{
  An RLI-plot is a boxplot that can be used to visualise unwanted variation in a data set. It is similar 
  to the relative log expression plot developed for microarray analysis - see Gandolfo and Speed (2018).
  Rather than examining gene expression, the RLI plot uses the MS intensities for each peptide or the 
  summarised protein intensities.
  
  The column provided to the "colourBy" argument will be used to colour the samples. The colours
  will be determined using the function \link{assignColours}, alternatively the user may specify a named
  vector of colours using the "sampleColours" argument. The names of the "sampleColours" vector should
  match the unique values in the colourBy column.
}
\references{Gandolfo LC, Speed TP (2018) RLE plots: Visualizing unwanted variation in high dimensional data. PLoS ONE 13(2): e0191629. \url{https://doi.org/10.1371/journal.pone.0191629}}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
rliPlot(MSnSet_data, title = "qPLEX_RIME_ER")

# custom colours
customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
rliPlot(MSnSet_data, title = "qPLEX_RIME_ER", sampleColours = customCols)
}
