\name{pcaPlot}
\alias{pcaPlot}
\title{PCA plot}
\description{
  PCA plots of the samples within MSnset
}
\usage{pcaPlot(MSnSetObj, omitIgG=FALSE, sampleColours=NULL, transFunc=log2xplus1, transform=TRUE, 
                    colourBy="SampleGroup", title="", labelColumn="BioRep", labelsize=4, pointsize=4, x.nudge=4, x.PC=1)}
\arguments{
  \item{MSnSetObj}{MSnSet; An object of class MSnSet}
  \item{omitIgG}{Logical: Whether to remove IgG from the PCA plot}
  \item{sampleColours}{Character: A named vector of colours for samples}
  \item{transFunc}{func: internal helper function for log transformation}
  \item{transform}{Logical: Whether to log transform intensities}
  \item{colourBy}{Character: column name to use for colouring samples from pData(MSnSetObj)}
  \item{title}{Character: title for the plot}
  \item{labelColumn}{Character: column name from pData(MSnSetObj) to use for labelling points on the plot}
  \item{labelsize}{Numeric: Size of the labels}
  \item{pointsize}{Numeric: Size of plotting points}
  \item{x.nudge}{Numeric: distance to move labels along the x-axis away from the plotting points}
  \item{x.PC}{Numeric: The principle component to plot on the x-axis; the following PC will be plotted on the y-axis}
}
\value{An instance of class "MSnSet".}
\details{
  The column provided to the "colourBy" argument will be used to colour the samples. The colours
  will be determined using the function \link{assignColours}, alternatively the user may specify a named
  vector of colours using the "sampleColours" argument. The names of the "sampleColours" vector should
  match the unique values in the colourBy column.
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
exprs(MSnSet_data) <- exprs(MSnSet_data)+0.01
pcaPlot(MSnSet_data, omitIgG = TRUE, labelColumn = "BioRep")
}
