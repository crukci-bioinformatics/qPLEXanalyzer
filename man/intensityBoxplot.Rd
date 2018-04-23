\name{intensityBoxplot}
\alias{intensityBoxplot}
\title{Intensity Distribution boxplot}
\description{
  Intensity distribution boxplot of all the samples}
\usage{intensityBoxplot(MSnSetObj, title="", sampleColours=NULL, colourBy="SampleGroup")}
\arguments{
   \item{MSnSetObj}{MSnSet; an object of class MSnSet}
   \item{title}{character; title of the plot}
   \item{sampleColours}{character: a named character vector of colors for samples}
   \item{colourBy}{character: column name from pData(MSnSetObj) to use for coloring samples}
}
\details{
  The column provided to the "colourBy" argument will be used to colour the samples. The colours
  will be determined using the function \link{assignColours}, alternatively the user may specify a named
  vector of colours using the "sampleColours" argument. The names of the "sampleColours" vector should
  match the unique values in the "colourBy" column.
}
\value{
An intensity distribution boxplot for quantitative proteomics data.
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
intensityBoxplot(MSnSet_data, title = "qPLEX_RIME_ER")

# custom colours
customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
intensityBoxplot(MSnSet_data, title = "qPLEX_RIME_ER", sampleColours = customCols)
}
