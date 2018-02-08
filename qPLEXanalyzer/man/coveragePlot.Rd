\name{coveragePlot}
\alias{coveragePlot}
\title{Plot protein sequence coverage}
\description{
  Computes and display protein sequence coverage in proteomics experiment
}
\usage{coveragePlot(data, ProteinID, name,fastaFile, col="brown")}
\arguments{
  \item{data}{MSnSet: An object of class MSnSet}
  \item{ProteinID}{character: Uniprot ID of the protein}
  \item{name}{character: Name of the protein}
  \item{fastaFile}{character: fasta file of protein sequence}
  \item{col}{character: Color for plotting}
}
\details{
  In the qPLEX-RIME experiment it is imperative for bait protein to have good sequence coverage. 
  This function plots the protein sequence coverage of the bait protein in the qPLEX-RIME
  experiment. It requires the fasta sequence file of bait protein as input to generate
  the plot.
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
mySequenceFile <- system.file('extdata', "P03372.fasta", package="qPLEXanalyzer")
coveragePlot(MSnset_data, ProteinID="P03372", name="ERa", fastaFile=mySequenceFile)
}
