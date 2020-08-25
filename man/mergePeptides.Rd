\name{mergePeptides}
\alias{mergePeptides}
\title{Merge identical peptides intensities}
\description{
  Merge identical peptides intensities belonging to same protein to single peptide intensity. This function is especially useful for phosphopeptide analysis.
}
\usage{mergePeptides(MSnSetObj, summarizationFunction, annotation, PosMasterProt=NULL)}
\arguments{
  \item{MSnSetObj}{MSnSet; an object of class MSnSet}
  \item{summarizationFunction}{function; method used to aggregate the peptides. Sum, mean or median}
  \item{annotation}{data.frame; a data.frame of protein annotation of four columns: 
                        "Accessions", "Gene", "Description" and "GeneSymbol"}
  \item{PosMasterProt}{numeric; Denote the column number of "Position in Master Protein" within fData(MSnSetObj)}
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
MSnset_P <- mergePeptides(MSnSet_data, sum, human_anno)
}
