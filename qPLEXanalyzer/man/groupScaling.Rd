\name{groupScaling}
\alias{groupScaling}
\title{Normalization by scaling within group}
\description{
  Performs scaling normalization on the intensities within group (median or mean)
}
\usage{groupScaling(data,func, Grp="Label")}
\arguments{
  \item{data}{MSnSet; An object of class MSnSet}
  \item{func}{function; median or mean}
  \item{Grp}{character; the feature on which groups would be based; example "Label"}
}
\details{
  In this normalization method the central tendencies (mean or median) of the samples within groups 
  are aligned. The argument "Grp" is used to define separate groups to normalize. The function takes
  one of the column of pData(data) as the variable for classifying group. The default variable is "Label". 
  It is imperative in qPLEX-RIME experiment to define IgG as a separate group and normalize it separately 
  from others. You could add a column into the metadata to define this classification.
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
MSnset_norm <- groupScaling(MSnset_data,func=median, Grp="Label")
}
