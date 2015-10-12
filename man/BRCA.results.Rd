\name{BRCA.results}
\alias{BRCA.results}
\docType{data}
\title{
TCGA Breast Cancer subset \code{\link{clusterRange}} output
}
\description{
Output of \code{\link{clusterRange}} on breast cancer gene expression test dataset (100 genes, 533 samples) (\code{\link{BRCA.100}}).
To be displayed with \code{\link{plotRange3D}} for visualization.
}
\usage{data("BRCA.results")}
\format{
  List of 2. First element is a list of 3 \code{\linkS4class{COMMUNAL}} objects. Second element is a vector of the corresponding number of data points that were clustered.
}
\examples{
## BRCA.results <- clusterRange(BRCA.100, ks=2:6, varRange=c(10,25,50,75,100))
data(BRCA.results)
}
\keyword{datasets}
