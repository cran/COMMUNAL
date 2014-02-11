\name{plotRange3D}
\alias{plotRange3D}
\title{
Visualize cluster quality
}
\description{
Given \code{\link{COMMUNAL}} outputs for a dataset, summarize the suitability of a
range of cluster numbers according to validation measures. Input can be the output of
\code{clusterRange}.
}
\usage{
plotRange3D(test_range, ks, goodAlgs, goodMeasures, filename=NULL, ...)
}
\arguments{
  \item{test_range}{
list of 2 items, where first element is list of \code{"\linkS4class{COMMUNAL}"} objects,
and second is a vector of the number of data points clustered for each element in the list. Can be the output of \code{\link{clusterRange}}.
}
  \item{ks}{
range of cluster number k to plot. Must be subset of ks supplied to \code{\link{COMMUNAL}}.
}
  \item{goodAlgs}{
which algorithms to use in summarizing validation measures. Must be subset of 
algorithms supplied to \code{\link{COMMUNAL}}.
}
  \item{goodMeasures}{
which validation measures to use in summarizing validation measures. Must be subset of 
measures supplied to \code{\link{COMMUNAL}}.
}
  \item{filename}{
optionally specify filename to save a snapshot of the 3D image.
}
  \item{\dots}{
other arguments to pass down to 2D plot of mean z-score against k.
}
}
\details{
A summarized validation measure value is computed for each value of k, for each dataset.
This is done by first subsetting the data to the measures, ks, and algorithms of interest,
and then computing averages of the measures for each dataset and k (number of clusters).

For some validation measures, a lower value implies better clustering, and for others a 
higher value is better. Prior to averaging, measures that favor a lower value are 
multiplied by negative one. Furthermore, each measure is scaled to have zero mean and
unit variance across all the datasets prior to averaging, so each measure has equal
weight, and we can compare the plot across datasets.

If the wrong number of clusters was identified
}
\value{
A 3D plot is generated, using the package "rgl".
A matrix of the plotted values is returned.
A 2D plot of average metric against k is also generated.
}

\author{
Albert Chen and Timothy E Sweeney\cr
Maintainer: Albert Chen acc2015@stanford.edu
}

\examples{
\dontrun{
data(BRCA.results) # clusterRange output for breast cancer dataset
goodAlgs <- c("hierarchical", "kmeans", "model", "agnes", "som")
goodMeasures <- c('wb.ratio', 'avg.silwidth', 'dunn')
ks <- 2:8
(values <- plotRange3D(BRCA.results, ks, goodAlgs, goodMeasures))
}
}
