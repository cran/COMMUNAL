\name{plotRange3D}
\alias{plotRange3D}
\title{
Visualize cluster stability
}
\description{
Given \code{\link{clusterRange}} output for a dataset, visualize the cluster optimality output for a range of K over a range of variable measures.
}
\usage{
plotRange3D(clusRange, ks=NULL, goodAlgs=NULL, goodMeasures=NULL, filename=NULL,
            colorbar=T, minSize = 3, plot3D=T, ...)
}
\arguments{
  \item{clusRange}{
The output of \code{\link{clusterRange}}.
}
  \item{ks}{
range of cluster number k to plot. If NULL, plots all ks in clusRange.
}
  \item{goodAlgs}{
which algorithms to use in summarizing validation measures. If NULL, plots all algorithms in clusRange.
}
  \item{goodMeasures}{
which validation measures to use in summarizing validation measures. If NULL, plots all validity measures in clusRange.
}
  \item{filename}{
optionally specify filename to save a snapshot of the 3D image.
}
  \item{colorbar}{
Whether to draw the front right color legend in the output.
}
  \item{minSize}{
The minimum acceptable size for a cluster. \code{plotRange3D} tests all algorithms for clusters <minSize and prints these to command line. The purpose of this is to alert the user that the passed algorithms may not be ideal. We recommend calling \code{\link{getGoodAlgs}} to decrease the number of non-robust algorithms. Having some clusters <minSize may be unavoidable, but should be minimized.
}
  \item{plot3D}{
Whether to do a 3D plot at all. If rgl() will not work locally, this allows the user to still see the 2D plot.
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

}
\value{
A 3D plot is generated, using the package "rgl".
A matrix of the plotted values is returned.
A 2D plot of average metric against k (i.e., the mean over varRange) is also generated.
}
\note{
rgl is a finnicky (but great!) package; we do not maintain it. Please address rgl questions to the rgl package maintainers.
}
\author{
Albert Chen and Timothy E Sweeney\cr
Maintainer: Albert Chen acc2015@stanford.edu
}

\examples{
\dontrun{
## clusterRange output for breast cancer dataset
data(BRCA.results) 

## automated selection of optimal algorithms and validity measures
goodAlgs <- getGoodAlgs(BRCA.results)
goodMeasures <- getNonCorrNonMonoMeasures(BRCA.results)

(values <- plotRange3D(BRCA.results, goodAlgs, goodMeasures))
}
}
