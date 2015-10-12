\name{COMMUNAL-package}
\alias{COMMUNAL-package}
\docType{package}
\title{
COmbined Mapping of Multiple clUsteriNg ALgorithms
}
\description{
This package allows for identification of optimal clustering for a data set. It provides a framework to run a wide range of clustering algorithms to determine the optimal number (k) of clusters in the data. It then provides a function to analyze the cluster assignments from each clustering algorithm to identify samples that repeatedly classify to the same group. We call these 'core clusters,' leading to optimal beds for later class discovery.
}
\details{
\tabular{ll}{
Package: \tab COMMUNAL\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-08-12\cr
License: \tab GPL-2\cr
Imports: \tab clValid, fpc, methods\cr
Depends: \tab R (>= 2.10), cluster\cr
Suggests: \tab RUnit, NMF, ConsensusClusterPlus, rgl\cr
}
Start with a matrix of data to cluster. Important functions are:\cr
\code{\link{COMMUNAL}} to run clustering algorithms once\cr
\code{\link{clusterRange}} to run COMMUNAL across increasing subsets of data \cr
\code{\link{getGoodAlgs}} to identify robust algorithms from \code{clusterRange} \cr
\code{\link{getNonCorrNonMonoMeasures}} to identify non-monotonic, non-correlated validity measures from \code{clusterRange} \cr
\code{\link{plotRange3D}} to pick k\cr
\code{\link{clusterKeys}} to identify core clusters\cr
\code{\link{returnCore}} to identify core clusters\cr

}
\author{
Albert Chen, Timothy E Sweeney, Olivier Gevaert\cr
Maintainer: Albert Chen acc2015@stanford.edu
}
\keyword{ package }
\examples{
\dontrun{
## create artificial data set with 3 distinct clusters
set.seed(1)
V1 = c(abs(rnorm(100, 2)), abs(rnorm(100, 50)), abs(rnorm(100, 140)))
V2 = c(abs(rnorm(100, 2, 8)), abs(rnorm(100, 55, 4)), abs(rnorm(100, 105, 1)))
data <- t(data.frame(V1, V2))
colnames(data) <- paste("Sample", 1:ncol(data), sep="")
rownames(data) <- paste("Gene", 1:nrow(data), sep="")

## run COMMUNAL
result <- COMMUNAL(data=data, ks=seq(2,5))  # result is a COMMUNAL object
k <- 3                                # suppose optimal cluster number is 3
clusters <- result$getClustering(k)   # method to extract clusters
mat.key <- clusterKeys(clusters)      # get core clusters
examineCounts(mat.key)                # help decide agreement.thresh
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters (all algs agree)
table(core) # the 'core' clusters

## Additional arguments are passed down to clValid, NMF, ConsensusClusterPlus
result <- COMMUNAL(data=data, ks=2:5,
                      clus.methods=c("diana", "ccp-hc", "nmf"), reps=20, nruns=2)

## To identify k, use clusterRange and plotRange3D to visualize validation measures
data(BRCA.100) # 533 tissues to cluster, with measurements of 100 genes each
varRange <- c(10,25,50,75,100)
meas <- c("Connectivity", "average.between",
          "ch", "sindex", "avg.silwidth", 
          "average.within", "dunn", "widestgap",
          "wb.ratio", "entropy", "dunn2", 
          "pearsongamma", "g3", "within.cluster.ss", 
          "min.separation", "max.diameter")
BRCA.results <- clusterRange(BRCA.100, ks=2:6, varRange=varRange, validation=meas)

goodMeasures <- getNonCorrNonMonoMeasures(BRCA.results)
goodAlgs <- getGoodAlgs(BRCA.results)

plot.data <- plotRange3D(BRCA.results, goodAlgs=goodAlgs, goodMeasures = goodMeasures)
}
}
