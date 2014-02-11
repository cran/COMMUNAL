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
Version: \tab 1.0\cr
Date: \tab 2015-01-05\cr
License: \tab GPL-2\cr
Imports: \tab clValid, fpc, methods\cr
Depends: \tab R (>= 2.10), cluster\cr
Suggests: \tab RUnit, NMF, ConsensusClusterPlus, rgl\cr
}
Start with a matrix of data to cluster. Important functions are:\cr
\code{\link{COMMUNAL}} to run clustering algorithms\cr
\code{\link{clusterRange}} to run clustering algorithms (harness for COMMUNAL) \cr
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
mat.key <- clusterKeys(clusters, k=k) # get core clusters
examineCounts(mat.key)                # help decide agreement.thresh
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters (all algs agree)
table(core) # the 'core' clusters

## Additional arguments are passed down to clValid, NMF, ConsensusClusterPlus
result <- COMMUNAL(data=data, ks=2:5,
                      clus.methods=c("diana", "ccp-hc", "nmf"), reps=20, nruns=2)

## To identify k, use clusterRange and plotRange3D to visualize validation measures
data(BRCA.100) # 533 tissues to cluster, with measurements of 100 genes each
varRange <- c(50, 75, 100)
clus.methods <- c("hierarchical", "kmeans")
validation <- c('wb.ratio', 'dunn', 'avg.silwidth')
range.results <- clusterRange(BRCA.100, varRange, ks=2:5, clus.methods=clus.methods,
                              validation=validation)
plot.data <- plotRange3D(range.results, ks=2:5, clus.methods, validation)
}
}
