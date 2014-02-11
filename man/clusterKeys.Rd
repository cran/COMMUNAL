\name{clusterKeys}
\alias{clusterKeys}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rekey cluster assignments.
}
\description{
Reindexes (rekeys) the cluster assignments to maximize overlap across algorithms.\cr
Ignores algorithms which could not find \code{k} clusters; i.e. when one of the clusters is smaller than the \code{min.size} argument.\cr Use this after determining the number of clusters.
}
\usage{
clusterKeys(clusters, k, min.size = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters}{
Data frame of cluster assignments, where rows are samples, columns are algorithms, assignments are integers. For example, the output of the \code{getClustering} method in \code{"\linkS4class{COMMUNAL}"}.
}
  \item{k}{
Number of clusters selected.
}
  \item{min.size}{
  Minimum cluster size. Algorithms that return clusters smaller than this (or that don't have \code{k} clusters) are tossed out.
}
}
\value{
Returns a matrix of rekeyed cluster assignments, such that cluster 'n' refers to the same cluster across all algorithms. Cluster 0 contains the samples for which no consistent 'core' cluster could be identified.
}
\author{
Albert Chen and Timothy E Sweeney\cr
Maintainer: Albert Chen acc2015@stanford.edu
}
\examples{
# reindexes cluster numbers to agree
k <- 3
clusters <- data.frame(
  alg1=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)),
  alg2=as.integer(c(1,1,1,1,1,3,3,3,3,3,2,2,2,2,2)),
  alg3=as.integer(c(3,3,3,3,3,1,1,1,1,1,2,2,2,2,2))
)
mat.key <- clusterKeys(clusters, k)
mat.key # cluster indices are relabeled
examineCounts(mat.key)
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters
table(core) # the 'core' clusters

# some clusters assignments are undetermined
k <- 3
clusters <- data.frame(
  alg1=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,1,1,2,2,3,3)),
  alg2=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,1,2,2,3,3,1)),
  alg3=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,2,3,1,1,2,3))
)
mat.key <- clusterKeys(clusters, k)
mat.key # last six samples have conflicting assignments
examineCounts(mat.key)
core <- returnCore(mat.key, agreement.thresh=66) # at least 2 of 3 algs must agree
table(core)
core <- returnCore(mat.key, agreement.thresh=99) # all algs must agree
table(core)
}