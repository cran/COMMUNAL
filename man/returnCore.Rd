\name{returnCore}
\alias{returnCore}
\title{
Assign points to core clusters
}
\description{
Takes keyed cluster assignments of each algorithm and makes final cluster assignmenmts
by majority vote. If the majority percentage is <= the specified threshold, the cluster cannot be 
determined and is set to 0.
}
\usage{
returnCore(mat.key, agreement.thresh = 50)
}
\arguments{
  \item{mat.key}{
Matrix of rekeyed cluster assignments. E.g., the output of \code{\link{clusterKeys}}.
}
  \item{agreement.thresh}{
Percent of algorithms required to agree for a cluster assignment to be accepted.
Otherwise, cluster is set to 0. By default at least half the algorithms must agree.
}
}
\details{
Can use \code{\link{examineCounts}} to help determine the \code{agreement.thresh} argument.
}
\author{
Albert Chen and Timothy E Sweeney\cr
Maintainer: Albert Chen acc2015@stanford.edu
}
\examples{
# reindexes cluster numbers to agree
clusters <- data.frame(
  alg1=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)),
  alg2=as.integer(c(1,1,1,1,1,3,3,3,3,3,2,2,2,2,2)),
  alg3=as.integer(c(3,3,3,3,3,1,1,1,1,1,2,2,2,2,2))
)
mat.key <- clusterKeys(clusters)
mat.key # cluster indices are relabeled
examineCounts(mat.key)
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters
table(core) # the 'core' clusters

# some clusters assignments are undetermined
clusters <- data.frame(
  alg1=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,1,1,2,2,3,3)),
  alg2=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,1,2,2,3,3,1)),
  alg3=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,2,3,1,1,2,3))
)
mat.key <- clusterKeys(clusters)
mat.key # last six samples have conflicting assignments
examineCounts(mat.key)
(core <- returnCore(mat.key, agreement.thresh=66)) # need at least 2 of 3 algs to agree
table(core)
(core <- returnCore(mat.key, agreement.thresh=99)) # need all algs to agree
table(core)
}