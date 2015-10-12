\name{examineCounts}
\alias{examineCounts}
\title{
Examine algorithm cluster agreement
}
\description{
Returns a data frame to help the user choose the \code{agreement.thresh} parameter for \code{\link{returnCore}}. Each row shows how many samples (\code{sample.count}) are at each agreement level (\code{percent.agreement}), and what percent of the data would be removed for that agreement threshold (\code{percent.remaining.if.removed}).
}
\usage{
examineCounts(mat.key)
}
\arguments{
  \item{mat.key}{
A matrix of cluster assignments, where rows are items and columns are algorithms. The output of \code{\link{clusterKeys}}.
}
}
\value{
A 3 column table, showing the number of samples at each agreement level and the percent of data above that agreement level.
The agreement level for a sample is the (highest) fraction of algorithms which agree on its cluster assignment.
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
core <- returnCore(mat.key, agreement.thresh=66) # need at least 2 of 3 algs to agree
table(core)
core <- returnCore(mat.key, agreement.thresh=99) # need all algs to agree
table(core)
}