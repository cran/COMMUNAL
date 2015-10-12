\name{COMMUNAL}
\alias{COMMUNAL}
\title{
Run clustering algorithms and evaluate validation metrics.
}
\description{
This functions runs various (user-specified) clustering algorithms on the data, for each potential number of clusters \code{k}. It then runs internal validation measures the quantify the fit of each clustering. The returned object is of class \code{"\linkS4class{COMMUNAL}"}, and can be used to identify 'core' clusters in the data. Currently supported clustering algorithms are those in packages "clValid", "NMF", and "ConsensusClusterPlus".

The COMMUNAL algorithm is designed to be run with \code{\link{clusterRange}}, via a call to COMMUNAL() (although this may still be useful to some researchers). After running \code{\link{clusterRange}}, use \code{\link{getGoodAlgs}} and \code{\link{getNonCorrNonMonoMeasures}} to get locally optimized clustering algorithms and validity measures. 

To determine the optimal number of clusters, use the \code{\link{plotRange3D}} function.

}
\usage{
COMMUNAL(data, ks, clus.methods = c("hierarchical", "kmeans", "diana",
                                    "som", "sota", "pam", "clara", "agnes"), 
         validation = c("Connectivity", "dunn", "wb.ratio", "g3", 
                        "g2", "pearsongamma", "avg.silwidth", "sindex"), 
         dist.metric = "euclidean", aggl.method = "ward", 
         neighb.size = 10, seed = NULL, parallel=F, gapBoot=20, 
         verbose=F, mc.cores=NULL, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data to cluster (numeric matrix or data frame). The columns are clustered, rows are features. If using cluster method \code{nmf}, all entries must be non-negative.
}
  \item{ks}{
A numeric vector of integers greater than 1, for the number of clusters to consider. For example, 2:4 tells the function to try clusterings with 2, 3, and 4 clusters.
}
  \item{clus.methods}{
Character vector of which clustering methods to use. Valid options: "\code{hierarchical}", "\code{kmeans}", "\code{diana}", "\code{fanny}", "\code{som}", "\code{model}", "\code{sota}", "\code{pam}", "\code{clara}","\code{agnes}", "\code{ccp-hc}","\code{ccp-km}", "\code{ccp-pam}", "\code{nmf}". In this list, "\code{nmf}" corresponds to "\code{nmf}" in package \pkg{NMF}, "\code{ccp-xx}" corresponds to "\code{xx}" in package 
pkg{ConsensusClusterPlus}, and the rest match to the method of the same name in package \pkg{clValid}.
}
  \item{validation}{
A character vector of the validation measures to consider. Valid options: "\code{Connectivity}", "\code{average.between}", "\code{g2}", "\code{ch}", "\code{sindex}","\code{avg.silwidth}", "\code{average.within}", "\code{dunn}", "\code{widestgap}", "\code{wb.ratio}", "\code{entropy}", "\code{dunn2}", "\code{pearsongamma}", "\code{g3}", "\code{within.cluster.ss}", "\code{min.separation}", "\code{max.diameter}", "\code{gapStatistic}". With the exception of "\code{Connectivity}", which is calculated by \code{clValid::connectivity}, and "\code{gapStatistic}", which is implemented by COMMUNAL based on cluster::clusGap(), these are calculated with \code{fpc::cluster.stats}.
}
  \item{dist.metric}{
Which metric to use when calculating the distance matrix. Used by \pkg{clValid} clustering algorithms, and in calculating validation measures. Available choices are "\code{euclidean}", "\code{correlation}", "\code{manhattan}".
}
  \item{aggl.method}{
The agglomeration method to use for "\code{hclust}" and "\code{agnes}" (if specified in \code{clus.methods}). Available choices are "\code{ward}", "\code{ward.D}", "\code{ward.D2}", "\code{single}", "\code{complete}", "\code{average}".
The ward methods have not been implemented in clValid as of this writing.
}
  \item{neighb.size}{
Numeric value. The neighborhood size used for calculating the \code{Connectivity} validation measure.
}
  \item{seed}{
Numeric value. Random seed to use in \pkg{ConsensusClusterPlus} and \pkg{NMF}.
}
  \item{parallel}{
Allows for parallel computation of the gap statistic bootstraps. WILL NOT WORK ON WINDOWS MACHINES (sorry).
}
  \item{gapBoot}{
The number of gap statistic bootstraps to perform. This recursively calls COMMUNAL for each bootstrap, though the other validation measures do not have to be calculated for each call.
}
  \item{verbose}{
Mostly output regarding clustering algorithms and gap statistic. 
}
  \item{mc.cores}{
If null, uses detectCores(). Ignored if parallel=F.
}
  \item{\dots}{
Other arguments to pass down to \pkg{ConsensusClusterPlus}, \pkg{NMF}, and \pkg{clValid}.
}
}
\value{
Return object is an object of class \code{\linkS4class{COMMUNAL}}. The class has a \code{getClustering} method to extract a data frame of cluster assignments. Alternatively, functions \code{\link{clusterKeys}} and \code{\link{returnCore}} are provided to identify core clusters. See examples below.
}
\author{
Albert Chen and Timothy E Sweeney\cr
Maintainer: Albert Chen acc2015@stanford.edu
}
\seealso{
Class \code{"\linkS4class{COMMUNAL}"}. Use functions \code{\link{clusterKeys}} and \code{\link{returnCore}} to identify core clusters.
}

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
mat.key <- clusterKeys(clusters) # get core clusters
examineCounts(mat.key)                # help decide agreement.thresh
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters (all algs agree)
table(core) # the 'core' cluster sizes
## Note: could try a different value for k to
##  see clusters with sub-optimal k

## Can specify clustering methods and validation measures
result <- COMMUNAL(data = data, ks=c(2,3),
                      clus.methods = c("diana", "som", "pam", "kmeans", "ccp-hc", "nmf"),
                      validation=c('pearsongamma', 'avg.silwidth'))
clusters <- result$getClustering(k=3)
mat.key <- clusterKeys(clusters)
examineCounts(mat.key)
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters
table(core) # the 'core' clusters

## Additional arguments are passed down to clValid, NMF, ConsensusClusterPlus
result <- COMMUNAL(data=data, ks=2:5,
                      clus.methods=c("diana", "ccp-hc", "nmf"), reps=20, nruns=2)
}
}
\keyword{ cluster }
\keyword{ robust }