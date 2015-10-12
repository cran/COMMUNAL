\name{Choose Optimal Clustering Algorithms}
\alias{testAlgsMinSize}
\alias{getGoodAlgs}
\title{
Choose Optimal Clustering Algorithms
}
\description{
These functions serve to find locally optimal clustering algorithms for an output of \code{clusterRange}. They measure the percentage of clusterings for which a given algorithm returned any clusters with minSize or fewer members. 
}
\usage{
testAlgsMinSize(clusRange, algs = "all", minSize = 3)

getGoodAlgs(clusRange, algs = "all", minSize = 3)
}
\arguments{
  \item{clusRange}{
The output from a call to \code{\link{clusterRange}}.
}
  \item{algs}{
The algorithms to examine. Defaults to "all" algorithms present in clusRange, but user can define a subset (character vector). 
}
  \item{minSize}{
The size at or below which a cluster is considered sub-optimal.
}
}
\details{
\code{testAlgsMinSize} iterates over the range of sub-datasets present in clusRange, and will print cluster assignment counts for all K where a cluster smaller than minSize has been returned (to help the user identify patterns in the data). This makes it quite verbose. It will then return the mean percentage of clusters < minSize for all algorithms.

\code{getGoodAlgs} is a wrapper for \code{testAlgsMinSize} that returns the algorithms that are at or under the mean for percent clusterings < minSize. 
}
\value{
\code{testAlgsMinSize} returns a named numeric vector of mean percentage of clusters < minSize.

\code{getGoodAlgs} returns a character vector of algorithms that are at or under the mean for percent clusterings < minSize, suitable for passing into \code{\link{plotRange3D}}
}
\author{
Sweeney
}
\seealso{
\code{\link{getNonCorrNonMonoMeasures}}
}
\examples{
## output from running \code{clusterRange} on data(BRCA.100)
data(BRCA.results)

## BRCA results does not have any clusters < minSize=3, so returns all.
testAlgsMinSize(BRCA.results)
getGoodAlgs(BRCA.results)

## force output, call minSize=50 (just to test here)
getGoodAlgs(BRCA.results, minSize=50)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
