\name{Measure Cluster Monotonicty}
\alias{getNonCorrNonMonoMeasures}
\alias{monotoneClusterRange}
\alias{measuresCorr}
\title{
Measure cluster monotonicty within the output from clusterRange. 
}
\description{
These three functions serve to automatically select non-correlated, non-monotonic measures from the output of a run of \code{clusterRange}. The user is only required to run \code{getNonCorrNonMonoMeasures}; natively, this will return 'numMeasures' measures, selected from all measures that are at or below the mean percentage monotonicity. 

The functions \code{measuresCorr} and \code{monotoneClusterRange} allow for other selection methods by advanced users.
}
\usage{
monotoneClusterRange(clusRange, goodAlgs = NULL)

measuresCorr(clusRange, goodMeasures)

getNonCorrNonMonoMeasures(clusRange, goodMeasures = "all", goodAlgs = NULL,
                          numMeasures = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusRange}{
The output from a call to \code{\link{clusterRange}}.
}
  \item{goodMeasures}{
A subset of validity measures pre-selected by the user. Defaults to "all". 
}
  \item{goodAlgs}{
A subset of clustering algorithms pre-selected by the user. The defauly (NULL) will use all algorithms present in the clusRange object. 
}
  \item{numMeasures}{
The number of clusters to be returned from the clustering of correlations. Literally, call 'cutree' at k=numMeasures on the hierarchical clustering of the non-monotonic measures. 
}
}
\details{
Advanced users may wish to select validation measures by a guideline of their choosing; these functions can assist in this process as well. In particular, the COMMUNAL manuscript shows a plot of correlations of measures; this was made using the corrplot::corrplot() function on the output from measuresCorr. Such visualization can be helpful in selecting numMeasures.
}
\value{
\code{monotoneClusterRange} outputs a named numeric vector of the percent monotonicity of all validity measures in clusRange.

\code{measuresCorr} outputs a correlation matrix of all goodMeasures in clusRange.

\code{getNonCorrNonMonoMeasures} outputs a character vector of validation measures, suitable for passing into \code{\link{plotRange3D}}
}
\references{
See COMMUNAL manuscript, Sweeney et al., 2015.
}
\author{
Sweeney
}
\seealso{
\code{\link{getGoodAlgs}}
}
\examples{
## output from running \code{clusterRange} on data(BRCA.100)
data(BRCA.results)

monotoneClusterRange(BRCA.results)

measuresCorr(BRCA.results, goodMeasures=c("avg.silwidth", "g3", "sindex"))

getNonCorrNonMonoMeasures(BRCA.results)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
