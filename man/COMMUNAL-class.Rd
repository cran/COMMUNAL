\name{COMMUNAL-class}
\Rdversion{1.1}
\docType{class}
\alias{COMMUNAL-class}

\title{Class \code{"COMMUNAL"}}
\description{
This is the class of the object returned by the \code{\link{COMMUNAL}} function. It contains all the clustering results, including the values of the validation measures, and the cluster assignments. Cluster assignments for a particular number of clusters can be conveniently extracted with the \code{getClustering(k)} method.
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\author{
Albert Chen and Timothy E Sweeney\cr
Maintainer: Albert Chen acc2015@stanford.edu
}

\seealso{
Vignette and \code{\link{COMMUNAL}}. The function \code{\link{COMMUNAL}} returns an object of this class.
}

\examples{
showClass("COMMUNAL")
\dontrun{
## create artificial data set with 3 distinct clusters
set.seed(1)
V1 = c(abs(rnorm(100, 2)), abs(rnorm(100, 50)), abs(rnorm(100, 140)))
V2 = c(abs(rnorm(100, 2, 8)), abs(rnorm(100, 55, 4)), abs(rnorm(100, 105, 1)))
data <- t(data.frame(V1, V2))
colnames(data) <- paste("Sample", 1:ncol(data), sep="")
rownames(data) <- paste("Gene", 1:nrow(data), sep="")

## run COMMUNAL with defaults
result <- COMMUNAL(data=data, ks=seq(2,5))  # result is a COMMUNAL object
k <- 3                                # suppose optimal cluster number is 3
clusters <- result$getClustering(k)   # extract clusters
mat.key <- clusterKeys(clusters)      # get core clusters
examineCounts(mat.key)                # all algorithms agree
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' cluster assignments
table(core) # the 'core' cluster sizes

result$measures # access validation measures
}
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{cluster.list}:}{Object of class \code{list} List of all cluster assignments from each algorithm. Use the \code{getClustering} method to conveniently extract clusters for a given value of k. }
    \item{\code{measures}:}{Object of class \code{array} The validation measures calculated }
    \item{\code{clus.methods}:}{Object of class \code{character} The clustering algorithms used. }
    \item{\code{ks}:}{Object of class \code{numeric} The range of cluster numbers tested.  }
    \item{\code{validation}:}{Object of class \code{character} The validation measures used. }
    \item{\code{dist.metric}:}{Object of class \code{character} The distance metric used to calculate validation scores. }
    \item{\code{item.names}:}{Object of class \code{character} The names of the clustered items. }
    \item{\code{call}:}{Object of class \code{call} The call to \code{\link{COMMUNAL}} used to create the object. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{getClustering(k)}:}{ For given value of \code{k}, extract the cluster assignments from each clustering algorithm. }
    \item{\code{show()}:}{ Default print method. }
  }
}
