## COMMUNAL package
## part 1: run COMMUNAL across a range of data

#############  Main function, apply COMMUNAL over a range of data  ###########################
clusterRange <- function(dataMtx, varRange, ...) {
#   The default settings for COMMUNAL:
#   ks = 2:10
#   clus.methods = c("hierarchical", "kmeans"),
#   validation = c("Connectivity", "dunn", "wb.ratio", "g3", "g2", "pearsongamma", "avg.silwidth", "sindex"),
#   dist.metric = "euclidean",
#   aggl.method = "average",
#   neighb.size = 10
  
  if (is.null(rownames(dataMtx))) {
    # add row names for indexing if necessary
    rownames(dataMtx) <- paste("Dim", 1:nrow(dataMtx), sep="")
    cat(paste("Added row names to dataMtx, from Dim1 to Dim", nrow(dataMtx), "\n", sep=""))
  }
  
  args <- list(...)  
  if ( !is.logical(args$reorder) || args$reorder) {
    # by default, sort by variance
    # but does not sort if reorder=FALSE is specified
    data.var <- apply(dataMtx, 1, var, na.rm=TRUE)
    row.order <- order(data.var, decreasing=T)
  } else {
    row.order <- 1:nrow(dataMtx)
  }
  
  genes.list <- list()
  for(i in 1:length(varRange)){
    genes.list[[i]] <- dataMtx[row.order[1:varRange[i]], , drop=FALSE]
  }
  ##  Run Function
  all.results <- lapply(genes.list, function(data) {
    cat(". \t")
    COMMUNAL(data=data, ...)
  })
  
  return(list(all.results=all.results, varRange=varRange))
}