## COMMUNAL package
## part 2: plot the cluster metrics

################   Functions   ####################################
#Here each measure is scaled across ALL data_range, ks, and algs, meaning the final plot
#is an overall measure

getScaledAverages <- function(data_range, geneRange, algs, metrics, ks){
  stopifnot(length(data_range)==length(geneRange)) #Number of cluster runs does not match number of variables
  lowerIsBetter <- c("widestgap", "wb.ratio", "max.diameter", "within.cluster.ss",
                     "average.within", "Connectivity", "entropy", "g3")
  
  #4D array to store all measures; will allow scaling
  mult_array <- array(data=NA, dim=c(length(metrics), length(data_range), length(ks), length(algs)),
                              dimnames=list(metrics, geneRange, ks, algs) )
  
  for(j in 1:length(data_range)){
    #get measures for one data_range instance
    measures <- data_range[[j]]$measures[metrics,as.character(ks),algs,drop=F]
    #make sure all measures in group are plotted  so that higher is better
    for (i in 1:length(metrics)) { 
      flip <- ifelse(metrics[i] %in% lowerIsBetter, -1, 1)
      mult_array[i,j, ,] <- flip*as.numeric(measures[i,,]) 
    }
  }; rm(i,j)
  
  mult_array_scaled <- array(NA, dim=dim(mult_array))
  #scale ALL values for single measure (across data range) together
  for(i in 1:length(metrics)){
    mult_array_scaled[i,,,] <- array(data=(mult_array[i,,,]-mean(mult_array[i,,,], na.rm=T))/sd(mult_array[i,,,], na.rm=T), 
                                     dim=dim(mult_array[i,,,])) 
  }
  dimnames(mult_array_scaled) <- dimnames(mult_array)

  #take average over all algs and measures
  average <- data.frame(matrix(NA, nrow=length(ks), ncol=length(data_range)))  
  rownames(average) <- ks; colnames(average) <- geneRange  
  for(j in 1:length(data_range)){
    ## if only one metric passed, then the first dim of array collapes, so need to decrement index
    ## if only one alg passed, 
    ## if only one of each passed, skip apply() and just take mean
    if(length(metrics) > 1) {
      average[,j] <- apply(mult_array_scaled[,j,,], 2, mean, na.rm=T)
    } else if (length(metrics)==1 & length(algs) > 1) {
      average[,j] <- apply(mult_array_scaled[,j,,], 1, mean, na.rm=T)
    } else if (length(metrics)==1 & length(algs) == 1) {
      average[,j] <- mult_array_scaled[,j,,]
    }
  }
  if(any(is.na(average))) {print(average); stop( "NAs in averages matrix (see above). Decrease range of K.")}
  
  return(average)
}


# Find the highest non-edge 'peak' in a vector (most concave point)
#  or, if no local 'peak', return highest value
findHighestPeak <- function(x){
  diff1 <- diff(diff(x))
  diff1[diff1>=0] <- NA
  diff2 <- which.min(diff1)+1 
  if(length(diff2)==0) {
    return(c(which.max(x), x[which.max(x)]))
  } else {
    return(c(diff2, x[diff2]))
  }
}

plotRange <- function(averageMtx, geneRange, ks, filename){
  geneRangeLabels <- factor(unlist(geneRange))
  names(geneRangeLabels) <- unlist(geneRange)
  rgl::open3d()
  for(i in 1:length(geneRangeLabels)){
    rgl::lines3d(x=ks, y=rep(geneRangeLabels[i], length(ks)), z=averageMtx[,i], lwd=3)
  }
  for(i in seq(min(ks), max(ks), 2)){
    rgl::lines3d(x=c(i,i), y=c(1, length(geneRangeLabels)), z=rep(min(averageMtx),2), lwd=0.5)
  }
  peaks <- apply(averageMtx, 2, findHighestPeak)
  rgl::points3d(x=ks[peaks[1,]], y=geneRangeLabels, z=peaks[2,], col=2, size=10)
  
  z <- data.matrix(averageMtx)
  zlim <- range(z)*2
  zlen <- zlim[2] - zlim[1] + 1
  colorlut <- rainbow(zlen,alpha=0) # height color lookup table
  col <- colorlut[ z*2-zlim[1]+1 ] # assign colors to heights for each point
  rgl::surface3d(x=ks, y=1:length(geneRangeLabels), z=z, color=col, alpha=0.5)
  
  rgl::box3d()
  rgl::mtext3d(text= "K", edge="x--", line=2)
  rgl::mtext3d(text= "Genes", edge="y-+", line=3)
  rgl::mtext3d(text= "measures", edge="z++", line=3)
  rgl::axis3d(edge="x--", at=ks, labels=ks);
  rgl::axis3d(edge="y-+", labels=geneRange, at=1:length(geneRange));
  rgl::axis3d(edge="z++", labels=signif(fivenum(z),2), at=signif(fivenum(z),2)); 
  if (!is.null(filename)) {
    rgl::rgl.snapshot(filename)
  }
}

plotMeans <- function(averageMtx, algs, metrics, ...){
  plot(apply(averageMtx, 1, mean), main=paste("Mean across all geneRange \n", paste(algs, collapse=", "), "\n", paste(metrics, collapse=", ")), 
        xlab="K", xaxt='n', ylab="Combined metric Z-score", type="b",  ...)
  axis(1, at=1:dim(averageMtx)[1], labels=rownames(averageMtx) )
}

#Wrapper
plotRange3D <- function(test_range, ks, goodAlgs, goodMeasures, filename=NULL, ...){
  average <- getScaledAverages(data_range=test_range[[1]], geneRange=test_range[[2]], 
                               algs=goodAlgs, metrics=goodMeasures, ks=ks)
  if (ncol(average) > 1) {
    plotRange(averageMtx=average, test_range[[2]], ks, filename)
  }
  plotMeans(averageMtx=average, algs=goodAlgs, metrics=goodMeasures, ...)
  average
}