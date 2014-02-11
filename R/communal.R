## COMMUNAL package
## part 1: run clustering algorithms - see clusterRange() for a convenience  harness

makeCommunal <- setRefClass( "COMMUNAL", 
             fields = list( 
               cluster.list = "list",
               measures     = "array",
               clus.methods = "character",
               ks           = "numeric",
               validation   = "character",
               dist.metric  = "character", 
               item.names   = "character",
               call         = "call"
             ), 
             methods = list( 
               getClustering = function(k) {
                 # Returns the cluster assignment of each clustering method
                 #   for the given value of k.
                 # Args:
                 #   k: The desired number of clusters
                 #
                 # Returns:
                 #   Data frame of the cluster assignments, one column for each method.
                 if (!k %in% ks) {
                   stop(paste("k must be one of the original ks: ", paste(ks, collapse=", ")))
                 }
                 clusters <- eval(parse(text=paste("cluster.list$'", k, "'", sep="")))
                 df <- as.data.frame(clusters)
                 data.frame(lapply(df, as.integer), row.names = item.names)
               },
               show = function() {
                 cat("Reference class", classLabel(class(.self)), "\n\n")
                 cat("\tks:", paste(ks, collapse=", "), "\n\n")
                 cat("\tCall: ")
                 methods::show(call)
                 cat("\n(Use the getClustering method to extract data frame of cluster assignments)")
               }
               )
             )

getValidation <- function(clusters, alg, data = NULL, distance = NULL, metric = "euclidean",
                          neighb.size = 10, validation = c("dunn", "avg.silwidth", "Connectivity"), 
                          num.clusters = -1) {
  # Computes the validation scores for a given clustering of the data.
  # Args:
  #   clusters: The cluster assignments
  #   alg: the algorithm name, used to print error messages
  #   data: Data matrix. Must specify either this or distance
  #   distance: Distance matrix
  #   metric: Used to calculate distance matrix
  #   neighb.size: used to calculate connectivity
  #   validation: types of validation metrics to evaluate
  #   num.clusters: the intended number of clusters
  # Returns:
  #   Data frame of the cluster assignments, one column for each method.
  result = numeric()
  metric <- match.arg(metric, c("euclidean", "correlation", "manhattan"))
  if (is.null(distance)) {
    if (is.null(data)) {
      stop("Missing data and distance args. Must specify at least one.")
    }
    distance <- dist(data, method=metric)
  }
  if (all(is.na(clusters))) {
    warning(paste0(alg, ": no cluster assignments available - validation measures will be NA for k=", num.clusters))
    return(rep(NA, length(validation)))
  }
  if ("Connectivity" %in% validation) {
    result <- c(result, clValid::connectivity(distance = distance, clusters = clusters, neighbSize = neighb.size))
  }
  run.silhouette <- "avg.silwidth" %in% validation
  run.g2 <- "g2" %in% validation
  run.g3 <- "g3" %in% validation
  run.sepindex <- "sindex" %in% validation
  stats <- fpc::cluster.stats(d = distance, clustering = clusters, silhouette = run.silhouette, G2 = run.g2, G3 = run.g3, sepindex = run.sepindex)
  for (val in validation) {
    result <- c(result, stats[[val]])
  }
  result
}

addLayer <- function(arr, layer, name) {
  # Adds another layer in the third dimension of the array.
  new.dim <- dim(arr)
  new.dim[3] <- new.dim[3] + 1
  new.names <- dimnames(arr)
  if (dim(arr)[3] == 0) {
    new.names[[3]] <- name
  } else {
    new.names[[3]] <- c(new.names[[3]], name)
  }
  new.array <- array(c(arr,as.vector(layer)), dim = new.dim, dimnames = new.names)
  if (dim(arr)[3] > 0) {
    for (i in 1:dim(arr)[3]) {
      stopifnot(identical(arr[, , i], new.array[, , i]))
    
    }
  }
  new.array
}

COMMUNAL <- function(data,
                     ks = 2:10,
                     clus.methods = c("hierarchical", "kmeans"),
                     validation = c("Connectivity", "dunn", "wb.ratio", "g3", "g2", "pearsongamma", "avg.silwidth", "sindex"),
                     dist.metric = "euclidean",
                     aggl.method = "average",
                     neighb.size = 10, 
                     seed = NULL,
                     ...) {
  # workhorse function to run clustering algorithms
  data <- as.matrix(data)
  if (is.null(colnames(data))) {
    colnames(data) <- paste("Sample", 1:ncol(data), sep="")
    cat(paste("Added column names to data, from Sample1 to Sample", ncol(data), "\n", sep=""))
  }
  
  clus.methods = unique(match.arg(clus.methods, c("hierarchical", "kmeans", "diana", "fanny", "som", 
                                                "model", "sota", "pam", "clara", "agnes", "ccp-hc",
                                                "ccp-km", "ccp-pam", "nmf"), several.ok = TRUE))
  validation = unique(match.arg(validation, c("Connectivity", "average.between", "g2", "ch", "sindex",
                                              "avg.silwidth", "average.within", "dunn",
                                              "widestgap", "wb.ratio", "entropy", "dunn2", "pearsongamma", "g3",
                                              "within.cluster.ss", "min.separation", "max.diameter"), several.ok = TRUE))
  dist.metric = match.arg(dist.metric, c("euclidean", "correlation", "manhattan"))
  aggl.method = match.arg(aggl.method, c("ward", "ward.D", "ward.D2", "single", "complete", "average"))
  cl.methods = intersect(clus.methods, c("hierarchical", "kmeans", "diana", "fanny", "som", 
                                       "model", "sota", "pam", "clara", "agnes"))
  ccp.methods = intersect(clus.methods, c("ccp-hc","ccp-km", "ccp-pam"))
  if ("nmf" %in% clus.methods) {
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop('NMF not installed.')
    }
    if (min(data) < 0) {
      stop("NMF requires all values in 'data' to be positive.")
    }
    
  }
  if (length(ccp.methods) > 0 && !requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    stop(paste("Package ConsensusClusterPlus is required for clustering with", paste(ccp.methods, collapse = ", ")))
  }
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol # from is.integer documentation
  }
  if (!all(sapply(ks, function(x) is.wholenumber(x) & x>1))) {
    stop(paste("'ks' is invalid. Must be whole numbers > 1:", paste(ks, collapse = ", ")))
  }
  if (any(duplicated(ks))) {
    stop("ks is invalid. Must be unique.")
  }
  
  getClus <- function(m, obj, k) {
    # access the clusters for each algorithm
    switch(m,
           "hierarchical" = cutree(obj,k),
           "diana"        = cutree(obj,k),
           "agnes"        = cutree(obj,k),
           "kmeans"       = obj[[as.character(k)]][["cluster"]],
           "clara"        = obj[[as.character(k)]][["clustering"]],
           "pam"          = obj[[as.character(k)]][["clustering"]],
           "fanny"        = obj[[as.character(k)]][["clustering"]],
           "sota"         = obj[[as.character(k)]][["clust"]],
           "model"        = obj[[as.character(k)]][["classification"]],
           "som"          = obj[[as.character(k)]][["unit.classif"]]
    )
  }
  
  cluster.list <- replicate(length(ks), list())
  names(cluster.list) <- as.character(ks)
  mes <- array(dim = c(length(validation), length(ks), 0), dimnames = list(validation, ks, c()))
  
  extractArgs <- function(params, valid.params) {
    # params is a list of parameters (...)
    # valid.params is a vector with valid param names, or example parameter list
    # Returns: subset of params that match valid.params
    if (is.list(valid.params)) {
      valid.params <- names(valid.params)
    }
    params[names(params) %in% valid.params]
  }
  
  # clValid ------------------------
  t.data <- t(data)  # clValid clusters rows, but we want to cluster columns of the input data
  if(dist.metric == "correlation") {
    # as in clValid-functions.R
    distance <- as.dist(1 - cor(data, use = "pairwise.complete.obs"))
  } else{
    distance <- dist(t.data, method = dist.metric)
  }
  if (length(cl.methods) > 0) {
    addl <- extractArgs(list(...), formals(clValid))
    if (length(addl) > 0 ) {
      cl <- clValid(t.data, ks, clMethods = cl.methods, validation = "internal", 
                    metric = dist.metric, method = aggl.method, neighbSize = neighb.size, addl)
    } else {
      cl <- clValid(t.data, ks, clMethods = cl.methods, validation = "internal", 
                    metric = dist.metric, method = aggl.method, neighbSize = neighb.size)
    }
    
    cobjs <- cl@clusterObjs
    for (k in ks) {
      k.clusts <- replicate(length(cl.methods), list())
      names(k.clusts) <- cl.methods
      for (m in cl.methods) {
        k.clusts[[m]] <- getClus(m, cobjs[[m]], k)
      }
      cluster.list[[as.character(k)]] <- k.clusts
    }
    for (alg in cl.methods) { # Compute validation measures.
      layer <- matrix(NA, nrow=dim(mes)[1], ncol = dim(mes)[2])
      for (i in 1:length(ks)) {
        k <- ks[i]
        clus <- cluster.list[[as.character(k)]][[alg]]
        layer[, i] <- getValidation(clusters = clus, distance = distance, metric = dist.metric, 
                                    validation = validation, neighb.size = neighb.size, num.clusters = k, alg=alg)
      }
      mes <- addLayer(mes, layer, name = alg)
    }
  }

  # ConsensusClusterPlus ------------
  if (length(ccp.methods) > 0) {    
    if (dist.metric == "manhattan") {
      warning("For ConsensusClusterPlus, 'manhattan' distance metric not available.
              Switching to 'euclidean' instead.")
    }
    d.metric <- switch(dist.metric, "correlation" = "pearson", "manhattan" = "euclidean", "euclidean" = "euclidean")
    addl <- extractArgs(list(...), formals(ConsensusClusterPlus::ConsensusClusterPlus))
    if (!'reps' %in% names(addl)){
      addl$reps <- 1000 # default number of reps
    }
    for (alg in ccp.methods) {
      alg1 <- sub("ccp-", "", alg)
      layer <- matrix(NA, nrow=dim(mes)[1], ncol = dim(mes)[2])
      
      rcc <- do.call(ConsensusClusterPlus::ConsensusClusterPlus, c(list(data, maxK = max(ks), clusterAlg = alg1, distance = d.metric, seed = seed), addl))
      for (i in 1:length(ks)) {
        k = ks[i]
        clus <- rcc[[k]][["consensusClass"]]
        cluster.list[[as.character(k)]][[alg]] <- clus
        layer[, i] <- getValidation(clusters = clus, distance = distance, metric = dist.metric, 
                                    validation = validation, neighb.size = neighb.size, num.clusters = k, alg=alg)
      }
      mes <- addLayer(mes, layer, name = alg)
    }
  }
  
  # NMF ---------------------------
  if ("nmf" %in% clus.methods) {
    layer <- matrix(NA, nrow = dim(mes)[1], ncol = dim(mes)[2])
    addl <- extractArgs(list(...), formals(NMF::nmf))
    for (i in 1:length(ks)) {
      k = ks[i]
      if (length(addl) > 0) {
        nmfRes <- NMF::nmf(data, rank = k, seed = seed, .options = "v", addl)
      } else {
        nmfRes <- NMF::nmf(data, rank = k, seed = seed, .options = "v")
      }
      NMFsilhouette <- silhouette(nmfRes)
      clus <- NMFsilhouette[,1]
      cluster.list[[as.character(k)]][["nmf"]] <- clus
      layer[,i] <- getValidation(clusters = clus, distance = distance, metric = dist.metric, 
                                 validation = validation, neighb.size = neighb.size,
                                 num.clusters = k, alg="NMF")
    }
    mes <- addLayer(mes, layer, name = "nmf")
  }

  makeCommunal$new(
                  cluster.list = cluster.list,
                  measures     = mes,
                  clus.methods = clus.methods,
                  ks           = ks,
                  validation   = validation,
                  dist.metric  = dist.metric, 
                  item.names   = colnames(data),
                  call         = match.call()
  )
}