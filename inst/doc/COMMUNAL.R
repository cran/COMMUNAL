### R code from vignette source 'COMMUNAL.Rnw'

###################################################
### code chunk number 1: COMMUNAL.Rnw:27-36
###################################################
## create artificial data set with 3 distinct clusters in two dimensions
set.seed(1)
V1 = c(abs(rnorm(20, 2, 20)), abs(rnorm(20, 65, 15)), abs(rnorm(20, 140, 20)))
V2 = c(abs(rnorm(20, 2, 20)), abs(rnorm(20, 65, 15)), abs(rnorm(20, 105, 20)))
data <- t(data.frame(V1, V2))
colnames(data) <- paste("Sample", 1:ncol(data), sep="")
rownames(data) <- paste("Gene", 1:nrow(data), sep="")
plot(V1, V2, col=rep(c("red", "blue", "black"), each=20), pch=rep(c(0,1,2), each=100),
     xlab="x", ylab="y")


###################################################
### code chunk number 2: COMMUNAL.Rnw:46-50
###################################################
## run COMMUNAL with defaults
library(COMMUNAL)
ks <- seq(2,5)
result <- COMMUNAL(data=data, ks=ks)


###################################################
### code chunk number 3: COMMUNAL.Rnw:54-55
###################################################
result


###################################################
### code chunk number 4: COMMUNAL.Rnw:59-60
###################################################
str(result$measures)


###################################################
### code chunk number 5: COMMUNAL.Rnw:68-72
###################################################
result.list <- list(list(result), ncol(data))
goodAlgs <- c('hierarchical', 'kmeans')
goodMeasures <- c('wb.ratio', 'avg.silwidth', 'dunn')
values <- plotRange3D(result.list, ks, goodAlgs, goodMeasures)


###################################################
### code chunk number 6: COMMUNAL.Rnw:80-82
###################################################
clusters <- result$getClustering(k=3)
table(clusters)


###################################################
### code chunk number 7: COMMUNAL.Rnw:86-92
###################################################
# find 'core' clusters
mat.key <- clusterKeys(clusters, k=3)
examineCounts(mat.key)
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters
table(core) # the 'core' clusters
head(core) # the cluster assignments


###################################################
### code chunk number 8: COMMUNAL.Rnw:96-102
###################################################
k <- 3
clusters <- data.frame(
  alg1=as.integer(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,1)),
  alg2=as.integer(c(1,1,1,1,1,3,3,3,3,3,2,2,2,2,1)),
  alg3=as.integer(c(3,3,3,3,3,1,1,1,1,1,2,2,2,2,2))
)


###################################################
### code chunk number 9: COMMUNAL.Rnw:106-108
###################################################
mat.key <- clusterKeys(clusters, k)
mat.key # cluster indices are relabeled


###################################################
### code chunk number 10: COMMUNAL.Rnw:112-113
###################################################
examineCounts(mat.key)


###################################################
### code chunk number 11: COMMUNAL.Rnw:117-119
###################################################
core <- returnCore(mat.key, agreement.thresh=50) # find 'core' clusters
table(core) # the 'core' clusters


###################################################
### code chunk number 12: COMMUNAL.Rnw:123-125
###################################################
core <- returnCore(mat.key, agreement.thresh=99)
table(core) # 0 is undetermined


###################################################
### code chunk number 13: COMMUNAL.Rnw:134-143
###################################################
data(BRCA.100)
algs <- c("hierarchical", "kmeans", "agnes")
measures <- c('wb.ratio', 'dunn', 'avg.silwidth')
varRange <- c(50, 70, 85, 100)
ks <- 2:5
range.results <- clusterRange(dataMtx=BRCA.100, varRange=varRange,
                              ks = ks,
                              clus.methods = algs,
                              validation = measures)


###################################################
### code chunk number 14: COMMUNAL.Rnw:148-149
###################################################
plot.data <- plotRange3D(range.results, ks, algs, measures, filename='snapshot.png')


