#library(fpc)
#library(Rtsne)
library(flowCore)
library(flowClust)

source_https <- function(u, unlink.tmp.certs = FALSE) {
  # load package
  require(RCurl)
 
  # read script lines from website using a security certificate
  if(!file.exists("cacert.pem")) download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.pem")
  script <- getURL(u, followlocation = TRUE, cainfo = "cacert.pem")
  if(unlink.tmp.certs) unlink("cacert.pem")
 
  # parase lines and evealuate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}

source_https('https://raw.githubusercontent.com/pontikos/FCS/master/fcs.R')


flowCore::read.FCS('beads_4815.fcs')->b

X <- b@exprs[,c('FSC-A','SSC-A')]

# 4 might not always be appropriate
res <- flowClust(X,K=4)

# gate out main cluster of cells on side and forward scatter
smoothPlot(X, classification=res@label, chulls=FALSE)
#cluster with the most events is the main one
i <- which.max(table(res@label))==res@label 
X <- b[i,]@exprs[,grep('-A$',colnames(b))]

X.trans <- apply(X,2,logicleTransform(w=.1))
X.trans <- X.trans[,-grep('SSC-A|FSC-A',colnames(X.trans))]

# i tried SNE dimensionality reduction: didn't look right
#tsne.out <- Rtsne(X.trans)
#smoothPlot(tsne.out$Y,outliers=TRUE)

# use PCA instead on all 18 colours
pca <- princomp(X.trans)
pca.X.trans <- pca$scores[,1:2]

# this didn't work because of very wide cluster
#res <- flowClust(pca.X.trans,K=8)

# dbscan is very slow
#res <- dbscan(pca.X.trans[,1],eps=.01)

# instead we will use a sliding window on the univariate
# density function
# returns the top K sliding window peaks
top.sliding.window.peaks <- function(d, K, span=40) {
  y <- d$y
  x <- d$x
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  peaks <- cbind(x=x[ind],y=y[ind])
  top.peaks <- peaks[order(peaks[,'y'],decreasing=TRUE)[1:K],]
  top.peaks <- top.peaks[order(top.peaks[,'x']),]
  return(top.peaks)
}

# find the top 8 peaks
p <- top.sliding.window.peaks(dens,K=8)[,1]
# looks right?
plot(dens)
abline(v=p)

# cluster using kmeans
res <- kmeans(pca.X.trans[,1],p)
plot(dens)
abline(v=res$centers)

# these are now the MFIs
MFIs <-do.call('rbind',by(X,res$cluster,colMeans))

# plot clustering results 
plotClusters(X.trans[,1:3],classification=res$cluster, chulls=FALSE, ellipses=TRUE)

plotClusters(X.trans[,4:6],classification=res$cluster, chulls=FALSE, ellipses=TRUE)

plotClusters(X.trans[,6:9],classification=res$cluster, chulls=FALSE, ellipses=TRUE)

plotClusters(X.trans[,10:13],classification=res$cluster, chulls=FALSE, ellipses=TRUE)


