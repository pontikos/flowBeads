flowBeads
=========

This is the public repository for the Bioconductor package for working with calibration beads in flow cytometry.
Based on flowCore package.

This package has received quite a lot of recent attention, in particular I have had a few questions about relative normalisation when no MEF is specified.
I think it is possible provided that the number of peaks is consistent across samples.
However I believe it easier to achieve this without using the package which can be quite rigid.

Here is some R code of how to go about it:

```R
require(flowCore)
# We will use
require(cluster)
# This a very useful Bioconductor package to do clustering by fitting a mixture of normal distributions.
# It also ignores outliers.
require(flowClust)
require(RCurl)
# This is another repository of mine containing some plotting functions for flow data.
download.file("https://raw.githubusercontent.com/pontikos/FCS/master/fcs.R", destfile = "fcs.R", method = "curl")
source('fcs.R')
```

```R
get.MFI <- function(X) {
  print(length(scatter.channels <- as.character(grep('FSC|SSC',colnames(X),value=TRUE))))
  print(length(fluo.channels <- as.character(grep('^SSC|FSC|Time',colnames(X),invert=T,value=T))))
  # use PCA to reduce the number of dimensions from 6 to 2
  pca <- princomp(X[,scatter.channels])
  pca.X <- pca$scores[,1:2]
  #smoothPlot(pca.X)
  res <- flowClust(pca.X,K=4)
  #plot to check result
  #plotClustRes(pca.X,res=res,outliers=FALSE)
  X <- X[which(res@label==which.max(res@w)),]
  X.trans <- apply(X[,fluo.channels],2,logicleTransform())
  res <- pam(X.trans[,1:8],8)
  #check results
  plotClusters(X.trans[,fluo.channels[1:4]],outliers=TRUE, classification=res$clustering,chulls=FALSE) 
  plotClusters(X.trans[,fluo.channels[4:8]],outliers=TRUE,classification=res$clustering,chulls=FALSE)
  # these are now the MFIs
  MFIs <-do.call('rbind',by(X[,fluo.channels],res$clustering,colMeans))
  MFIs <- MFIs[order(MFIs[,1]),]
  return(MFIs)
}
```
 
```R
beads1 <- flowCore::read.FCS(file.path("lucas","QC8PeaksBeads_ After Capture Beads_SAS_ARIAIII_CORDOBA_19112014.fcs"))
beads2 <- flowCore::read.FCS(file.path("lucas","QC8PeaksBeads_32140219_SAS_ARIA_19MAR2015_19MAR2015.fcs"))
MFI.beads1 <- get.MFI(beads1@exprs)
MFI.beads2 <- get.MFI(beads2@exprs)
```

Next it's simply a question of defining the linear transform which maps the peaks between samples and applying that transform to each channel.

```R
trans <-do.call('rbind', lapply( fluo.channels, function(chan) coefficients(lm(logicleTransform()(MFI.beads1[,chan])  ~ logicleTransform()(MFI.beads2[,chan]))) ) )
rownames(trans) <- fluo.channels
colnames(trans) <- c('a','b')
```

Now trans contains the transform to compare samples analysed on the same day as ```beads2``` with those analysed at the same time as ```beads1```.
Note that the transforms are applied on the data after the ```logicleTransform``` so you can use the ```inverseLogicleTransform``` to get back the raw data after having applied the transform.
It might simpler to use a ```log10``` transform instead... more soon

(to be continued...)







