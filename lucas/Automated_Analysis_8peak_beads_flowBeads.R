## Automated opencyto ##
# Author: Lucas Le Lann
# Date: 20/05/2015
# 
setwd('lucas')

## loading packages ##
suppressPackageStartupMessages(
{
  library(xtable)	
  library(flowBeads)
  library(flowCore)
  library(flowStats)
  library(flowViz)
  })

## Automated gate ##

beads1 <- BeadFlowFrame("QC8PeaksBeads_ After Capture Beads_SAS_ARIAIII_CORDOBA_19112014.fcs")
beads2 <- BeadFlowFrame("QC8PeaksBeads_32140219_SAS_ARIA_19MAR2015_19MAR2015.fcs")

gbeads1 <- gateBeads(beads1, K = 8)
gbeads2 <- gateBeads(beads2, K = 8)

plot(gbeads1)


relative.transforms <- relativeNormalise(gbeads1, gbeads2)
names(relative.transforms)
name <- names(relative.transforms)
relative.transforms

## FITC ##
fun <- relative.transforms$FITC$fun
mfi1 <- getTransformFunction(gbeads1)(getClusteringStats(gbeads1)['mean.fi','FITC',])
mfi2 <- getTransformFunction(gbeads2)(getClusteringStats(gbeads2)['mean.fi','FITC',])
fun.mfi1 <- fun(mfi1)
fun.mfi2 <- fun(mfi2)

par(mfrow = c(2, 1))
plot(mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

plot(fun.mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

## PE ##
fun <- relative.transforms$PE$fun
mfi1 <- getTransformFunction(gbeads1)(getClusteringStats(gbeads1)['mean.fi','PE',])
mfi2 <- getTransformFunction(gbeads2)(getClusteringStats(gbeads2)['mean.fi','PE',])
fun.mfi1 <- fun(mfi1)
fun.mfi2 <- fun(mfi2)

par(mfrow = c(2, 1))
plot(mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

plot(fun.mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

## PC7 ##
fun <- relative.transforms$PC7$fun
mfi1 <- getTransformFunction(gbeads1)(getClusteringStats(gbeads1)['mean.fi','PC7',])
mfi2 <- getTransformFunction(gbeads2)(getClusteringStats(gbeads2)['mean.fi','PC7',])
fun.mfi1 <- fun(mfi1)
fun.mfi2 <- fun(mfi2)

par(mfrow = c(2, 1))
plot(mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

plot(fun.mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

## PC5.5 ##
fun <- relative.transforms$PC5.5$fun
mfi1 <- getTransformFunction(gbeads1)(getClusteringStats(gbeads1)['mean.fi','PC5.5',])
mfi2 <- getTransformFunction(gbeads2)(getClusteringStats(gbeads2)['mean.fi','PC5.5',])
fun.mfi1 <- fun(mfi1)
fun.mfi2 <- fun(mfi2)

par(mfrow = c(2, 1))
plot(mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

plot(fun.mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)


## APC ##
fun <- relative.transforms$APC$fun
mfi1 <- getTransformFunction(gbeads1)(getClusteringStats(gbeads1)['mean.fi','APC',])
mfi2 <- getTransformFunction(gbeads2)(getClusteringStats(gbeads2)['mean.fi','APC',])
fun.mfi1 <- fun(mfi1)
fun.mfi2 <- fun(mfi2)

par(mfrow = c(2, 1))
plot(mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

plot(fun.mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

## APC.AF750 ##

fun <- relative.transforms$APC.AF750$fun
mfi1 <- getTransformFunction(gbeads1)(getClusteringStats(gbeads1)['mean.fi','APC.AF750',])
mfi2 <- getTransformFunction(gbeads2)(getClusteringStats(gbeads2)['mean.fi','APC.AF750',])
fun.mfi1 <- fun(mfi1)
fun.mfi2 <- fun(mfi2)

par(mfrow = c(2, 1))
plot(mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

plot(fun.mfi1, mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

plot(fun.mfi2, mfi1,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)
plot(mfi2, fun.mfi1,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

plot(mfi1, fun.mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)
plot(fun.mfi1, fun.mfi2,xlim=c(0,5),ylim=c(0,5))
abline(b=1,a=0)

mfi1
fun.mfi1
mfi2
fun.mfi2

#generateReport(gbeads1, output.file="report.html")
