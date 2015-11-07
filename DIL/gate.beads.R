#!/usr/bin/env Rscript
#
#single colour calibration
#

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("flowCore"))
suppressPackageStartupMessages(library("flowClust"))


option_list <- list( 
    make_option(c("-f","--fcs"), help = "fcsfile to parse"),
    make_option(c("-p","--parameter"), help="parameter on which to gate"),
    make_option(c("--bead.type"), default="sixbeads-mef.csv", help = ""),
    make_option(c("--report"), help="report")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (is.null(opt$fcs)) {
    stop("No FCS file specified on command line!")
}
if (is.null(opt$parameter)) {
    stop("No gating parameter specified on command line!")
}

#bead.data <- BeadFlowFrame(fcs.filename=opt$fcs)
#gated.bead.data <- gateBeads(bead.data)

d <- read.FCS(opt$fcs)

x <- d@exprs
x <- x[apply(x,1,function(y)sum(y<=0)==0),]

#singlet filtering
#colMeans(x[,c('FSC-A','SSC-A')])

x <- apply(x,2,log10)

#apply(x,2,quantile,probs=seq(0,1,1/6))

param <- colnames(x)[8:15]
pca.x <- princomp(x[,param])

res <- flowClust( x[,param], K=6 )

do.call('rbind', by(x[,param],res@label,colMeans ))


bead.stats <- function(x, labels) {
    res <- list()
    res$N <- tapply(x,labels, length)
    res$mean.fi <- tapply(x,labels, mean)
    res$count <- tapply(x,labels, length)
    res$min.fi <- tapply(x,labels, min)
    res$max.fi <- tapply(x,labels, max)
    res$p95.fi <- tapply(x,labels, function(y) quantile(y,probs=seq(0,1,.05))[['95%']])
    res$mean.fi <- tapply(x,labels, mean)
    res$sd.fi <- tapply(x,labels, sd)
    res$snr <- res$mean.fi/res$sd.fi
    res$cv <- 100*res$sd.fi/res$mean.fi
    res$mean.fi.2 <- res$mean.fi[-1]
    return(res)
}

bead.stats(x[,'APC-A'], res@label)

print(res)

exit()


lm( log10(mef) ~ log10(mean.fi.2) )->m
a=round(m$coefficients[1], digits=3)
b=round(m$coefficients[2], digits=3)
rse=summary(m)$sigma
residuals=m$residuals

if (!is.null(opt$report)) 
    generate.report(output.file=opt$report)

#print header
cat(">>file.name,fcs.version,date,channel,mef,")
cat(
paste("N",1:opt$NumC, sep="."),
paste("mfi",1:opt$NumC, sep="."),
paste("sdfi",1:opt$NumC, sep="."),
paste("minfi",1:opt$NumC, sep="."),
paste("maxfi",1:opt$NumC, sep="."),
paste("p95fi",1:opt$NumC, sep="."),
"mef.alpha","mef.beta","mef.rse",sep=",")
cat("\n")

#print data
cat(sprintf(">>%s,%s,%s,%s,%s,", toupper(file.name), fcs.version, experiment.date, gating.parameter, opt$mef))
cat(b@N,sep=",")
cat(',')
cat(b@mean.fi,sep=",")
cat(',')
cat(b@sd.fi,sep=",")
cat(',')
cat(b@min.fi,sep=",")
cat(',')
cat(b@max.fi,sep=",")
cat(',')
cat(b@p95.fi,sep=",")
cat(sprintf(",%.4f,%.4f,%.4f\n",b@alpha,b@beta,b@rse))

warnings()


