
#' gateBeads 
#' 
#' @description
#' \code{gateBeads} gates on all channels, apply scatter gate first.
#' Find parameters in MEF data.frame which are also present in BeadFlowFrame
#' Use k-medoids to identify bead populations, the number of bead populations
#' expected depends on the type of beads
#' @param \code{bead.data}
#' @param \code{verbose}
#' @return \link{GatedBeadFlowFrame}
#' @export
#' @docType methods
#' @examples
#' data(beads1)
#' gateBeads(beads1, verbose=T)
setMethod('gateBeads',
          signature=signature(bead.data='BeadFlowFrame'),
          definition=function(bead.data, verbose=F, ...) {
        #gating on scatter is important to eliminate doublets
        #remove all negative SSC and FSC
        if (verbose) cat('Removing all negative FSC and SSC...\n')
        bead.data <- bead.data[ bead.data@exprs[,'FSC'] > 0 & bead.data@exprs[,'SSC'] > 0, ]
        #gate on FSC SSC first to retrieve singlets
        if (verbose) cat('Gating on FSC and SSC with norm2Filter...\n')
        bead.data <- Subset(bead.data, filter(bead.data, flowCore::norm2Filter(filterId="ScatterFilter", x=c("FSC", "SSC"), scale.factor=1)))
        K <- dim(bead.data@beads.mef)[1]
        all.parameters <- dimnames(bead.data@exprs)[[2]]
        bead.parameters <- grep('^(SSC|FSC)', all.parameters, value=T, invert=T)
        gated.bead.data <- as(bead.data, 'GatedBeadFlowFrame')
        #gating.parameter <- all.parameters[grep(gating.parameter, all.parameters, ignore.case=T)]
        #for all parameters except forward and side scatter
        gated.bead.data@clustering <- matrix(nrow=nrow(bead.data@exprs), ncol=length(bead.parameters), dimnames=list(NULL, bead.parameters))
        stats.fun <- c('count', 'mean.fi', 'sd.fi', 'cv')
        gated.bead.data@clustering.stats <- array(dim=c(length(stats.fun),length(bead.parameters),K), dimnames=list(stats.fun,bead.parameters))
        #parameters which overlap between MEF file and bead.data
        trans <- bead.data@trans
        inv.trans <- bead.data@inv.trans
        for (gating.parameter in bead.parameters) {
            x <- bead.data@exprs[,gating.parameter]
            labels <- numeric(length(x))
            #for FCS 2 we need to log10 the data before clustering
            if (bead.data@description$FCS == '2') {
                #work with a subset of x which only contains >1
                i <- x>1
                intensity <- x[i]
            #for FCS 3 no need to linearise the data 
            } else if (bead.data@description$FCS == '3') {
                i <- TRUE
                intensity <- x
            }
            if(verbose) cat('Gating on', gating.parameter, 'for', K, 'bead populations with pam...\n')
            res <- cluster::pam(trans(intensity), k=K)
            cluster <- res$clustering
            labels[i] <- cluster
            min.cluster <- cluster[which.min(x[i])] 
            #add back the points that were removed before clustering if any
            labels[!i] <- min.cluster
            #sort the labels
            labels <- factor(labels, levels=names(sort(tapply(x, labels, min))))
            levels(labels) <- 1:K
            #store clustering results in a matrix of same dimension as exprs
            gated.bead.data@clustering[,gating.parameter] <- labels
            #
            gated.bead.data@clustering.stats['count',gating.parameter,] <- as.numeric(tapply(x,labels, length))
            gated.bead.data@clustering.stats['mean.fi',gating.parameter,] <- as.numeric(tapply(x,labels, mean))
            gated.bead.data@clustering.stats['sd.fi',gating.parameter,] <- as.numeric(tapply(x,labels, sd))
            gated.bead.data@clustering.stats['cv',gating.parameter,] <- as.numeric(100*as.numeric(tapply(x,labels, sd))/as.numeric(tapply(x,labels, mean)))
        }
        #we can compute the MEF transform for these
        common.params <- intersect(bead.parameters, getMEFparams(bead.data))
        if (verbose) cat('Common bead fluorochromes and channels:', common.params, '\n')
        for (cp in common.params) {
            mef <- bead.data@beads.mef[,cp]
            mfi <- gated.bead.data@clustering.stats['mean.fi',cp,]
            m <- lm( trans(mef) ~ trans(mfi[-1]) )
            a <- as.numeric(round(m$coefficients[1], digits=3))
            b <- as.numeric(round(m$coefficients[2], digits=3))
            rse <- as.numeric(round(summary(m)$sigma, digits=3))
            f <- function(a,b) {force(a); force(b); function(x) inv.trans(b*trans(x)+a)}
            gated.bead.data@mef.transform[[cp]] <- list(alpha=a, beta=b, m=m, rse=rse, fun=f(a,b))
        }
        return(gated.bead.data)
    } 
)


#' toMEF
#' 
#' @description
#' Given bead.data and a flow.data apply the MEF transform to matching channels in \code{flow.data}.
#' @param bead.data
#' @param flow.data
#' @export
#' @docType methods
setMethod('toMEF',
          signature=signature(bead.data='GatedBeadFlowFrame', flow.data='flowFrame'),
          definition=function(bead.data, flow.data) {
            #standardise file name
            colnames(flow.data@exprs) <- gsub('.A$','', toupper(colnames(flow.data@exprs)))
            flow.data@parameters@data['name'] <- colnames(flow.data@exprs)
            for (p in names(bead.data@mef.transform)) {
              normalise <- bead.data@mef.transform[[p]]$fun
              x <- flow.data@exprs[,p]
              flow.data@exprs[,p] <- normalise(x)
            }
            return(flow.data)
          }
)


#' absoluteNormalise
#' 
#' @description
#' Absolute normalise to align peaks of bead.data to MEF.
#' @return A list of affine functions from transformed MFI relative coordinates to transformed MEF absolute coordinates.
#' @param \code{bead.data} \code{\link{GatedBeadFlowFrame}}
#' @param mef.data \link{data.frame}
#' @export
#' @docType methods
setMethod('absoluteNormalise',
          signature=signature(bead.data='GatedBeadFlowFrame', mef.data='data.frame'),
          definition=function(bead.data, mef.data) {
            params <- intersect(getParams(bead.data), names(mef.data))
            absolute.normalise <- list()
            for (p in params) {
              mfi <- bead.data@trans(bead.data@clustering.stats['mean.fi', p,])
              mef <- bead.data@trans(mef.data[,p])
              m <- lm( mef ~ mfi[-1] )
              a <- as.numeric(round(m$coefficients[1], digits=3))
              b <- as.numeric(round(m$coefficients[2], digits=3))
              rse <- as.numeric(round(summary(m)$sigma, digits=3))
              f <- function(a, b) {force(a); force(b); function(x) b*x+a}
              absolute.normalise[[p]] <- list( mfi=mfi, mef=mef, alpha=a, beta=b, m=m, rse=rse, fun=f(a,b) )
            }
            return(absolute.normalise)
          }
)

#' relativeNormalise
#' 
#' @description
#' Relative normalise to align peaks of bead.data1 to those of bead.data2
#' Returns a list of affine functions from transformed MFI day one coordinates to transformed MFI day two coordinates.
#' This permits comparison of channels across two days, provided the detector is stable, even in the absence of absolute MEF values.
#' @return A list of affine functions from MFI day one coordinates to MFI day two coordinates.
#' @param bead.data1: \code{\link{GatedBeadFlowFrame}} object with MFIs from day one
#' @param bead.data2: \code{\link{GatedBeadFlowFrame}} object with MFIIs from day two
#' @export
#' @docType methods
setMethod('relativeNormalise',
          signature=signature(bead.data1='GatedBeadFlowFrame', bead.data2='GatedBeadFlowFrame'),
          definition=function(bead.data1, bead.data2) {
            params <- intersect(getParams(bead.data1), getParams(bead.data2))
            relative.normalise <- list()
            for (p in params) {
              mfi1 <- bead.data1@trans(bead.data1@clustering.stats['mean.fi', p,])
              mfi2 <- bead.data2@trans(bead.data2@clustering.stats['mean.fi', p,])
              m <- lm( mfi2 ~ mfi1 )
              a <- as.numeric(round(m$coefficients[1], digits=3))
              b <- as.numeric(round(m$coefficients[2], digits=3))
              rse <- as.numeric(round(summary(m)$sigma, digits=3))
              f <- function(a, b) {force(a); force(b); function(x) b*x+a}
              relative.normalise[[p]] <- list( mfi1=mfi1, mfi2=mfi2, alpha=a, beta=b, m=m, rse=rse, fun=f(a,b))
            }
            return(relative.normalise)
          }
)


#' generateReport
#' 
#' Generate an HTML report from a Markdown template using \link{knitr}.
#' @seealso knitr
#' @param bead.data \code{\link{GatedBeadFlowFrame}}
#' @param output.file name of the file to which to output the HTML report.
#' @export
#' @docType methods
setMethod('generateReport',
          signature=signature(bead.data='GatedBeadFlowFrame', output.file='character'),
          definition=function(bead.data, output.file, template=system.file("markdownTemplates/bead-report.Rmd", package = "flowBeads")) {
            #knitr::knit2html(template, options=c("use_xhtml","smartypants","mathjax","highlight_code"), out=output.file)
            knitr::knit2html(template, out=output.file)
          }
)



