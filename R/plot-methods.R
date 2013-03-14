
#' Plot the results of the clustering.
#' Plot only the requested channel
#' which should have a corresponding entry in the MEF files
setMethod('plot',
    signature=signature(x='GatedBeadFlowFrame',y='character'),
    definition=function(x, y, ...) {
        bead.data <- x
        bead.parameters <- y
        trans <- bead.data@trans
        trans.name <- bead.data@trans@transformationId

        old.par <- par(no.readonly=T)
        par(mfrow=c(length(bead.parameters), 1), oma=c(0,0,4,0), mar=c(4,3,2,0), mgp=c(2, 1, 0))

        for (p in bead.parameters) {
            labels <- bead.data@clustering[,p]
            x <- trans(bead.data@exprs[,p])
            mfi <- trans(bead.data@clustering.stats['mean.fi',p,])[-1]
            h <- hist(x, breaks=1000, plot=F)
            col.mids <- list()
            col.mids[as.character(h$mids)] <- 0
            col.mids[as.character(h$mids[findInterval(x, h$breaks)])] <- labels
            xlab <- substitute(trans.name(p), list(trans.name=trans.name, p=p) )
            if (hasMEF(bead.data, p)) {
                #mef draw histogram and regression line
                mef <- trans(bead.data@beads.mef[,p])
                ymax <- max(mef)+1
                alpha <- bead.data@mef.transform[[p]][['alpha']]
                beta <- bead.data@mef.transform[[p]][['beta']]
                equation <- substitute( transname(MEF) == b %*% transname(p) + a, list(b=beta, a=alpha, p=paste(p, 'MFI'), transname=trans.name) )
                plot( mfi, mef,
                    xlim=c(0, max(x)),
                    ylim=c(0, ymax),
                    xlab=equation,
                    ylab=substitute(trans.name(MEF), list(trans.name=trans.name)),
                    main=xlab)
                abline(b=beta, a=alpha) 
                #the mfis
                for (i in sort(unique(labels))) {
                    segments(x0=mean(x[labels==i]),y0=0, y1=beta*mean(x[labels==i])+alpha, lty=2)
                }
                h$counts <- ymax*h$counts/max(h$counts)
                #text(2,5, equation, cex=.75)
            } else {
                #no mef just draw histogram
                plot(h,
                     xlim=c(0, max(x)),
                     col='white',
                     xlab='',
                     ylab=NULL,
                     main=xlab)
            }
            #the coloured histogram
            segments(x0=h$mids, y0=0, x1=h$mids, y1=h$counts, col=as.numeric(col.mids[as.character(h$mids)]))
        }
        title( sprintf('Date: %s', getDate(bead.data)), line=1, outer=T)
        title( sprintf('Number of beads: %d', length(bead.data)), line=2, outer=T)
        par(old.par)
    }
)

#'
#' Ungated bead data, simply draw all channels individually (no colours).
#'
#'
setMethod('plot',
          signature=signature(x='BeadFlowFrame',y='character'),
          definition=function(x, y, ...) {
            bead.data <- x
            bead.parameters <- y
            trans <- bead.data@trans
            trans.name <- bead.data@trans@transformationId            
            old.par <- par(no.readonly=T)
            par(mfrow=c(length(bead.parameters), 1), oma=c(0,0,4,0), mar=c(4,3,2,0), mgp=c(2, 1, 0))           
            for (p in bead.parameters) {
              x <- trans(bead.data@exprs[,p])
              h <- hist(x, breaks=1000, plot=F)
              xlab <- substitute(trans.name(p), list(trans.name=trans.name, p=p) )
              plot(h,
                     xlim=c(0, max(x)),
                     col='white',
                     xlab='',
                     ylab=NULL,
                     main=xlab)
              segments(x0=h$mids, y0=0, x1=h$mids, y1=h$counts)
            }
            title( sprintf('Date: %s', getDate(bead.data)), line=1, outer=T)
            title( sprintf('Number of beads: %d', length(bead.data)), line=2, outer=T)
            par(old.par)
          }
)



#'  Plot function for \code{BeadFlowFrame}
#'
#'  If no argument specified then plot all parameters
setMethod('plot',
    signature=signature(x='BeadFlowFrame',y='missing'),
    definition=function(x, ...) {
        plot(x, getParams(x))
    }
)

