
#' BeadFlowFrame
#'
#' 
#' Extension of \code{\link{flowFrame}} specific for bead data.
#' 
#' 
#' @section Slots:
#'  \describe{
#'   \item{\code{fcs.filename}:}{The file name of the FCS file from which to read.}
#'   \item{\code{bead.filename}:}{The file name of the bead config file.}
#'   \item{\code{beads.mef}:}{The \code{\link{data.frame}} containing the MEF of the bead populations on different channels.}
#'   \item{\code{trans}:}{The transform \eqn{f} to linearise the fluorescence.}
#'   \item{\code{inv.trans}:}{The inverse transform of \eqn{f^{-1}}.}
#' }
#' 
#' @export
#' @name BeadFlowFrame-class
#' @aliases dBeadFlowFrame
#' @rdname BeadFlowFrame-class
#' @docType class
#' @import flowCore
setClass(
         "BeadFlowFrame",
         contains='flowFrame',
         representation(
                        fcs.filename='character',
                        bead.filename='character',
                        beads.mef='data.frame',
                        trans='transform',
                        inv.trans='transform'
                        ),
         prototype(
                   bead.filename='',
                   trans=new('transform'),
                   inv.trans=new('transform')
                   )
         )

#' The constructor take as arguments the FCS file and the file containing the MEF values of the beads on the different detector channels
#' 
#' @param fcs.filename The file name of the FCS to load.  File is loaded with the \code{\link[flowCore:flowCore-package]{read.FCS}} function.
#' @param bead.filename The file name of the MEF configuration files indicating the type of beads in the FCS file. The bead.file is read with \link{read.csv}.
#' @export
#' @rdname BeadFlowFrame-class
BeadFlowFrame <- function(fcs.filename, bead.filename) {
  #if no name is given will try to guess the number of bead populations
    if (missing(bead.filename)) {
        #data(package='flowBeads', dakomef)
        #beads.mef <- dakomef
        bead.filename <- ''
        beads.mef<-data.frame()
    } else {
        beads.mef <- read.csv2(bead.filename)
        colnames(beads.mef) <- toupper(colnames(beads.mef))
    }
    #only keep columns which finish with .A (ignore .W and .H)
    flow.frame <- read.FCS(fcs.filename, alter.names=TRUE, column.pattern=".A")
    #rename the parameters
    colnames(flow.frame@exprs) <- gsub('.A$','', toupper(colnames(flow.frame@exprs)))
    flow.frame@parameters@data['name'] <- colnames(flow.frame@exprs)
    #trans <- gated.bead.data@trans
    if (flow.frame@description$FCSversion == '3') {
        trans <- logicleTransform()
        inv.trans <- inverseLogicleTransform(trans=trans)
    } else if (flow.frame@description$FCSversion == '2') {
        trans <- function (transformationId = "Log[10]") { k <- new('transform', .Data=function(x) x<-log10(x)); k@transformationId <- transformationId; k }
        trans <- trans()
        inv.trans <- function (transformationId = "10^") { k <- new('transform', .Data=function(x) x<-10**x); k@transformationId <- transformationId; k }
        inv.trans <- inv.trans()
    }
    new(
        'BeadFlowFrame',
        fcs.filename=fcs.filename,
        bead.filename=bead.filename,
        beads.mef=beads.mef,
        exprs=flow.frame@exprs,
        parameters=flow.frame@parameters,
        description=flow.frame@description,
        trans=trans,
        inv.trans=inv.trans
    ) 
}


#' GatedBeadFlowFrame
#' 
#' @param labels The resulting labels of the clustering assigning each event to a different bead population.
#' @param clustering.stats Three dimensional array summarising the stats per channel and population.
#' @param mef.tranform The list of MEF transforms
#' @export
#' @name GatedBeadFlowFrame-class
#' @aliases GatedBeadFlowFrame
#' @rdname GatedBeadFlowFrame-class
#' @docType class
setClass(
         "GatedBeadFlowFrame",
         contains='BeadFlowFrame',
         representation(
                        labels="factor",
                        clustering.stats='array',
                        mef.transform='list'
                        )
         )


#' Logicle transformation constructor
#' 
#' Input parameters are to be provided in decades
#' @param transformationId The name of the transformation.
#' @param alpha The intercept of the MEF transform.
#' @param beta The slope of the MEF transform.
#' @export
mefTransform <- function(transformationId="mefTransform", alpha, beta) {
    k <- new("transform", .Data=function(x) x <- beta*x + alpha)            
    k@transformationId <- transformationId
    k
}


