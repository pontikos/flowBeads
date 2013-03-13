#' BeadFlowFrame
#'
#' 
#' Extension of \code{\link{flowFrame}} specific for bead data.
#' 
#' 
#' @section Slots:
#' \describe{
#'  \item{\code{beads.mef:} lakdf
#'  \item{\code{trans}:} aiodnf
#'  \item{\code{inv.trans}:} adionf
#'  }
#' 
#' @keywords beads flowcytometry
#' @export
#' @docType class
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

#' constructor
#' 
#' @param fcs.filename
#' @param bead.filename
BeadFlowFrame <- function(fcs.filename, bead.filename) {
    if (missing(bead.filename)) {
        data(package='flowBeads', dakomef)
        beads.mef <- dakomef
        bead.filename <- ''
    } else {
        beads.mef <- read.csv(bead.filename)
    }
    colnames(beads.mef) <- toupper(colnames(beads.mef))
    #only keep columns which finish with .A (ignore .W and .H)
    flow.frame <- read.FCS(fcs.filename, alter.names=TRUE, column.pattern=".A")
    #rename the parameters
    colnames(flow.frame@exprs) <- gsub('.A$','', toupper(colnames(flow.frame@exprs)))
    flow.frame@parameters@data['name'] <- colnames(flow.frame@exprs)
    #trans <- gated.bead.data@trans
    if (flow.frame@description$FCSversion == '3') {
        trans <- logicleTransform()
        inv.trans <- inverseLogicleTransform(trans=trans)
    } else if (flow.frame@description$FCSVersion == '2') {
        trans <- new('transform', transformationId='Log[10]', .Data=function(x) log10(x))
        inv.trans <- new('transform', transformationId='10^', .Data=function(x) 10**x)
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
#' @param labels
#' @param clustering
#' @param clustering.stats
#' @param mef.tranform
#' @export
#' @docType class
setClass(
         "GatedBeadFlowFrame",
         contains='BeadFlowFrame',
         representation(
                        labels="factor",
                        clustering='matrix',
                        clustering.stats='array',
                        mef.transform='list'
                        )
         )

#' Logicle transformation constructor
#' 
#' Input parameters are to be provided in decades
#' @param alpha
#' @param beta
mefTransform <- function(transformationId="mefTransform", alpha, beta) {
    k <- new("transform", .Data=function(x) x <- beta*x + alpha)            
    k@transformationId <- transformationId
    k
}


