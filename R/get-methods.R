
#' getDate
#' 
#' @param flow.frame \code{\link{flowFrame}} object on which to get the date field
#' @rdname getDate-methods
#' @export
#' @docType methods
setMethod('getDate',
          signature(flow.frame='flowFrame'),
          definition=function(flow.frame) {
            as.Date(flow.frame@description[["$DATE"]],format="%d-%b-%Y")
          } )


#' length
#' 
#' @description
#' Returns the number of events in a \code{\link{flowFrame}} object.
#' @param flow.frame \code{\link{flowFrame}} object on which to get number of beads
#' @rdname length-methods
#' @export
#' @docType methods
setMethod('length',
          signature(x='flowFrame'),
          definition=function(x) {
              return(nrow(x))
          } )

#' getParams
#' 
#' @description
#' Returns all the parameter names except the scatter channels.
#' @rdname getParams-methods
#' @export
#' @docType methods
setMethod('getParams',
          signature(flow.frame='flowFrame'),
          definition=function(flow.frame) {
              return(grep('^(SSC|FSC)', colnames(flow.frame), value=T, invert=T))
          } )

#' getMEFparams
#' 
#' @description
#' Returns all the MEF parameter names.
#' @rdname getMEFparams-methods
#' @export
#' @docType methods
setMethod('getMEFparams',
          signature(bead.data='BeadFlowFrame'),
          definition=function(bead.data) {
            return( names(bead.data@beads.mef) )
          } )

#' hasMEF
#' 
#' @description
#' Checks whether we have the MEF for a channel name.
#' @rdname hasMEF-methods
#' @param bead.data \code{\link{BeadFlowFrame}}
#' @param parameter \code{\link{character}}
setMethod('hasMEF',
          signature(bead.data='BeadFlowFrame', parameter='character'),
          definition=function(bead.data, parameter) {
            return(parameter %in% names(bead.data@beads.mef))
          } )

#' getClusteringStats
#' 
#' @description
#' Returns clustering stats as a 3-dimensional array.
#' @rdname getClusteringStats-methods
#' @export
#' @docType methods
setMethod('getClusteringStats',
          signature(bead.data='GatedBeadFlowFrame'),
          definition=function(bead.data) {
            return(bead.data@clustering.stats)
          } )

#' getMEFtransform
#' 
#' @description
#' Returns MEF transform function.
#' @rdname getMEFtransform-methods
#' @export
#' @docType methods
setMethod('getMEFtransform',
          signature(bead.data='GatedBeadFlowFrame'),
          definition=function(bead.data) {
            return(bead.data@mef.transform)
          } )


