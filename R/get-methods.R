#' getDate
#' 
#' @param flow.frame \code{\link{flowFrame}} object on which to get the date field
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
setMethod('length',
          signature(x='flowFrame'),
          definition=function(x) {
              return(dim(exprs(x))[1])
          } )

#' getParams
#' 
#' @description
#' Returns all the parameter names except the scatter channels.
#' 
setMethod('getParams',
          signature(flow.frame='flowFrame'),
          definition=function(flow.frame) {
              return(grep('^(SSC|FSC)', dimnames(flow.frame@exprs)[[2]], value=T, invert=T))
          } )

#' getMEFparams
#' 
#' @description
#' Returns all the MEF parameter names.
#' 
setMethod('getMEFparams',
          signature(bead.data='BeadFlowFrame'),
          definition=function(bead.data) {
            return( names(bead.data@beads.mef) )
          } )

#' hasMEF
#' 
#' @description
#' Checks whether we have the MEF for a channel name.
#' @param bead.data \code{\link{BeadFlowFrame}}
#' @param parameter \code{\link{character}}
setMethod('hasMEF',
          signature(bead.data='BeadFlowFrame', parameter='character'),
          definition=function(bead.data, parameter) {
            return(parameter %in% names(bead.data@beads.mef))
          } )


