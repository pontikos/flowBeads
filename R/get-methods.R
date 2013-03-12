###
setMethod('getDate',
          signature(flow.frame='flowFrame'),
          definition=function(flow.frame) {
            as.Date(flow.frame@description[["$DATE"]],format="%d-%b-%Y")
          } )


###
setMethod('length',
          signature(x='BeadFlowFrame'),
          definition=function(x) {
              return(dim(exprs(x))[1])
          } )

###
setMethod('getParams',
          signature(flow.frame='flowFrame'),
          definition=function(flow.frame) {
              return(grep('^(SSC|FSC)', dimnames(flow.frame@exprs)[[2]], value=T, invert=T))
          } )

###
setMethod('getMEFparams',
          signature(bead.data='BeadFlowFrame'),
          definition=function(bead.data) {
            return( names(bead.data@beads.mef) )
          } )

###
setMethod('hasMEF',
          signature(bead.data='BeadFlowFrame', parameter='character'),
          definition=function(bead.data, parameter) {
            return(parameter %in% names(bead.data@beads.mef))
          } )


