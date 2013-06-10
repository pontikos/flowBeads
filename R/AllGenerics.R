# get-methods.R

#' @export
#' @docType methods
#' @rdname length-methods
#' @aliases length,flowFrame-method
setGeneric(name='length')

#' @export
#' @docType methods
#' @rdname getDate-methods
#' @aliases getDate,flowFrame-method
setGeneric(
           name='getDate',
           def=function(flow.frame) standardGeneric('getDate')
           )

#' @export
#' @docType methods
#' @rdname getParams-methods
#' @aliases getParams,flowFrame-method
setGeneric(
           name='getParams',
           def=function(flow.frame) standardGeneric('getParams')
           )

#' @export
#' @docType methods
#' @rdname getMEFparams-methods
#' @aliases getMEFparams,BeadFlowFrame-method
setGeneric(
           name='getMEFparams',
           def=function(bead.data) standardGeneric('getMEFparams')
           )

#' @export
#' @docType methods
#' @rdname hasMEF-methods
#' @aliases hasMEF,BeadFlowFrame,character-method
setGeneric(
          name='hasMEF',
          def=function(bead.data, parameter) standardGeneric('hasMEF')
         )

#' @export
#' @docType methods
#' @rdname getClusteringStats-methods
#' @aliases getClusteringStats,GatedBeadFlowFrame-method
setGeneric(
          name='getClusteringStats',
          def=function(bead.data) standardGeneric('getClusteringStats')
         )

#' @export
#' @docType methods
#' @rdname getMEFtransform-methods
#' @aliases getMEFtransform,GatedBeadFlowFrame-method
setGeneric(
          name='getMEFtransform',
          def=function(bead.data) standardGeneric('getMEFtransform')
         )

# plot-methods.R
#' @export
#' @docType methods
#' @rdname plot-methods
setGeneric(name='plot')

# beads.R

#' @export
#' @docType methods
#' @rdname gateBeads-methods
#' @aliases gateBeads,BeadFlowFrame-method
setGeneric(
           name='gateBeads',
           def=function(bead.data, ...) standardGeneric('gateBeads')
           )

#' @export
#' @docType methods
#' @rdname toMEF-methods
#' @aliases toMEF,GatedBeadFlowFrame,flowFrame-method
setGeneric(
           name='toMEF',
           def=function(bead.data, flow.data) standardGeneric('toMEF')
           )

#' @export
#' @docType methods
#' @rdname absoluteNormalise-methods
#' @aliases absoluteNormalise,GatedBeadFlowFrame,data.frame-method
setGeneric(
           name='absoluteNormalise',
           def=function(bead.data, mef.data) standardGeneric('absoluteNormalise')
           )

#' @export
#' @docType methods
#' @rdname relativeNormalise-methods
#' @aliases relativeNormalise,GatedBeadFlowFrame,GatedBeadFlowFrame-method
setGeneric(
           name='relativeNormalise',
           def=function(bead.data1, bead.data2) standardGeneric('relativeNormalise')
           )

#' @export
#' @docType methods
#' @rdname generateReport-methods
#' @aliases generateReport,GatedBeadFlowFrame,character-method
setGeneric(
           name='generateReport',
           def=function(bead.data, output.file, ...) standardGeneric('generateReport')
           )

