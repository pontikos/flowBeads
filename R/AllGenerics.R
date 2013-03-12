# get-methods.R
setGeneric(name='length')

setGeneric(
           name='getDate',
           def=function(flow.frame) standardGeneric('getDate')
           )

setGeneric(
           name='getParams',
           def=function(flow.frame) standardGeneric('getParams')
           )

setGeneric(
           name='getMEFparams',
           def=function(bead.data) standardGeneric('getMEFparams')
           )

setGeneric(
          name='hasMEF',
          def=function(bead.data, parameter) standardGeneric('hasMEF')
         )

# plot-methods.R
setGeneric(name='plot')

# beads.R
setGeneric(
           name='gateBeads',
           def=function(bead.data, ...) standardGeneric('gateBeads')
           )

setGeneric(
           name='toMEF',
           def=function(bead.data, flow.data) standardGeneric('toMEF')
           )

setGeneric(
           name='absoluteNormalise',
           def=function(bead.data, mef.data) standardGeneric('absoluteNormalise')
           )

setGeneric(
           name='relativeNormalise',
           def=function(bead.data1, bead.data2) standardGeneric('relativeNormalise')
           )

setGeneric(
           name='generateReport',
           def=function(bead.data, output.file, ...) standardGeneric('generateReport')
           )

