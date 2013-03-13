## ==========================================================================
## show and print methods display details about an object, and print usually
## allows for a bit mor fine-grained control.
## ==========================================================================


#' 
#' BeadFlowFrame
#' 
#' @rdname show-methods
setMethod("show",
          signature=signature(object="BeadFlowFrame"),
          definition=function(object)
      {
          dm <- dim(exprs(object))
          cat(paste("BeadFlowFrame object '",
                    identifier(object),
                    "'\nfrom ", getDate(object),
                    "\nwith ", length(object), " beads and ", 
                    dm[2], " observables:\n", sep=""))
          show(pData(parameters(object)))
          cat(paste(length(description(object)), " keywords are stored in the ",
                    "'description' slot\n", sep = ""))
          cat('Beads MEF\n')
          show(object@beads.mef)
          return(invisible(NULL))
      })


#'
#' GatedBeadFlowFrame
#'
#'@rdname show-methods
setMethod("show",
          signature=signature(object="GatedBeadFlowFrame"),
          definition=function(object)
      {
          dm <- dim(exprs(object))
          cat(paste("GatedBeadFlowFrame object '",
                    identifier(object),
                    "'\nfrom ", getDate(object),
                    "\nwith ", length(object), " beads and ", 
                    dm[2], " observables:\n", sep=""))
          show(pData(parameters(object)))
          cat(paste(length(description(object)), " keywords are stored in the ",
                    "'description' slot\n", sep = ""))
          cat('Beads MEF\n')
          show(object@beads.mef)
          return(invisible(NULL))
      })

