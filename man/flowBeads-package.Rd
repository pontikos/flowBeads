\name{flowBeads-package}
\alias{flowBeads-package}
\alias{flowBeads}
\docType{package}
\title{
flowBeads: Analysis of Bead data for Flow Cytometry
}
\description{
    Provides additional functionality to the \code{\link{flowCore}} package to deal with normalisation beads in flow cytometry.
}
\details{
\tabular{ll}{
Package: \tab flowBeads \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-02-28\cr
License: \tab Artistic-2.0\cr
Depends: \tab methods\cr
}
Beads are routinely ran in flow cytometry for instrument calibration and for normalisation purposes.
}
\author{
Author and Maintainer: Nikolas Pontikos <n.pontikos@gmail.com>
}

\keyword{ package }
\seealso{
\code{\link[flowCore:flowCore-package]{flowCore}}
}
\examples{
library(flowBeads)
b <- BeadFlowFrame(fcs.filename=system.file('extdata', 'beads1.fcs', package='flowBeads'))
gb <- gateBeads(b)
}
