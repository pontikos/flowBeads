\name{flowBeads-package}
\alias{flowBeads-package}
\alias{flowBeads}
\docType{package}
\title{
flowBeads: Analysis of Bead data for Flow Cytometry
}
\description{
    Provides additional functionality to the flowCore package to deal with normalisation beads in flow cytometry.
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
b <- BeadFlowFrame(fcs.filename='~nikolas/CAD54_2008FEB11_TREG_6BEADS_110208.FCS')
gb <- gateBeads(b, scatter.gate=T, gating.parameter='APC')
}
