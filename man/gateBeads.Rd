\docType{methods}
\name{gateBeads}
\alias{gateBeads}
\title{gateBeads}
\arguments{
  \item{\code{bead.data}}{}

  \item{\code{verbose}}{}
}
\value{
  \link{GatedBeadFlowFrame}
}
\description{
  \code{gateBeads} gates on all channels, apply scatter
  gate first. Find parameters in MEF data.frame which are
  also present in BeadFlowFrame Use k-medoids to identify
  bead populations, the number of bead populations expected
  depends on the type of beads
}
\examples{
data(beads1)
gateBeads(beads1, verbose=T)
}

