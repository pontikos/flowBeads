\docType{methods}
\name{relativeNormalise}
\alias{relativeNormalise}
\title{relativeNormalise}
\arguments{
  \item{bead.data1}{\code{\link{GatedBeadFlowFrame}}}

  \item{bead.data2}{\code{\link{GatedBeadFlowFrame}}}
}
\description{
  Relative normalise to align peaks of bead.data1 to those
  of bead.data2 Returns a list of affine functions from
  transformed MFI day 1 coordinates to transformed MFI day
  2 coordinates. This permits comparison of channels across
  two days, provided the detector is stable, even in the
  absence of absolute MEF values.
}

