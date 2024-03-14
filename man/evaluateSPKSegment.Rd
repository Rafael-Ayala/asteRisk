\name{evaluateSPKSegment}
\alias{evaluateSPKSegment}
\title{Evaluate a given SPK segment
}
\description{
SPK (Spacecraft and Planet Kernel) is a binary file format developed by NAIF to 
store ephemerides (trajectory) of celestial bodies and spacecraft in the Solar 
System. A detailed description of the SPK file format can be found in NAIF's 
documentation (https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html).

Each SPK file contains several segments, with each segment comprising a summary
(or descriptor), a name and an array of double precision elements. Each segment
is conceptually equivalent to an array in the context of generic DAF files.
There are several types of SPK segments defined by NAIF, each identified by an SPK type
code currently ranging from 1 to 21 (some intermediate values are not used or not
available for general public use). Each segment type provides ephemerides information
in a different way. Note that the segments stored in a single SPK file can be
of different types. A detailed description of the organization of the arrays for
each SPK type can be found at https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html

SPK files can be read with \link{readBinSPK}. This function allows the evaluation of
one segment of an SPK file to any target times. The evaluated segment should be
one of the elements in the list of segments returned by calling \link{readBinSPK}
on an SPK file. 

Evaluation at target times of the SPK segment produces state vectors comprising
X, Y and Z components for position and velocity for the target body. The target
body, center body and frame of reference are those indicated in the corresponding
summary of the SPK segment. 

The specific algorithm through which an SPK segment is evaluated depends on the
type of SPK segment. For a detailed description, see the documentation of \link{readBinSPK}
or NAIF's documentation.

The output values for all types of segments have been verified to match those
of CSPICE to a precision of \code{sqrt(.Machine$double.eps)}.

Note that, in addition to position and velocity components, \code{evaluateSPK}
will also produce X, Y and Z components of acceleration for segments of types
1, 2, 3, 8, 9, 12, 13, 14, 18, 19, 20 and 21. This is possible because all of
these segments provide some sort of interpolation polynomial coefficients, which
can be differentiated to obtain acceleration values. However, note that said
coefficients are, in principle, not intended to be used to calculate acceleration
values, and SPICE does not return these. Therefore, proceed with caution if 
using acceleration values obtained in this way. In the specific case of SPK segments
of types 1 and 21, the provided coefficents are for an interpolation polynomial 
for acceleration, and therefore are probably the most reliable. In particular,
for target times exactly matching the reference epoch of any of the data points
included in the segment, the acceleration value should match the original value
used to fit the interpolation polynomial.
}
\usage{
evaluateSPKSegment(segment, targetEpochs)
}
\arguments{
  \item{segment}{Single SPK segment to be evaluated. The segment should have the
  same structure as that of segments returned by the \link{readBinSPK} function.
  Note said function returns all segments contained in the read SPK file; only
  one of those should be selected and provided here.
  }
  \item{targetEpochs}{Numeric vector indicating the target epochs at which the
  segment should be evaluated. The epochs should be provided in TDB seconds since
  J2000, also known as ephemeris time in SPICE. Note that all target epochs should
  be larger than the start epoch of the segment, and smaller than the end epoch of
  the segment, both of which are specified in the corresponding segment summary.
  Segments of type 1 and 21 are an exception. There are segments of these types
  where the final record covers epochs up to a time larger than the global end
  epoch of the corresponding segment; in cases where epochs larger than the global
  end epoch of the segment, but smaller than the end epoch of the last data record,
  these epochs will also be evaluated.
  }
}
\value{
A matrix with 7 (for SPK segments of types 5, 10 and 15) or 10 (for all other
types) columns. Each row in the matrix represents a target time of evaluation.
Column 1 provides the epochs of evaluation. Columns 2-4 provide position components.
Columns 5-7 provide velocity components. Column 8-10, when present, provide
acceleration components.
}
\references{
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
Shampine, L. F. and Gordon, M. K., Computer Solution of Ordinary Differential Equations: The Initial Value Problem, 1975
Robert Werner, SPICE spke01 math, 2022.
https://doi.org/10.5270/esa-tyidsbu
}
\examples{
# The file vgr2_jup230.bsp provided with the package includes information for the
# Jupiter flyby of Voyager 2

testSPK <- readBinSPK(paste0(path.package("asteRisk"), "/vgr2_jup230.bsp"))
length(testSPK$segments)

# It contains a single segment.

testSegment <- testSPK$segments[[1]]

# Check the initial and end epochs of the interval covered by the segment

testSegment$segmentSummary$initialEpoch
testSegment$segmentSummary$finalEpoch

# Evaluate at target epochs

evaluateSPKSegment(testSegment, c(-649364400, -649364374.68,  -647364600, -645364800))
}
