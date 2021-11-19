\name{LATLONtoITRF}
\alias{LATLONtoITRF}
\title{Convert coordinates from geodetic latitude, longitude and altitude to ITRF
}
\description{
The ITRF (International Terrestrial Reference Frame) is an ECEF (Earth Centered, 
Earth Fixed) frame of reference, i.e., a non-inertial frame of reference where 
the origin is placed at the center of mass of Earth, and the frame rotates 
with respect to the stars to remain fixed with respect to the Earth surface as 
it rotates. The Z-axis extends along the true North as defined by the IERS
reference pole, and the X-axis extends towards the intersection between the
equator and the Greenwich meridian at any time. This function converts geodetic 
latitude, longitude and altitude to Cartesian coordinates in the ITRF frame. The 
WGS84 Earth ellipsoid model is used.
}
\usage{
LATLONtoITRF(position_LATLON, degreesInput=TRUE)
}
\arguments{
  \item{position_LATLON}{Vector with the latitude, longitude and altitude of the
  object. Latitude and longitude can be provided in sexagesimal degrees or in 
  radians (by default, sexagesimal degrees are asumed). Altitude must be 
  provided in meters.
  }
  \item{degreesInput}{Logical indicating if the input latitude and longitude are
  in sexagesimal degrees. If \code{degreesInput=FALSE}, the input will be
  considered to be in radians.
  }
}
\value{
A vector with three elements, corresponding to the  X, Y and Z components of 
position in meters in the ITRF frame, in this order.
}
\references{
https://apps.dtic.mil/sti/pdfs/ADA280358.pdf
}
\examples{
latitude <- 37.3891
longitude <- -5.9845
altitude <- 20000

# Let´s calculate the corresponding coordinates in ECEF frame

coordinates_ITRF <- LATLONtoITRF(c(latitude, longitude, altitude))
}
