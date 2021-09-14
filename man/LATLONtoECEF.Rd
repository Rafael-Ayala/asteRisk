\name{LATLONtoECEF}
\alias{LATLONtoECEF}
\title{Convert coordinates from geodetic latitude, longitude and altitude to ECEF
}
\description{
The ECEF (Earth Centered, Earth Fixed) frame of reference is a non-inertial
coordinate frame, where the origin is placed at the center of mass of the Earth
and the frame rotates with respect to the stars to remain fixed with respect to 
the Earth surface as it rotates. This function converts geodetic latitude, 
longitude and altitude to Cartesian coordinates in the ECEF frame. The WGS84
Earth ellipsoid model is used.
}
\usage{
LATLONtoECEF(position_LATLON, degreesInput=TRUE)
}
\arguments{
  \item{position_LATLON}{Vector with the latitude, longitude and altitude of the
  object. Latitude and longitude can be provided in sexagesimal degrees or in 
  radians (by default, decimal degrees are asumed). Altitude must be provided in 
  meters.
  }
  \item{degreesInput}{Logical indicating if the input latitude and longitude are
  in sexagesimal degrees. If \code{degreesInput=FALSE}, the input will be
  considered to be in radians.
  }
}
\value{
A vector with three elements, corresponding to the  X, Y and Z components of 
position in meters in the ECEF frame, in this order.
}
\references{
https://apps.dtic.mil/sti/pdfs/ADA280358.pdf
}
\examples{
latitude <- 37.3891
longitude <- -5.9845
altitude <- 20000

# Let´s calculate the corresponding coordinates in ECEF frame

coordinates_ECEF <- LATLONtoECEF(c(latitude, longitude, altitude))
}
