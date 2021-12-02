\name{LATLONtoGCRF}
\alias{LATLONtoGCRF}
\title{Convert coordinates from geodetic latitude, longitude and altitude to GCRF
}
\description{
The GCRF (Geocentric Celestial Reference Frame) frame of reference is an 
Earth-centered inertial coordinate frame, where the origin is also placed at the 
center of mass of Earth and the coordinate frame is fixed with respect to the 
stars (and therefore not fixed with respect to the Earth surface in its 
rotation). The X-axis is aligned with the mean equinox of Earth at 12:00
Terrestrial Time on the 1st of January, 2000, and the Z-axis is aligned with
the Earth´s rotation axis. This function converts geodetic latitude, 
longitude and altitude to Cartesian coordinates in the GCRF frame. The WGS84
Earth ellipsoid model is used.
}
\usage{
LATLONtoGCRF(position_LATLON, dateTime, degreesInput=TRUE)
}
\arguments{
  \item{position_LATLON}{Vector with the latitude, longitude and altitude of the
  object. Latitude and longitude can be provided in sexagesimal degrees or in 
  radians (by default, sexagesimal degrees are asumed). Altitude must be 
  provided in meters.
  }
  \item{dateTime}{Date-time string with the date and time in UTC corresponding
  to the provided geodetic coordinates.
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
if(requireNamespace("asteRiskData", quietly = TRUE)) {
latitude <- 37.3891
longitude <- -5.9845
altitude <- 20000

# Let´s calculate the corresponding coordinates in GCRF frame, assuming that
# the provided coordinates were valid the 20th of October, 2020 at 15:00:00 UTC

coordinatesGCRF <- LATLONtoGCRF(c(latitude, longitude, altitude), 
                                dateTime="2020-10-20 15:00:00")
}
}
