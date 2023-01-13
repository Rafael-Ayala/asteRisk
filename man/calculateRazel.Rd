\name{calculateRazel}
\alias{calculateRazel}
\title{Calculates azimuth, elevation and range of a given object
}
\description{
The horizontal coordinate system, also called azimuth-elevation system, uses the 
local horizon of an observer as its fundamental plane. In it, a given point is
defined by 2 main angles: azimuth and elevation. Azimuth defines the angle
of the point around the horizon in the X-Y plane, measured from the true North 
and usually increasing towards the East. Elevation is the angle between the 
object and the X-Y plane. Finally, the range defines the distance between the
observer and the point.

This function calculates the azimuth, elevation and range given the coordinates
of an observed satellite and of an observer. Both sets of coordinates must be
provided as Cartesian geocentric coordinates in ITRF.
}
\usage{
calculateRazel(geocentricObserver, geocentricSatellite, degreesOutput=TRUE)
}
\arguments{
  \item{geocentricObserver}{Vector with the X, Y and Z components of the position
  of the observer in ITRF frame.
  }
  \item{geocentricSatellite}{Vector with the X, Y and Z components of the position
  of the satellite in ITRF frame.
  }
  \item{degreesOutput}{Logical indicating if the output should be in sexagesimal
  degrees. If \code{degreesOutput=FALSE}, the output will be in radians.
  }
}
\value{
A vector with three elements, corresponding to the azimuth and elevation in
degrees (or radians if specified) and the range in the same unit as the provided
Cartesian coordinates.
}
\references{
https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
}
\examples{
# The following were the coordinates of Italsat-2 in ITRF the 27th of June, 2006
# at 00:58:29.34 UTC, in meters.

italsat_ITRF <- c(-37325542.8, 19152438.3, 138384.5)

# Let us calculate its azimuth, elevation and range for an observer from Tokyo.
# The latitude and longitude of the city are 35.6762 degrees North, 139.6503 
# degrees East. Let's assume an observer placed at sea level (0 m)
# We first convert these coordinates to ITRF:

observer_ITRF <- LATLONtoITRF(c(35.6762, 139.6503, 0), degreesInput=TRUE)

# We can now calculate the azimuth and elevation:

razel <- calculateRazel(observer_ITRF, italsat_ITRF, degreesOutput=TRUE)
razel[1] # Azimuth
razel[2] # Elevation
}
