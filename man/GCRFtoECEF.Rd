\name{GCRFtoECEF}
\alias{GCRFtoECEF}
\title{Convert coordinates from GCRF to ECEF
}
\description{
The GCRF (Geocentric Celestial Reference Frame) frame of reference is an 
Earth-centered inertial coordinate frame, where the origin is placed at the 
center of mass of Earth and the coordinate frame is fixed with respect to the 
stars (and therefore not fixed with respect to the Earth surface in its 
rotation). The X-axis is aligned with the mean equinox of Earth at 12:00
Terrestrial Time on the 1st of January, 2000, and the Z-axis is aligned with
the Earth´s rotation axis.

It is equivalent to the J2000 frame of reference (also called EME2000), and in 
some contexts it is also referred to as ICRF frame (although in its strict 
definition, the origin of coordinates is placed at the barycenter of the Solar 
System).

In the ECEF frame, the origin is also placed at the center of mass of Earth, 
but the frame rotates with respect to the stars to remain fixed with respect to 
the Earth surface as it rotates.

The coordinates and velocities input and calculated with the high-precision
orbital propagator (\link{hpop}) are in the GCRF frame of reference.

This function requires the \code{asteRiskData} package, which can be installed 
by running \code{install.packages('asteRiskData', repos='https://rafael-ayala.github.io/drat/')}
}
\usage{
GCRFtoECEF(position_GCRF, velocity_GCRF, dateTime)
}
\arguments{
  \item{position_GCRF}{Vector with the X, Y and Z components of the position
  of an object in GCRF frame, in m.
  }
  \item{velocity_GCRF}{Vector with the X, Y and Z components of the velocity
  of an object in GCRF frame, in m/s.
  }
  \item{dateTime}{Date-time string with the date and time in UTC corresponding
  to the provided position and velocity vectors.
  }
}
\value{
A list with two elements representing the position and velocity of the satellite
in the ECEF (Earth Centered, Earth Fixed) frame of reference. Position values
are in m, and velocity values are in m/s. Each of the two elements contains
three values, corresponding to the X, Y and Z components of position and
velocity in this order.
}
\references{
https://celestrak.com/columns/v02n01/
}
\examples{
if(requireNamespace("asteRiskData", quietly = TRUE)) {
# The following were the position and velocity of satellite MOLNIYA 1-83 
# the 25th of June, 2006 at 00:33:43 UTC in the GCRF frame (in m and m/s).

position_GCRF <- c(-14471729.582, -4677558.558, 9369.461)
velocity_GCRF <- c(-3251.691, -3276.008, 4009.228)

# Let´s convert them into the ECEF frame

coordinates_ECEF <- GCRFtoECEF(position_GCRF, velocity_GCRF, "2006-06-27 00:58:29.34")
}
}
