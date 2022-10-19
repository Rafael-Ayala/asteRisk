\name{ITRFtoGCRF}
\alias{ITRFtoGCRF}
\title{Convert coordinates from ITRF to GCRF
}
\description{
The ITRF (International Terrestrial Reference Frame) is an ECEF (Earth Centered, 
Earth Fixed) frame of reference, i.e., a non-inertial frame of reference where 
the origin is placed at the center of mass of Earth, and the frame rotates 
with respect to the stars to remain fixed with respect to the Earth surface as 
it rotates. The Z-axis extends along the true North as defined by the IERS
reference pole, and the X-axis extends towards the intersection between the
equator and the Greenwich meridian at any time.

The GCRF (Geocentric Celestial Reference Frame) frame of reference is an 
Earth-centered inertial coordinate frame, where the origin is also placed at the 
center of mass of Earth and the coordinate frame is fixed with respect to the 
stars (and therefore not fixed with respect to the Earth surface in its 
rotation). The X-axis is aligned with the mean equinox of Earth at 12:00
Terrestrial Time on the 1st of January, 2000, and the Z-axis is aligned with
the Earth´s rotation axis.

This function requires the \code{asteRiskData} package, which can be installed 
by running \code{install.packages('asteRiskData', repos='https://rafael-ayala.github.io/drat/')}
}
\usage{
ITRFtoGCRF(position_ITRF, velocity_ITRF, dateTime)
}
\arguments{
  \item{position_ITRF}{Vector with the X, Y and Z components of the position
  of an object in ITRF frame, in m.
  }
  \item{velocity_ITRF}{Vector with the X, Y and Z components of the velocity
  of an object in ITRF frame, in m/s.
  }
  \item{dateTime}{Date-time string with the date and time in UTC corresponding
  to the provided position and velocity vectors. This specifies the time for 
  which the conversion from ITRF to GCRF coordinates will be performed. It is
  important to provide an accurate value, since the point over the surface of
  Earth to which a set of GCRF coordinates refers varies with time due to
  the motion of Earth.
  }
}
\value{
A list with two elements representing the position and velocity of the satellite
in the GCRF (Earth-centered non-intertial) frame of reference. Position values
are in m, and velocity values are in m/s. Each of the two elements contains
three values, corresponding to the X, Y and Z components of position and
velocity in this order.
}
\references{
https://celestrak.org/columns/v02n01/
}
\examples{
if(requireNamespace("asteRiskData", quietly = TRUE)) {
# The following were the position and velocity of satellite MOLNIYA 1-83 
# the 25th of June, 2006 at 00:33:43 UTC in the ECEF frame (in m and m/s).

position_ITRF <- c(1.734019e+06, -1.510972e+07, 39.08228)
velocity_ITRF <- c(1468.832, -3962.464, 4007.039)

# Let´s convert them into the GCRF frame

coordinates_GCRF <- ITRFtoGCRF(position_ITRF, velocity_ITRF, "2006-06-25 00:33:43")
}
}
