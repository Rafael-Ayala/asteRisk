\name{ECEFtoGCRF}
\alias{ECEFtoGCRF}
\title{Convert coordinates from ECEF to GCRF
}
\description{
The ECEF (Earth Centered, Earth Fixed) is a non-inertial frame of reference
where the origin is placed at the center of mass of Earth, and the frame rotates 
with respect to the stars to remain fixed with respect to the Earth surface as 
it rotates.

The GCRF (Geocentric Celestial Reference Frame) frame of reference is an 
Earth-centered inertial coordinate frame, where the origin is also placed at the 
center of mass of Earth and the coordinate frame is fixed with respect to the 
stars (and therefore not fixed with respect to the Earth surface in its 
rotation). The X-axis is aligned with the mean equinox of Earth at 12:00
Terrestrial Time on the 1st of January, 2000, and the Z-axis is aligned with
the Earth´s rotation axis.
}
\usage{
ECEFtoGCRF(position_ECEF, velocity_ECEF, dateTime)
}
\arguments{
  \item{position_ECEF}{Vector with the X, Y and Z components of the position
  of an object in ECEF frame, in m.
  }
  \item{velocity_ECEF}{Vector with the X, Y and Z components of the velocity
  of an object in ECEF frame, in m/s.
  }
  \item{dateTime}{Date-time string with the date and time in UTC corresponding
  to the provided position and velocity vectors.
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
https://celestrak.com/columns/v02n01/
}
\examples{
# The following were the position and velocity of satellite MOLNIYA 1-83 
# the 25th of June, 2006 at 00:33:43 UTC in the ECEF frame (in m and m/s).

position_ECEF <- c(1.734019e+06, -1.510972e+07, 39.08228)
velocity_ECEF <- c(1468.832, -3962.464, 4007.039)

# Let´s convert them into the GCRF frame

coordinates_GCRF <- ECEFtoGCRF(position_ECEF, velocity_ECEF, "2006-06-25 00:33:43")
}
