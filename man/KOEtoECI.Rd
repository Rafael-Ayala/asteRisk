\name{KOEtoECI}
\alias{KOEtoECI}
\title{Calculate ECI coordinates from Keplerian orbital elements
}
\description{
Keplerian orbital elements are a set of six parameters used to described the
orbits of celestial objects, including satellites. While satellites do not
follow a perfectly Keplerian orbit, their state at any point can be defined
by the orbital parameters that they would have if they were located at the same
position with the same velocity following a perfectly Keplerian orbit (i.e.,
if perturbations were absent). These are called osculating orbital elements.

A complete set of six Keplerian elements defines unequivocally the position
and velocity of the satellite in a given frame of reference, and therefore can
be used to calculate its cartesian coordinates. This function calculates the
coordinates of a satellite in an ECI (Earth-centered inertial) frame of
reference from a set of Keplerian orbital elements. The exact ECI frame of the
resulting coordinates is the same used to define the supplied orbital elements.
}
\usage{
KOEtoECI(a, e, i, M, omega, OMEGA, keplerAccuracy=10e-8, maxKeplerIterations=100)
}
\arguments{
  \item{a}{Semi-major axis of orbital ellipse in meters. 
  }
  \item{e}{Numerical eccentricity of the orbit. Eccentricity measures how much
  the orbit deviates from being circular.
  }
  \item{i}{Inclination of the orbital plane in radians. Inclination is the
  angle between the orbital plane and the equator.
  }
  \item{M}{Mean anomaly of the orbit in radians. Mean anomaly indicates where
  the satellite is along its orbital path, and is defined as the angle between
  the direction of the perigee and the hypothetical point where the object would 
  be if it was moving in a circular orbit with the same period as its true orbit 
  after the same amount of time since it last crossed the perigee had ellapsed.
  }
  \item{omega}{Argument of perigee in radians. This is the angle between the 
  direction of the ascending node and the direction of the perigee (the point of 
  the orbit at which the object is closest to the Earth).
  }
  \item{OMEGA}{Right ascension of the ascending node in radians. This is the 
  angle between the direction of the ascending node (the point where the 
  satellite crosses the equatorial plane moving north) and the direction of the 
  First Point of Aries (which indicates the location of the vernal equinox).
  }
  \item{keplerAccuracy}{Accuracy to consider Kepler's equation solved when
  calculating eccentric anomaly from mean anomaly. If two consecutive solutions 
  differ by a value lower than this accuracy, integration is considered to have 
  converged.
  }
  \item{maxKeplerIterations}{Maximum number of iterations after which 
  fixed-point integration of Kepler's equation will stop, even if convergence 
  according to the accuracy criterion has not been reached.
  }
}
\value{
A list with two elements representing the position and velocity of the satellite
in the same ECI (Earth Centered, Earth Fixed) frame of reference into which
the provided orbital elements were defined. Position values are in m, and 
velocity values are in m/s. Each of the two elements contains three values,
corresponding to the X, Y and Z components of position and velocity in this 
order.
}
\references{
https://www.gsc-europa.eu/system-service-status/orbital-and-technical-parameters
https://celestrak.org/columns/v02n01/
https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
}
\examples{
# Let's calculate the ECI coordinates from the orbital elements provided by a
# TLE. It should be noted that this is often not recommended, since the orbital
# elements supplied in a TLE are not osculating orbital elements, but instead
# mean orbital elements set to fit a range of actual observations. The
# recommended procedures are to use TLE only in conjunction with the SGP4/SDP4
# models, and viceversa.
# The following orbital parameters correspond to an object with NORAD catalogue
# number 24208 (Italsat 2) the 26th of June, 2006 at 00:58:29.34 UTC.

n0 <- 1.007781*((2*pi)/(86400))  # Multiplication by 2pi/86400 to convert to radians/s
e0 <- 0.002664 # mean eccentricity at epoch
i0 <- 3.8536*pi/180 # mean inclination at epoch in radians
M0 <- 48.3*pi/180 # mean anomaly at epoch in radians
omega0 <- 311.0977*pi/180 # mean argument of perigee at epoch in radians
OMEGA0 <- 80.0121*pi/180 # mean longitude of ascending node at epoch in radians

# The semi-major axis can be calculated from the mean motion in radians/s 
# as follows: (mu is the standard gravitational parameter of Earth)

mu <- 3.986004418e14 # in units of m3 s-2
a0 <- (mu^(1/3))/(n0^(2/3))

# The ECI coordinates can then be calculated. In this case, they will be in TEME
# frame, since the original orbital elements are derived from a TLE
coordinates_ECI <- KOEtoECI(a0, e0, i0, M0, omega0, OMEGA0)
}
