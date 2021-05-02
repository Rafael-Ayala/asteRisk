\name{ECItoKOE}
\alias{ECItoKOE}
\title{Calculate ECI coordinates from Keplerian orbital elements
}
\description{
Keplerian orbital elements are a set of six parameters used to described the
orbits of celestial objects, including satellites. While satellites do not
follow a perfectly Keplerian orbit, their state at any point can be defined
by the orbital parameters that they would have if they were located at the same
position with the same velocity following a perfectly Keplerian orbit (i.e.,
if perturbations were absent). These are called osculating orbital elements.

Keplerian orbital elements can be unequivocally determined from a satellite
if its position and velocity are known. This function calculates orbital
elements from the position and velocity of a satellite in an ECI (Earth-centered
inertial) frame of reference. The elements (such as the equatorial plane) with 
respect to which the resulting orbital elements will be defined are the same as 
those used for the ECI frame of reference. The function calculates the six
standard orbital elements, plus some alternative elements useful for the 
characterization of special orbits, such as circular ones or orbits with no
inclination.
}
\usage{
ECItoKOE(position_ECI, velocity_ECI)
}
\arguments{
  \item{position_ECI}{Vector with the X, Y and Z components of the position
  of an object in an ECI frame, in m.
  }
  \item{velocity_ECI}{Vector with the X, Y and Z components of the velocity
  of an object in an ECI frame, in m/s.
  }
}
\value{
A list with the following standard and alternative orbital elements:
  \item{semiMajorAxis}{Semi-major axis of orbital ellipse in meters.
  }
  \item{eccentricity}{Numerical eccentricity of the orbit. Eccentricity measures how much
  the orbit deviates from being circular.
  }
  \item{inclination}{Inclination of the orbital plane in radians. Inclination is the
  angle between the orbital plane and the equator.
  }
  \item{meanAnomaly}{Mean anomaly of the orbit in radians. Mean anomaly indicates where
  the satellite is along its orbital path, and is defined as the angle between
  the direction of the perigee and the hypothetical point where the object would 
  be if it was moving in a circular orbit with the same period as its true orbit 
  after the same amount of time since it last crossed the perigee had ellapsed.
  }
  \item{argumentPerigee}{Argument of perigee in radians. This is the angle between the 
  direction of the ascending node and the direction of the perigee (the point of 
  the orbit at which the object is closest to the Earth).
  }
  \item{longitudeAscendingNode}{Longitude of the ascending node (also called
  right ascension of the ascending node) in radians. This is the angle between 
  the direction of the ascending node (the point where thesatellite crosses the 
  equatorial plane moving north) and the direction of the First Point of Aries 
  (which indicates the location of the vernal equinox).
  }
  \item{trueAnomaly}{True anomaly of the orbit in radians. Unlike mean anomaly,
  true anomaly is the angle between the direction of the perigee and the actual
  position of the satellite.
  }
  \item{argumentLatitude}{Argument of latitude of the orbit in radians. Defined
  as the angle between the equator and the position of the satellite. It is
  useful to define the position of satellites in circular orbits, where the
  argument of perigee and true anomaly are not well defined.
  }
  \item{longitudePerigee}{Longitude of perigee of the orbit in radians. Defined
  as the angle between the vernal equinox and the perigee. It is useful for
  cases of orbits with 0 inclination, where the longitude of the ascending node
  and the argument of perigee are not well defined.
  }
  \item{trueLongitude}{Longitude of perigee of the orbit in radians. Defined
  as the angle between the vernal equinox and the position of the satellite. 
  It is useful for cases of circular orbits with 0 inclination, where the 
  longitude of the ascending node, the argument of perigee and true anomaly are
  not well defined.
  }
}
\references{
https://www.gsc-europa.eu/system-service-status/orbital-and-technical-parameters
https://celestrak.com/columns/v02n01/
https://www.faa.gov/about/office_org/headquarters_offices/avs/offices/aam/cami/library/online_libraries/aerospace_medicine/tutorial/media/iii.4.1.4_describing_orbits.pdf
}
\examples{
# The following were the position and velocity of satellite MOLNIYA 1-83 
# the 25th of June, 2006 at 00:33:43 UTC in the GCRF frame (in m and m/s).

position_GCRF <- c(-14471729.582, -4677558.558, 9369.461)
velocity_GCRF <- c(-3251.691, -3276.008, 4009.228)

# Let's calculate the orbital elements of the satellite at that time

orbital_elements <- ECItoKOE(position_GCRF, velocity_GCRF)
}
