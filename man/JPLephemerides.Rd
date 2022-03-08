\name{JPLephemerides}
\alias{JPLephemerides}
\title{Calculate JPL main celestial objects ephemerides for a given Modified Julian Date
}
\description{
NASA's Jet Propulsion Laboratory (JPL) provides mathematical models of the Solar
System known as Development Ephemerides (DE). The models are given as sets of 
Chebyshev coefficients, which cam be used to calculate the position (and its 
derivatives) of the Sun, the eight major planets, Pluto and the Moon.
This function employes JPL DE440 to calculate the position (and optionally
velocities also) of the mentioned celestial objects, in ICRF frame. 
JPL DE440 covers the period from 1550 to 2650 AC. In addition to the position of 
celestial objects, lunar libration angles are also calculated. Internally, 
calculations are performed by employing Clenshaw's algorithm together with the 
Chebyshev coefficients provided by JPL DE440.
The target time should be specified as a Modified Julian Date (MJD). MJD in different
time systems can be used. Currently, UTC, UT1, TT and TDB are supported.
Additionally, a central body with respect to which positions and velocities are
calculated should be specified. By default, the Solar System Barycenter (SSB) is
used, but additionally Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus,
Neptune or Pluto can be selected.
Note that this function requires the additional package asteRiskData, which
provides the Chebyshev coefficients, and can be installed by running
\code{install.packages("asteRiskData", repos="https://rafael-ayala.github.io/drat/")}
}
\usage{
JPLephemerides(MJD, timeSystem="UTC", centralBody="SSB", includeVelocities=TRUE)
}
\arguments{
  \item{MJD}{Modified Julian Date of the time for which celestial object ephemerides
  should be calculated. MJD are fractional number of days since midnight of the
  17th of November, 1858. The MJD of a date-time string can be obtained with
  function \link{dateTimeToMJD}.
  }
  \item{timeSystem}{Time system into which the MJD is provided. Should be one 
  from "UTC" (Coordinated Universal Time; default), "UT1" (Universal Time), 
  "TT" (Terrestrial Time) and "TDB" (Barycentric Dynamical Time).
  }
  \item{centralBody}{String indicating the celestial object that will be taken
  as the center of coordinates to which positions and velocities are referred.
  Must be one of "SSB" (Solar System Barycenter), "Mercury", "Venus", "Earth",
  "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune" or "Pluto".
  }
  \item{includeVelocities}{Logical indicating if velocities of the celestial
  objects should be calculated in addition to their positions. By default,
  they are calculated. Additionally, if set to TRUE, the first order derivatives
  of the libration angles of the Moon are also calculated.
  }
}
\value{
A list of vectors providing the positions (in meters), velocities (in m/s; only
if requested), Moon libration angles (in radians) and first derivatives of Moon
libration angles (in radians/s; only if velocities were requested) of celestial
objects with respect to the specified central body. For position and velocity
vectors, X, Y and Z components are given in this order. For Moon libration angles
and their derivatives, they are given in the following order: phi, theta and psi.
}
\references{
https://gssc.esa.int/navipedia/index.php/Julian_Date
https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems
}
\examples{
if(requireNamespace("asteRiskData", quietly = TRUE)) {
# Let's calculate the MJD of the 12th of June, 2000 at 10:00:00 UTC time, in UTC

MJD_UTC <- dateTimeToMJD("2000-06-12 10:00:00", timeSystem = "UTC")

# Let's now calculate the JPL ephemerides using Earth as the central body:

ephemerides <- JPLephemerides(MJD_UTC, timeSystem = "UTC", centralBody="Earth")

# We can now calculate, for example, the exact distance between the barycenters
# of Earth and Moon

sqrt(sum(ephemerides$positionMoon^2))
}
}
