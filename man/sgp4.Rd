\name{sgp4}
\alias{sgp4}
\title{Propagate an orbital state vector with the SGP4 model
}
\description{
Given an orbital state vector of a satellite, applies the SGP4 model to
propagate its orbit to the desired time point. This allows the calculation of
the position and velocity of the satellite at different times, both before
and after the time corresponding to the known state vector (referred to as
"epoch"). Kepler's equation is solved through fixed-point integration.
The SGP4 model can only accurately propagate the orbit of objects near Earth
(with an orbital period shorter than 225 minutes, corresponding approximately to
an altitude lower than 5877.5 km). For propagation of objects in deep space, the
SDP4 model should be used, available through the \code{\link{sdp4}} function.
This implementation is based on the theory and implementation described in
Space Track Report #3, and includes the corrections summarized in Revisiting
Space Track Report #3.
}
\usage{
sgp4(n0, e0, i0, M0, omega0, OMEGA0, Bstar, initialDateTime=NULL, targetTime,
     keplerAccuracy=10e-12, maxKeplerIterations=10)
}
\arguments{
  \item{n0}{Mean motion of the satellite at epoch in radians/min.
  }
  \item{e0}{Mean eccentricity of the orbit of the satellite at epoch.
  Eccentricity ranges from 0 (perfectly circular orbit) to 1 (parabolic
  trajectory).
  }
  \item{i0}{Mean orbital inclination of the satellite at epoch in radians.
  }
  \item{M0}{Mean anomaly of the satellite at epoch.
  }
  \item{omega0}{Mean argument of perigee of the satellite at epoch.
  }
  \item{OMEGA0}{Mean longitude of the ascending node of the satellite at epoch.
  Also known as right ascension of the ascending node.
  }
  \item{Bstar}{Drag coefficient of the satellite in units of (earth radii)^-1^.
  Bstar is an adjusted value of the ballistic coefficient of the satellite, and
  it indicates how susceptible it is to atmospheric drag.
  }
  \item{initialDateTime}{Optional date-time string in UTC indicating the time
  corresponding to the known state vector of the satellite. It must be provided
  if \code{targetTime} is provided as a date-time string.
  }
  \item{targetTime}{Time at which the position and velocity of the satellite
  should be calculated. It can be provided in two different ways: either as
  a number corresponding to the time in minutes counting from epoch at which
  the orbit should be propagated, or as a date-time string in UTC, in which case
  the date-time string for epoch must be provided through \code{initialDateTime}.
  }
  \item{keplerAccuracy}{Accuracy to consider Kepler´s equation solved. If two
  consecutive solutions differ by a value lower than this accuracy, integration
  is considered to have converged.
  }
  \item{maxKeplerIterations}{Maximum number of iterations after which fixed-point
  integration of Kepler's equation will stop, even if convergence according to
  the accuracy criterion has not been reached.
  }
}
\value{
A list with three elements. The first two elements represent the position and
velocity of the satellite at the target time, in the TEME (True Equator, Mean
Equinox) frame of reference. Position values are in km, and velocity values are
in km/s. Each of these two elements contains three values, corresponding to the
X, Y and Z components of position and velocity in this order. The third element
indicates the algorithm used to propagate the orbit (sgp4).
}
\references{
https://celestrak.org/NORAD/documentation/spacetrk.pdf
http://www.celestrak.org/publications/aiaa/2006-6753/AIAA-2006-6753.pdf
}
\examples{
# The following orbital parameters correspond to an object with NORAD catalogue
# number 88888 the 1st of October, 1980 at 23:41:24 UTC.

n0 <- 16.05824518*((2*pi)/(1440))  # Multiplication by 2pi/1440 to convert to radians/min
e0 <- 0.0086731 # mean eccentricity at epoch
i0 <- 72.8435*pi/180 # mean inclination at epoch in radians
M0 <- 110.5714*pi/180 # mean anomaly at epoch in radians
omega0 <- 52.6988*pi/180 # mean argument of perigee at epoch in radians
OMEGA0 <- 115.9689*pi/180 # mean longitude of ascending node at epoch in radians
Bstar <- 0.66816e-4 # drag coefficient

# Calculation of the orbital period

2*pi/n0

# The period is lower than 225 min, and therefore the SGP4 model is valid.
# Let´s calculate the position and velocity of the satellite 40 minutes after
# epoch

new_state <- sgp4(n0=n0, e0=e0, i0=i0, M0=M0, omega0=omega0, OMEGA0=OMEGA0,
                  Bstar=Bstar,targetTime = 40)
new_state
}
