\name{TEMEtoLATLON}
\alias{TEMEtoLATLON}
\title{Convert coordinates from TEME to geodetic latitude, longitude and altitude
}
\description{
The TEME (True Equator, Mean Equinox) frame of reference is an Earth-centered
inertial coordinate frame, where the origin is placed at the center of mass
of Earth and the coordinate frame is fixed with respect to the stars (and
therefore not fixed with respect to the Earth surface in its rotation). The
coordinates and velocities calculated with the SGP4 and SDP4 models are in the
TEME frame of reference. This function converts positions and velocities in
TEME to geodetic latitude, longitude and altitude, which can be considered to be
a non-inertial, Earth-centered frame of reference.
}
\usage{
TEMEtoLATLON(position_TEME, dateTime)
}
\arguments{
  \item{position_TEME}{Vector with the X, Y and Z components of the position
  of an object in TEME frame, in km.
  }
  \item{dateTime}{Date-time string with the date and time in UTC corresponding
  to the provided position and velocity vectors.
  }
}
\value{
A vector with three elements, corresponding to the latitude and longitude in
degrees and the altitude in km.
}
\references{
https://arc.aiaa.org/doi/10.2514/6.2006-6753
}
\examples{
# The following orbital parameters correspond to an object with NORAD catalogue
# number 24208 (Italsat 2) the 26th of June, 2006 at 00:58:29.34 UTC.

n0 <- 1.007781*((2*pi)/(1440))  # Multiplication by 2pi/1440 to convert to radians/min
e0 <- 0.002664 # mean eccentricity at epoch
i0 <- 3.8536*pi/180 # mean inclination at epoch in radians
M0 <- 48.3*pi/180 # mean anomaly at epoch in radians
omega0 <- 311.0977*pi/180 # mean argument of perigee at epoch in radians
OMEGA0 <- 80.0121*pi/180 # mean longitude of ascending node at epoch in radians
Bstar <- 1e-04 # drag coefficient
epochDateTime <- "2006-06-26 00:58:29.34"

# Let´s calculate the position and velocity of the satellite 1 day later

state_1day_TEME <- sgdp4(n0=n0, e0=e0, i0=i0, M0=M0, omega0=omega0, OMEGA0=OMEGA0,
                         Bstar=Bstar, initialDateTime=epochDateTime, targetTime=1440)

# We can now convert the results in TEME frame to ECEF frame

state_1day_ECEF <- TEMEtoECEF(state_1day_TEME$position, state_1day_TEME$velocity,
                              "2006-06-27 00:58:29.34")
                              
# Finally, we can convert the ECEF coordinates to geodetic latitude, longitude
# and altitude

state_1day_geodetic <- ECEFtoLATLON(state_1day_ECEF$position)
}
