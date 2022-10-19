\name{TEMEtoGCRF}
\alias{TEMEtoGCRF}
\title{Convert coordinates from TEME to GCRF
}
\description{
The TEME (True Equator, Mean Equinox) and GCRF (Geocentric Celestial Reference
Frame) are both ECI frames of reference, i.e., Earth-centered inertial 
coordinate frames, where the origin is placed at the center of mass
of Earth and the coordinate frame is fixed with respect to the stars (and
therefore not fixed with respect to the Earth surface in its rotation). 

The difference between the two resides in the fact that in the GCRF frame,
the X-axis and Z-axis are aligned respectively with the mean equinox and 
rotation axis of Earth at 12:00 Terrestrial Time on the 1st of January, 2000,
while in the TEME frame they are aligned with the mean equinox and rotation
axis at the time of the corresponding TLE. Due to the change of the direction
of the vernal equinox and the rotation axis over time, coordinates in the two
frames differ slightly.

This function requires the \code{asteRiskData} package, which can be installed 
by running \code{install.packages('asteRiskData', repos='https://rafael-ayala.github.io/drat/')}
}
\usage{
TEMEtoGCRF(position_TEME, velocity_TEME, dateTime)
}
\arguments{
  \item{position_TEME}{Vector with the X, Y and Z components of the position
  of an object in TEME frame, in m.
  }
  \item{velocity_TEME}{Vector with the X, Y and Z components of the velocity
  of an object in TEME frame, in m/s.
  }
  \item{dateTime}{Date-time string with the date and time in UTC corresponding
  to the provided position and velocity vectors. This specifies the time for 
  which the conversion from TEME to GCRF coordinates will be performed. It is
  required due to the change in the exact position of the rotation axis of 
  Earth due to precesion, nutation and polar motion.
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
https://celestrak.org/columns/v04n03/#FAQ01
}
\examples{
if(requireNamespace("asteRiskData", quietly = TRUE)) {
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

# We can now convert the results in TEME frame to GCRF frame, previously
# multiplying by 1000 to convert the km output of sgdp4 to m

state_1day_GCRF <- TEMEtoGCRF(state_1day_TEME$position*1000, 
                              state_1day_TEME$velocity*1000,
                              "2006-06-27 00:58:29.34")
}
}
