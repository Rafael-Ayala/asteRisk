\name{GCRFtoLATLON}
\alias{GCRFtoLATLON}
\title{Convert coordinates from GCRF to geodetic latitude, longitude and altitude
}
\description{
The GCRF (Geocentric Celestial Reference Frame) frame of reference is an 
Earth-centered inertial coordinate frame, where the origin is placed at the 
center of mass of Earth and the coordinate frame is fixed with respect to the 
stars (and therefore not fixed with respect to the Earth surface in its 
rotation). The X-axis is aligned with the mean equinox of Earth at 12:00
Terrestrial Time on the 1st of January, 2000, and the Z-axis is aligned with
the Earth´s rotation axis. This function converts position in GCRF to geodetic 
latitude, longitude and altitude, which can be considered to be a non-inertial, 
Earth-centered frame of reference.

This function requires the \code{asteRiskData} package, which can be installed 
by running \code{install.packages('asteRiskData', repos='https://rafael-ayala.github.io/drat/')}
}
\usage{
GCRFtoLATLON(position_GCRF, dateTime, degreesOutput=TRUE)
}
\arguments{
  \item{position_GCRF}{Vector with the X, Y and Z components of the position
  of an object in TEME frame, in m.
  }
  \item{dateTime}{Date-time string with the date and time in UTC corresponding
  to the provided position vector.
  }
  \item{degreesOutput}{Logical indicating if the output should be in sexagesimal
  degrees. If \code{degreesOutput=FALSE}, the output will be in radians.
  }
}
\value{
A vector with three elements, corresponding to the latitude and longitude in
degrees and the altitude in m.
}
\references{
https://arc.aiaa.org/doi/10.2514/6.2006-6753
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

# Finally, we convert the results in GCRF frame to geodetic latitude, longitude
# and altitude

state_1day_geodetic <- GCRFtoLATLON(state_1day_GCRF$position, "2006-06-27 00:58:29.34")
}
}
