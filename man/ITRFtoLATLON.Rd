\name{ITRFtoLATLON}
\alias{ITRFtoLATLON}
\title{Convert coordinates from ITRF to geodetic latitude, longitude and altitude
}
\description{
The ITRF (International Terrestrial Reference Frame) is an ECEF (Earth Centered, 
Earth Fixed) frame of reference, i.e., a non-inertial frame of reference where 
the origin is placed at the center of mass of Earth, and the frame rotates 
with respect to the stars to remain fixed with respect to the Earth surface as 
it rotates. The Z-axis extends along the true North as defined by the IERS
reference pole, and the X-axis extends towards the intersection between the
equator and the Greenwich meridian at any time. This function converts Cartesian 
coordinates in the ECEF frame to geodetic latitude, longitude and altitude.
}
\usage{
ITRFtoLATLON(position_ITRF, degreesOutput=TRUE)
}
\arguments{
  \item{position_ITRF}{Vector with the X, Y and Z components of the position
  of an object in ITRF frame, in m.
  }
  \item{degreesOutput}{Logical indicating if the output should be in sexagesimal
  degrees. If \code{degreesOutput=FALSE}, the output will be in radians.
  }
}
\value{
A vector with three elements, corresponding to the latitude and longitude in
degrees (or radians if specified) and the altitude in m.
}
\references{
https://arc.aiaa.org/doi/10.2514/6.2006-6753
}
\examples{
coordinates_ITRF <- c(5062040.1, -530657.4, 3863936.5)

# Let's calculate the geodetic latitude, longitude and altitude

geodetic <- ITRFtoLATLON <- (coordinates_ITRF)

# A longer application example follows, useful to represent the groundtrack of a 
# satellite after propagaton with SGP4/SDP4

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

# We can now convert the results in TEME frame to ITRF frame, previously
# multiplying by 1000 to convert the km output of sgdp4 to m

state_1day_ITRF <- TEMEtoITRF(state_1day_TEME$position, state_1day_TEME$velocity,
                              "2006-06-27 00:58:29.34")
                              
# Finally, we can convert the ECEF coordinates to geodetic latitude, longitude
# and altitude

state_1day_geodetic <- ITRFtoLATLON(state_1day_ITRF$position)
}
}
