\name{hpop}
\alias{hpop}
\title{High-precision numerical orbital propagator
}
\description{
Given the position and velocity of a satellite at a given time (in the ICRF
system of coordinates centered on the Solar System Barycenter, any of the main planets,
Earth's Moon or Pluto), propagates its position by calculating its acceleration 
(based on a force model) and solving the resulting second-order ODE through 
numerical integration. This allows propagation of orbits with considerably 
higher accuracy than other propagators such as SGP4 and SDP4, but at the expense 
of a much higher computational cost. The forces and effects currently considered 
are gravitational attraction by the  Earth (using the GGM05C gravity model, with
spherical harmonics up to degree and order of 360); effects of Earth ocean and 
solid tides; gravitational attraction by the Moon (using the GRGM1200B gravity
model with spherical harmonics up to degree and order of 1200), Sun and planets 
(considered as point masses); solar radiation pressure; atmospheric drag, and relativistic 
effects. The forcefield is based on the forces described in Satellite Orbits: 
Models, Methods and Applications (Oliver Montenbruck and Eberhard Gill) and 
Fundamentals of Astrodynamics and Applications (David Vallado).
The NRLMSISE-00 model is used to calculate atmospheric density for the
calculation of atmospheric drag. The FES2014 model is used to calculate Earth 
geopotential model corrections due to ocean tides.
As mentioned before, the central body for the frame of reference can be any of the
Solar System Barycenter (SSB), any of the main planets, Earth's Moon or Pluto.
By default, it is assumed to be Earth, corresponding to GCRF (Geocentric ICRF).
The initial position will be checked against the position of said celestial bodies,
to identify if it falls under the Laplacian gravitational sphere of influence of any of them.
If this is the case, and it differs from the specified central body, the coordinate
system will be changed to be centered on the celestial body whose sphere of influence
includes the object of interest. This avoids instability in propagation.
The high-precision numerical orbital propagator requires the \code{asteRiskData}
package, which provides the data and coefficients required for calculation of
the modeled forces. \code{asteRiskData} can be installed by running 
\code{install.packages('asteRiskData', repos='https://rafael-ayala.github.io/drat/')}
}
\usage{
hpop(position, velocity, dateTime, times, satelliteMass, dragArea, 
     radiationArea, dragCoefficient, radiationCoefficient, 
     earthSphericalHarmonicsDegree=130, solidEarthTides=TRUE,
     oceanTides=TRUE, moonSphericalHarmonicsDegree=150, centralBody="Earth", 
     autoCentralBodyChange=TRUE, ...)
}
\arguments{
  \item{position}{Initial position of the satellite in the GCRF system of
  coordinates. Should be provided as a numeric vector with 3 components that
  indicate the X, Y and Z components of the position in meters.
  }
  \item{velocity}{Initial velocity of the satellite in the GCRF system of
  coordinates. Should be provided as a numeric vector with 3 components that
  indicate the X, Y and Z components of the position in meters/second.
  }
  \item{dateTime}{Date time string in the YYYY-MM-DD HH:MM:SS format indicating
  the time corresponding to the initial position and velocity, in UTC time.
  }
  \item{times}{Vector with the times at which the position and velocity of the
  satellite should be calculated, in seconds since the initial time.
  }
  \item{satelliteMass}{Mass of the satellite in kilograms.
  }
  \item{dragArea}{Effective area of the satellite for atmospheric drag in
  squared meters. If the way that a satellite will orient with respect to its
  velocity is not known, a mean cross-sectional area should be calculated
  assuming that the orientation of the satellite with respect to its velocity
  will vary uniformly. A decent estimate can be obtained with a flat-plate model,
  where the satellite is considered to be parallelepiped-shaped. The mean
  effective area can then be calculated as CSA = (S1 + S2 + S3 (+S4))/2, where
  S1, S2 and S3 are the areas of the three perpendicular surfaces of the model
  and S4 is an optional term to account for the area of solar panels (potential
  masking between the solar panels and the main surfaces is not considered; this
  might be partially accounted for by introducing a factor to reduce the 
  calculated effective area). 
  }
  \item{radiationArea}{Effective area of the satellite subject to the effect of
  radiation pressure in squared meters.
  }
  \item{dragCoefficient}{Drag coefficient (Cd) used for the calculation of
  atmospheric drag. For low Earth-orbiting satellites, a value of 2.2 is
  frequently employed if a better approximation is not available.
  }
  \item{radiationCoefficient}{Coefficient for the force resulting from radiation
  pressure. This parameter is usually referred to as reflectivity coefficient
  (Cr) and the value varies for different satellites and orbits. If unknown,
  a value of 1.2 is usually a decent approximation.
  }
  \item{earthSphericalHarmonicsDegree}{Maximum degree and order that should be 
  considered when calculating the Earth geopotential model. The model will be 
  complete up to the specified degree/order, i.e., all zonal, sectorial and
  tesseral spherical harmonics will be calculated. The maximum possible value
  is 360, since that is the highest degree and order of the Stokes' coefficients
  provided in the GGM05C model. Note that spherical harmonics for Earth gravity
  field will only be used if Earth is the central body for propagation; otherwise,
  only a point-mass attraction will be calculated.
  }
  \item{solidEarthTides}{Logical indicating if corrections of the Cnm and Snm
  Stokes' coefficients for the geopotential model due to solid Earth tides 
  should be performed, following IERS 2010 procedures and considering anelasticity
  of the Earth.
  }
  \item{oceanTides}{Logical indicating if corrections of the Cnm and Snm
  Stokes' coefficients for the geopotential model due to ocean tides should be 
  performed, using the FES2014 oceanic tides model.
  }
  \item{moonSphericalHarmonicsDegree}{Maximum degree and order that should be 
  considered when calculating the Moon gravity model. The model will be 
  complete up to the specified degree/order, i.e., all zonal, sectorial and
  tesseral spherical harmonics will be calculated. The maximum possible value
  is 1200, since that is the highest degree and order of the Stokes' coefficients
  provided in the GRGM1200B model. Note that spherical harmonics for Moon gravity
  field will only be used if Moon is the central body for propagation; otherwise,
  only a point-mass attraction will be calculated.
  }
  \item{centralBody}{Character string indicating the celestial body on which the
  supplied initial position (in ICRF) are centered. Should be one of "SSB" (meaning
  Solar System Barycenter), "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter",
  "Saturn", "Uranus", "Neptune" or "Pluto". The initial position will be checked 
  against the position of said celestial bodies, to identify if it falls under the 
  Laplacian gravitational sphere of influence of any of them. If this is the case, 
  and it differs from the specified central body, the coordinate system will be 
  changed to be centered on the celestial body whose sphere of influence includes 
  the object of interest.
  }
  \item{autoCentralBodyChange}{Logical indicating if the celestial object used as
  the center of coordinates should be automatically updated during propagation
  based on the radii of the spheres of influence of the main planets, the Moon
  and Pluto. By default, \code{autoCentralBodyChange=TRUE}.
  }
  \item{...}{Additional parameters to be passed to \link{ode} to control how
  numerical integration is performed. By default, the RADAU5 solver is used.
  }
}
\value{
A data frame with the results of the numerical integration at the requested times.
Each row contains the results for one of the requested times. The data frame 
contains 11 columns: time (indicating the time for the corresponding row
in seconds since the initial time), positionX, positionY, positionZ (indicating 
the X, Y and Z components of the position for that time in meters), velocityX, 
velocitydY and velocitydZ (indicating the X, Y and Z components of the velocity 
for that time in meters/second), accelerationX, accelerationY, accelerationZ 
(indicating the X, Y and Z components of the acceleration for that time in 
meters/second^2) and centralBody, indicating the central body of the frame of 
reference for the results for the corresponding time.
Positions and velocities are returned in the ICRF frame of reference, centered
in the celestial body specified in column centralBody. If \code{autoCentralBodyChange=TRUE},
the celestial body whose sphere of influence includes the object of interest
will be automatically used as the central body. Additionally, if transitions in
or out of the spheres of influence of the main celestial bodies are detected
during propagation of the trajectory, the central body will be automatically 
modified accordingly. If \code{autoCentralBodyChange=FALSE}, such automatic
changes of the central body will not be performed, and instead the user-specified
central body will be used at all times. Note, however, that it is not recommended
to perform propagation in a frame center at an object different than the celestial
body whose sphere of influence includes the target of propagation, since this
can lead to a substantial loss of accuracy. For details, see \href{https://etd.auburn.edu/bitstream/handle/10415/1133/Vautier_Mana_34.pdf}{M. Vautier, 2008}. 
Note that, if none of the spheres of influence of the 
planets, Moon or Pluto  included the object of interest, the center of the ICRF 
frame will be placed at the Solar System Barycenter.
}
\references{
Satellite Orbits: Models, Methods and Applications. Oliver Montenbruck and Eberhard Gill.
Fundamentals of Astrodynamics and Applications. David Vallado.
https://www.mathworks.com/matlabcentral/fileexchange/55167-high-precision-orbit-propagator
https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1144&context=smallsat
https://iopscience.iop.org/article/10.1088/1742-6596/911/1/012009/pdf
https://www.sciencedirect.com/science/article/pii/S1110016821000016
https://etd.auburn.edu/bitstream/handle/10415/1133/Vautier_Mana_34.pdf?sequence=1
}
\examples{
if(requireNamespace("asteRiskData", quietly = TRUE)) {
# The following are the position and velocity in the GCRF frame of satellite
# MOLNIYA 1-83 the 25th of June, 2006 at 00:33:43 UTC.

initialPosition <-c(-14568679.5026116, -4366250.78287623, 9417.9289105405)
initialVelocity <- c(-3321.17428902497, -3205.49400830455, 4009.26862308806) 
initialTime <- "2006-06-25 00:33:43"

# Molniya satellites have a mass of approximately 1600 kg and a cross-section of
# 15 m2. Additionally, let´s use 2.2 and 1.2 as approximately values of the
# drag and reflectivity coefficients, respectively.

molniyaMass <- 1600
molniyaCrossSection <- 15
molniyaCr <- 1.2
molniyaCd <- 2.2

# Let´s calculate the position and velocity of the satellite for each minute of
# the following 10 minutes.

targetTimes <- seq(0, 600, by=60)
hpop_results <- hpop(initialPosition, initialVelocity, initialTime, targetTimes, 
                     molniyaMass, molniyaCrossSection, molniyaCrossSection,
                     molniyaCr, molniyaCd)
}
}
