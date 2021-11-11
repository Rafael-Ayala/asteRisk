\name{hpop}
\alias{hpop}
\title{High-precision numerical orbital propagator
}
\description{
Given the position and velocity of a satellite at a given time (in the GCRF
system of coordinates), propagates its position by calculating its acceleration 
(based on a force model) and solving the resulting second-order ODE through 
numerical integration. This allows propagation of orbits with considerably 
higher accuracy than other propagators such as SGP4 and SDP4, but at the expense 
of a much higher computational cost. The forces and effects currently considered 
are gravitational attraction by the  Earth (using a geopotential model based on 
spherical harmonics); effects of Earth ocean and solid tides; gravitational 
attraction by the Moon, Sun and planets (considered as point masses); solar 
radiation pressure; atmospheric drag, and relativistic effects. The force
field is based on the forces described in Satellite Orbits: Models, Methods and 
Applications (Oliver Montenbruck and Eberhard Gill) and Fundamentals of
Astrodynamics and Applications (David Vallado). Parts of this implementation
are based on a previous MATLAB implementation by Meysam Mahooti.
The NRLMSISE-00 model is used to calculate atmospheric density for the
calculation of atmospheric drag.
The high-precision numerical orbital propagator requires the \code{asteRiskData}
package, which provides the data and coefficients required for calculation of
the modeled forces. \code{asteRiskData} can be installed by running 
\code{install.packages('asteRiskData', repos='https://rafael-ayala.github.io/drat/')}
}
\usage{
hpop(position, velocity, dateTime, times, satelliteMass, dragArea, 
     radiationArea, dragCoefficient, radiationCoefficient, ...)
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
  \item{...}{Additional parameters to be passed to \link{ode} to control how
  numerical integration is performed. By default, the RADAU5 solver is used.
  }
}
\value{
A matrix with the results of the numerical integration at the requested times.
Each row contains the results for one of the requested times. The matrix 
contains seven columns: time (indicating the time for the corresponding row
in seconds since the initial time), X, Y, Z (indicating the X, Y and Z
components of the position for that time in meters), dX, dY and dZ (indicating
the X, Y and Z components of the velocity for that time in meters/second).
Positions and velocities are returned in the GCRF frame of reference.
}
\references{
Satellite Orbits: Models, Methods and Applications. Oliver Montenbruck and Eberhard Gill.
Fundamentals of Astrodynamics and Applications. David Vallado.
https://www.mathworks.com/matlabcentral/fileexchange/55167-high-precision-orbit-propagator
https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1144&context=smallsat
https://iopscience.iop.org/article/10.1088/1742-6596/911/1/012009/pdf
https://www.sciencedirect.com/science/article/pii/S1110016821000016
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
