\name{readBinSPK}
\alias{readBinSPK}
\title{Read a binary SPK file
}
\description{
SPK (Spacecraft and Planet Kernel) is a binary file format developed by NAIF to 
store ephemerides (trajectory) of celestial bodies and spacecraft in the Solar 
System. The file format is based on the DAF architecture (see \link{readBinDAF}).
A detailed description of the SPK file format can be found in NAIF's documentation
(https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html).

Each SPK file contains several segments, with each segment comprising a summary
(or descriptor), a name and an array of double precision elements. Each segment
is conceptually equivalent to an array in the context of generic DAF files.
There are several types of SPK segments defined by NAIF, each identified by an SPK type
code currently ranging from 1 to 21 (some intermediate values are not used or not
available for general public use). Each segment type provides ephemerides information
in a different way. Note that the segments stored in a single SPK file can be
of different types. A detailed description of the organization of the arrays for
each SPK type can be found at https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html

This function allows to read SPK binary files of all types except 4, 6, 7 and 16.
The data will be presented properly formatted and the meaning of each element is
assigned during the reading process of the file. It should be noted that this 
function just reads SPK kernels; it does not provide any evaluation of ephemerides
at arbitrary target times. The process of performing such evaluation differs between
SPK types. For example, types 2, 3, 14 and 20 require Chebyshev interpolation;
typees 8 and 9 require Lagrange interpolation; type 10 requires the aplication of
SGP4/SDP4, etc. Nevertheless, it is still possible to obtain direct ephemerides 
information just by reading the SPK kernels since many of the types often include
reference state vectors between which interpolation is applied, but that can be
directly used at the epochs corresponding to said reference state vectors.
}
\usage{
readBinSPK(filename)
}
\arguments{
  \item{filename}{Path to the binary SPK file.
  }
}
\value{
A list with two elements. The first element, named \code{comments}, is a character 
vector where each element is a line of comments.

The second element is named \code{segments}, and is a nested list where each top-level
element represents one of the segments stored in the SPK file and its associated 
metadata. Each of the top-level elements is itself a list with the following 3 
elements:
  \item{segmentName }{String with the name of the segment}
  \item{segmentSummary }{A list with the multiple doubles and integers that are stored
  in each array summary of summary records and which provide metadata describing
  each array. In the case of SPK files, these are always the following 9 elements:
  \item{SPKType }{A description of the type of SPK segment}
  \item{initialEpoch }{The initial epoch for the interval for which ephemeris data
  are provided in the segment, in ephemeris seconds (equivalent to TDB seconds) 
  past Julian year 2000}
  \item{finalEpoch }{The initial epoch for the interval for which ephemeris data
  are provided in the segment, in ephemeris seconds (equivalent to TDB seconds) 
  past Julian year 2000}
  \item{targetNAIFCode }{The NAIF integer code for the object for which the segment
  provides ephemerides. For details, see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
  }
  \item{targetNAIFCode }{The NAIF integer code for the central body of the reference
  frame in which the segment provides ephemerides. For details, see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
  }
  \item{frameNAIFCode }{The NAIF frame code for the reference frame in which the 
  segment provides ephemerides. For details, see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
  }
  \item{SPKType }{The SPK type code for the segment}
  \item{initialArrayAddress }{The initial address of the array elements corresponding 
  to this segment within the SPK file, in double precision numbers (in order to obtain
  byte address, multiply by 8 and subtract 7)}
  \item{finalArrayAddress }{The final address of the array elements corresponding 
  to this segment within the SPK file, in double precision numbers (in order to obtain
  byte address, multiply by 8 and subtract 7)}
  }
  \item{segmentData }{A list with the actual ephemeris data contained in the segment,
  as well as some type-specific additional metadata}

The contents of the last element, \code{segmentData}, are different for each SPK
type. Here a summary is provided for each one, but for more detailed descriptions
see NAIF's documentation.

For type 1, which provide Modified Difference Arrays (MDAs), a list where each
element is one of the records of the segment. Each of these elements contains
the following elements:
  \item{finalEpoch }{The final epoch (and also reference point) for which this 
  MDA should be used to compute a state vector}
  \item{stepsizeFunctionVector }{A vector of differences between the final epoch/reference
  point and the epochs of all the other data points used to fit the interpolation
  polynomials from which the MDA is derived. This is of length 15}
  \item{referencePosition }{A numeric vector of length 3 containing the X, Y and Z
  components of the position at the final/reference epoch, in km}
  \item{referencePosition }{A numeric vector of length 3 containing the X, Y and Z
  components of the velocity at the final/reference epoch, in km/s}
  \item{MDA }{A matrix with 15 rows and 3 columns providing the constant coefficients
  to interpolate position, velocity and acceleration at a target epoch. The 1st, 
  2nd and 3rd columns give the coefficients for the X, Y and Z components respectively. 
  The given coefficients basically are the coefficients of the interpolation polynomial
  when expressed in Modified Divided Differences (MDD) form. For details, see
  Shampine and Gordon, 1975. Note that even though the matrix will always have
  15 rows (corresponding to 15 coefficients for each components, or an interpolation
  polynomial degree of 14), some of these can have a value of 0 up to the 15th row, 
  effectively leading to an interpolation polynomial of degree lower than 14.}
  \item{maxIntegrationOrderP1 }{The maximum order for the interpolation polynomial
  that should be applied amongst all 3 components plus 1. For example, if the interpolation
  polynomial orders for the X, Y and Z components are 6, 8 and 7 respectively, this
  element will have a value of 9 (\code{(max(c(6,8,7))+1)}).}
  \item{integrationOrderArray }{A vector of 3 integers indicating the order of
  the interpolation polynomials that should be applied for the X, Y and Z components
  of acceleration respectively.}
A brief description is provided here for the meaning of the coefficients, based
on the 'SPICE spke01 math' monograph by Robert Werner.
MDAs are a modified version of the more standard coefficients obtained through the
divided differences method, which represent an interpolation polynomial in Newton
form. In order to calculate a position and velocity, first the basis functions
must be computed for both position and velocity. Then each of the components of 
the basis functions (each of an increasing degree) must be multiplied by the 
corresponding MDA coefficient, and all the resulting terms are summed.
  
For type 2, which provide Chebyshev coefficients for position only and at equally 
spaced time steps, a list with the following elements:
  \item{polynomialDegree }{An integer indicating the order of the interpolation 
  polynomial that should be applied for all the components.}
  \item{chebyshevCoefficients }{A matrix where each row corresponds to an interpolation
  interval, and with the following columns: 
  \item{initialEpoch}{Initial epoch of the interpolation intervals, in ephemeris 
  (TDB) seconds since J2000}
  \item{midPoint}{Epoch for the midpoint of the interpolation intervals, in ephemeris 
  (TDB) seconds since J2000}
  \item{intervalRadius}{Radius of the interpolation intervals, in seconds)}
  \item{positionXCoeffi}{A set of N columns, with i ranging from 1 to N, providing
  the Chebyshev coefficients for the X component of the position. N is the number 
  of coefficients for each component, which is equal to \code{polynomialDegree + 1}}
  \item{positionYCoeffi}{As \code{positionXCoeffi}, but for the Y component of position.}
  \item{positionZCoeffi}{As \code{positionXCoeffi}, but for the Z component of position.}
  }

For type 3, which provide Chebyshev coefficients for position and velocity and at 
equally spaced time steps, the same as for type 2, but \code{chebyshevCoefficients}
contains the following additional elements:
  \item{velocityXCoeffi}{A set of N columns, with i ranging from 1 to N, providing
  the Chebyshev coefficients for the X component of the velocity. N is the number 
  of coefficients for each component, which is equal to \code{polynomialDegree + 1}}
  \item{velocityYCoeffi}{As \code{velocityXCoeffi}, but for the Y component of velocity.}
  \item{velocityZCoeffi}{As \code{velocityXCoeffi}, but for the Z component of velocity.}

For type 5, which provide discrete state vectors to be propagated following the
laws of two-body motion, a list with the following elements:
  \item{centralBodyGM }{The GM parameter (gravitational constant) for the central body, in
  cube kilometers per square seconds}
  \item{stateVectors }{A matrix where each row corresponds to a state vector, and
  with the following columns: 
  \item{epoch}{Epoch of the state vectors intervals, in ephemeris (TDB) seconds since J2000}
  \item{positionX}{X component of the position, in km}
  \item{positionY}{Y component of the position, in km}
  \item{positionZ}{Z component of the position, in km}
  \item{velocityX}{X component of the velocity, in km/s}
  \item{velocityY}{Y component of the velocity, in km/s}
  \item{velocityZ}{Z component of the velocity, in km/s}
  }

For type 8, which provide discrete state vectors at equally spaced time steps to 
which Lagrange interpolation should be applied to obtain state vectors at arbitrary 
target times, a list with the following elements:
  \item{polynomialDegree }{The degree of the interpolation polynomial that should
  be applied}
  \item{stateVectors }{A matrix where each row corresponds to a state vector, and
  with the following columns: 
  \item{epoch}{Epoch of the state vectors intervals, in ephemeris (TDB) seconds since J2000}
  \item{positionX}{X component of the position, in km}
  \item{positionY}{Y component of the position, in km}
  \item{positionZ}{Z component of the position, in km}
  \item{velocityX}{X component of the velocity, in km/s}
  \item{velocityY}{Y component of the velocity, in km/s}
  \item{velocityZ}{Z component of the velocity, in km/s}
  }

For type 9, which provide discrete state vectors at unequally spaced time steps to 
which Lagrange interpolation should be applied to obtain state vectors at arbitrary 
target times, the same as for type 8.

For type 10, which provide TLE that should be propagated with SGP4/SDP4, a list
with the following elements:
  \item{constants }{A list with the following constants used by SGP4/SDP4:
  \item{J2}{J2 parameter (dimensionless)}
  \item{J3}{J3 parameter (dimensionless)}
  \item{J4}{J4 parameter (dimensionless)}
  \item{sqrtGM}{Square root of the GM parameter, where GM is in cube Earth radii
  per square minutes}
  \item{highAltBound}{High altitude boundary for atmospheric model (in km)}
  \item{lowAltBound}{Low altitude boundary for atmospheric model (in km)}
  \item{earthRadius}{Equatorial radius of Earth (in km)}
  \item{distUnitsPerRadius}{Distance units per Earth radius. This is usually 1.
  If different than 1, interprete results with caution}
  }
  \item{TLEs }{A matrix where each row corresponds to a TLE, and with the following 
  columns (nutation angles and their rates are not present in some old SPK files, 
  in which case they will have NULL values): 
  \item{epoch}{Epoch of the TLE, in ephemeris (TDB) seconds since J2000}
  \item{meanMotionDerivative}{First derivative of the mean motion in radians/min^2}
  \item{meanMotionSecondDerivative}{Second derivative of the mean motion in radians/min^3}
  \item{Bstar}{Drag coefficient of the satellite in units of (earth radii)^-1. 
  Bstar is an adjusted value of the ballistic coefficient of the satellite, and 
  it indicates how susceptible it is to atmospheric drag.}
  \item{inclination}{Mean orbital inclination of the satellite in radians}
  \item{ascension}{Mean longitude of the ascending node of the satellite at epoch, 
  also known as right ascension of the ascending node, in radians}
  \item{eccentricity}{Mean eccentricity of the orbit of the object}
  \item{perigeeArgument}{Mean argument of the perigee of the object in radians} 
  \item{meanAnomaly}{Mean anomaly of the orbit of the object in radians}
  \item{meanMotion}{Mean motion of the satellite at epoch in radians/min}
  \item{deltaPsi}{Obliquity (psi angle) of the nutation at epoch, in radians}
  \item{deltaEpsilon}{Longitude (epsilon angle) of the nutation at epoch, in radians}
  \item{deltaPsiDerivative}{Derivative of the obliquity (psi angle) of the nutation 
  at epoch, in radians/second}
  \item{deltaEpsilonDerivative}{Derivative of the longitude (epsilon angle) of the 
  nutation at epoch, in radians/second}
  }

For type 12, which provide discrete state vectors at equally spaced time steps to 
which Hermite interpolation should be applied to obtain state vectors at arbitrary 
target times, the same as for type 8, but the list contains an additional element:
  \item{windowSize}{The window size that should be applied during interpolation}
  
For type 13, which provide discrete state vectors at unequally spaced time steps to 
which Hermite interpolation should be applied to obtain state vectors at arbitrary 
target times, the same as for type 12.

For type 14, which provide Chebyshev coefficients for position and velocity and at 
unequally spaced time steps, the same as for type 3.

For type 15, which provide elements for calculation of ephemerides through the application
of a precessing conic propagation model, a list with the following elements:
  \item{epochPeriapsis}{Epoch of the periapsis passage in ephemeris (TDB) seconds since J2000}
  \item{unitVectorTrajectoryPoleX}{X component of the unit trajectory pole vector, in km} 
  \item{unitVectorTrajectoryPoleY}{Y component of the unit trajectory pole vector, in km}
  \item{unitVectorTrajectoryPoleZ}{Z component of the unit trajectory pole vector, in km}
  \item{unitVectorPeriapsisX}{X component of the unit periapsis vector, in km} 
  \item{unitVectorPeriapsisY}{Y component of the unit periapsis vector, in km}
  \item{unitVectorPeriapsisZ}{Z component of the unit periapsis vector, in km}
  \item{semiLatusRectum}{Semi-latus rectum, in km} 
  \item{eccentricity}{Eccentricity of the orbit}
  \item{J2ProcessingFlag}{Flag indicating what J2 corrections should be applied
  when propagating. If 1, regress line of nodes only. If 2, precess line of apsides
  only. If 3, don't use any corrections. For any other values, regress line of nodes
  and precess line of apsides}
  \item{unitVectorCentralBodyPoleX}{X component of the unit central body pole vector, in km} 
  \item{unitVectorCentralBodyPoleY}{Y component of the unit central body pole vector, in km} 
  \item{unitVectorCentralBodyPoleZ}{Z component of the unit central body pole vector, in km} 
  \item{centralBodyGM}{The GM parameter (gravitational constant) for the central body, in
  cube kilometers per square seconds}
  \item{centralBodyJ2}{The J2 parameter for the central body (dimensionless)}
  \item{centralBodyRadius}{Radius of the central body, in km}

For type 17, which provide equinoctial elements modelling an object following an
elliptic orbit with precessing line of nodes and argument of periapse relative 
to the equatorial frame of a central body, a list with the following elements:
  \item{epochPeriapsis}{Epoch of the periapsis passage in ephemeris (TDB) seconds since J2000}
  \item{semiMajorAxis}{Semi-major axis of the orbit, in km} 
  \item{equinoctialH}{Value of the equinoctial parameter H at epoch} 
  \item{equinoctialK}{Value of the equinoctial parameter K at epoch} 
  \item{meanLongitude}{Mean longitude of the orbit at epoch, in radians} 
  \item{equinoctialP}{Value of the equinoctial parameter P} 
  \item{equinoctialQ}{Value of the equinoctial parameter Q} 
  \item{longitudePeriapsisDerivative}{Derivative of the longitude of periapse at 
  epoch (but it is assumed to be constant at other times), in radians/s} 
  \item{meanLongitudeDerivative}{Derivative of the mean longitude at epoch
  (but it is assumed to be constant at other times), in radians/s} 
  \item{longitudeAscendingNodeDerivative}{Derivative of the longitude of the ascending
  node at epoch, in radians/s}
  \item{equatorialPoleRightAscension}{Right ascension of the pole of the orbital 
  reference system relative to the reference frame of the corresponding SPK segment,
  in radians}
  \item{equatorialPoleDeclination}{Declination of the pole of the orbital reference
  system relative to the reference frame of the corresponding SPK segment, in radians}

For type 18, which provide ephemerides in the format used by ESA on the Mars Express,
Rosetta, SMART-1 and Venus Express missions (although applicable to any other object),
there are 2 different subtypes: subtype 0 and subtype 1. Subtype 0 should be used
to perform sliding-window Hermite interpolation of position and velocity independently.
Subtype 1 should be used to perform sliding-window Lagrange interpolation of
position and velocity independently. In both cases, \code{segmentData} is a list 
with the following elements:
  \item{subTypeCode}{Subtype code for the type 18 SPK segment}
  \item{polynomialDegree }{An integer indicating the order of the interpolation 
  polynomial that should be applied for all the components.}
  \item{interpolationType }{Type of the interpolation that should be applied.
  Hermite for subtype 0, and Lagrange for subtype 1} 
  \item{windowSize}{The window size that should be applied during interpolation}
  \item{meanLongitude}{Mean longitude of the orbit at epoch, in radians} 
  \item{stateVectors}{A matrix where each row corresponds to a state vector. The
  columns differ depending on the subtype. The following will always be present:
  \item{epoch}{Epoch of the state vectors intervals, in ephemeris (TDB) seconds since J2000}
  \item{positionX}{X component of the position, in km}
  \item{positionY}{Y component of the position, in km}
  \item{positionZ}{Z component of the position, in km}
  The following are present for subtype 0:
  \item{firstVelocityX}{X component of the first velocity value, in km/s. The first
  velocity value should be used together with the reference position to interpolate
  position values}
  \item{firstVelocityY}{Y component of the first velocity value, in km/s}
  \item{firstVelocityZ}{Z component of the first velocity value, in km/s}
  \item{secondVelocityX}{X component of the second velocity value, in km/s. The second
  velocity value should be used together with the reference acceleration to interpolate
  velocity values}
  \item{secondVelocityY}{Y component of the second velocity value, in km/s}
  \item{secondVelocityZ}{Z component of the second velocity value, in km/s}
  \item{accelerationX}{X component of the acceleration, in km/s^2}
  \item{accelerationY}{Y component of the acceleration, in km/s^2}
  \item{accelerationZ}{Z component of the acceleration, in km/s^2}
  The following are present for subtype 1:
  \item{velocityX}{X component of the velocity, in km/s}
  \item{velocityY}{Y component of the velocity, in km/s}
  \item{velocityZ}{Z component of the velocity, in km/s}
  }

For type 19, which provides the same data as type 18 but condensing multiple
type 18 segments into a single 19 segments (only possible if all the segments
have the same target object, central body and reference frame; additionally,
the coverage of the segments must overlap only at endpoints and leave no gaps),
a list with the following 2 elements:
  \item{boundaryChoiceFlag}{A flag indicating which minisegment should be used
  for the epochs at which 2 minisegments overlap. If 0, the earlier minisegment 
  that ends at that epoch is used. If 1, the later minisegment that begins at
  that epoch is used}
  \item{minisegments}{A nested list where each top-level element represents a type 19 
  subsegment (called minisegments in NAIF's documentation), each of which is a list 
  with the same elements as a \code{segmentData} for a type 18 SPK segment, plus 
  the following 2 additional elements:
  \item{intervalStartEpoch}{Beginning of the interpolation interval covered by 
  this minisegment, in ephemeris (TDB) seconds since J2000}
  \item{intervalEndEpoch}{End of the interpolation interval covered by this minisegment, 
  in ephemeris (TDB) seconds since J2000}
  Furthermore, a third subtype can be found in type 19 minisegments which is not
  found for type 18 segments (at least according to NAIF's documentation). This
  has subtype code 2, and should be used to perform sliding-window Hermite 
  interpolation of position and velocity together (note the difference with subtype 0,
  where position and velocity are interpolated independently). In this case,
  element \code{interpolationType} has a value of "Hermite-joint", and the element
  describing the minisegment contains the same elements as a type 18 SPK segment of
  subtype 1 (and not of subtype 0).
  }

For type 20, which provide which provide Chebyshev coefficients for velocity only 
and at equally spaced time steps together with a reference position (so that
position can be interpolated by integration of velocity), a list with the following
elements:
  \item{polynomialDegree }{An integer indicating the order of the interpolation 
  polynomial that should be applied for all the components.}
  \item{dScale }{Distance scale used for both position and velocity, in km. For 
  example, if dScale has a value of 149597870.7 (the length of an astronomical
  unit, AU, in km), it means the distance units are AU.}
  \item{tScale }{Time scale used for velocity, in TDB seconds. For example, a value 
  of 1 means we are using velocity directly in TDB seconds. A value of 86400 means 
  the time units of velocity would be TDB Julian days, etc.}
  \item{chebyshevCoefficients }{A matrix where each row corresponds to an interpolation
  interval, and with the following columns: 
  \item{initialEpoch}{Initial epoch of the interpolation intervals, in ephemeris 
  (TDB) seconds since J2000}
  \item{midPoint}{Epoch for the midpoint of the interpolation intervals, in ephemeris 
  (TDB) seconds since J2000}
  \item{intervalRadius}{Radius of the interpolation intervals, in seconds)}
  \item{velocityXCoeffi}{A set of N columns, with i ranging from 1 to N, providing
  the Chebyshev coefficients for the X component of the velocity. N is the number 
  of coefficients for each component, which is equal to \code{polynomialDegree + 1}}
  \item{velocityYCoeffi}{As \code{velocityXCoeffi}, but for the Y component of velocity.}
  \item{velocityZCoeffi}{As \code{velocityXCoeffi}, but for the Z component of velocity.}
  \item{midPointPositionX}{X component of the reference position, valid at the midpoint
  of the interpolation interval. Should be used to interpolate position through 
  integration.}
  \item{midPointPositionY}{As \code{midPointPositionX}, but for the Y component of velocity.}
  \item{midPointPositionZ}{As \code{midPointPositionX}, but for the Z component of velocity.}
  }

For type 21, which, like type 1, also provide MDAs, the same nested list as for 
type 1, with 2 differences. Firstly, unlike for type 1, MDAs of type 21 segments 
are not limited to a maximum of 15 coefficients per component. Second, each top-level 
element of the nested list contains an additional element:
  \item{numberCoefficients}{The number of coefficients provided for the component
  with the highest interpolation order}
}
\references{
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
Shampine, L. F. and Gordon, M. K., Computer Solution of Ordinary Differential Equations: The Initial Value Problem, 1975
Robert Werner, SPICE spke01 math, 2022.
}
\examples{
# The file vgr2_jup230.bsp provided with the package includes information for the
# Jupiter flyby of Voyager 2

testSPK <- readBinSPK(paste0(path.package("asteRisk"), "/vgr2_jup230.bsp"))
length(testSPK$segments)
# It contains a single segment
testSPK$segments[[1]]$segmentSummary$SPKType
testSPK$segments[[1]]$segmentSummary$SPKTypeCode
# The segment is of type 1, containing Modified Difference Arrays
length(testSPK$segments[[1]]$segmentData)
# It contains 566 MDAs
}
