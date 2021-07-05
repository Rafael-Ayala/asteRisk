\name{readTLE}
\alias{readTLE}
\title{Read a TLE file
}
\description{
TLE (Two-/Three- Line Element) is a standard format for representing orbital
state vectors. This function reads a TLE file containing one or more TLEs. The
TLE file can contain either Two Line Elements or Three Line Elements, but all
the TLE in a single file must be of the same type. The two lines of a Two Line
Element contain all the ephemeris information. The additional line in a Three
Line Element is optional, and contains just the satellite name. For a detailed
description of the TLE format, see https://celestrak.com/columns/v04n03/#FAQ01.
}
\usage{
readTLE(filename)
}
\arguments{
  \item{filename}{Path to the TLE file.
  }
}
\value{
A list with the following elements that define the orbital state vector of the
satellite (or, if the file contained multiple TLE, a nested list, where each
element of the top level list represents an orbital state vector, and comprises
the following elements):
  \item{NORADcatalogNumber }{NORAD Catalog Number, also known as Satellite
  Catalog Number, assigned by United States Space Command to each artificial
  object orbiting Earth}
  \item{classificationLevel }{Classification level of the information for the
  orbiting object. Can be unclassified, classified, secret or unknown}
  \item{internationalDesignator }{International Designator, also known as COSPAR
  ID, of the object. It consists of the launch year, separated by a hyphen from
  a three-digit number indicating the launch number for that year and a set of
  one to three letters indicating the piece for a launch with multiple pieces.}
  \item{launchYear }{The launch year of the object}
  \item{launchNumber }{The launch number of the object during its launch year}
  \item{launchPiece }{The piece for the launch of the object, if it was a launch
  with multiple pieces}
  \item{dateTime }{Date time string to which the orbital state vector
  corresponds}
  \item{elementNumber }{Element number for the object. In principle, every time
  a new TLE is generated for an object, the element number is incremented, and
  therefore element numbers could be used to assess if all the TLEs for a
  certain object are available. However, in practice it is observed that this is
  not always the case, with some numbers skipped and some numbers repeated.}
  \item{inclination }{Mean orbital inclination of the satellite in degrees.
  This is the angle between the orbital plane of the satellite and the
  equatorial plane}
  \item{ascension }{Mean longitude of the ascending node of the satellite at
  epoch, also known as right ascension of the ascending node, in degrees. This
  is the angle between the direction of the ascending node (the point where the
  satellite crosses the equatorial plane moving north) and the direction
  of the First Point of Aries (which indicates the location of the vernal
  equinox)}
  \item{eccentricity }{Mean eccentricity of the orbit of the object.
  Eccentricity is a measurement of how much the orbit deviates from a circular
  shape, with 0 indicating a perfectly circular orbit and 1 indicating an
  extreme case of parabolic trajectory}
  \item{perigeeArgument }{Mean argument of the perigee of the object in degrees.
  This is the angle between the direction of the ascending node and the
  direction of the perigee (the point of the orbit at which the object is
  closest to the Earth)}
  \item{meanAnomaly }{Mean anomaly of the orbit of the object in degrees. This
  indicates where the satellite is along its orbital path. It is provided as
  the angle between the direction of the perigee and the hypothetical point
  where the object would be if it was moving in a circular orbit with the same
  period as its true orbit after the same amount of time since it last crossed
  the perigee had ellapsed. Therefore, 0 denotes that the object is at the
  perigee}
  \item{meanMotion}{Mean motion of the satellite at epoch in revolutions/day
  }
  \item{meanMotionDerivative}{First time derivative of the mean motion of the
  satellite in revolutions/day^2^}
  \item{meanMotionSecondDerivative}{Second time derivative of the mean motion of
  the satellite in revolutions/day^3^.
  }
  \item{Bstar}{Drag coefficient of the satellite in units of (earth radii)^-1^.
  Bstar is an adjusted value of the ballistic coefficient of the satellite, and
  it indicates how susceptible it is to atmospheric drag.
  }
  \item{ephemerisType}{Source for the ephemeris (orbital state vector). Most
  commonly, it is distributed data obtained by combaining multiple observations
  with the SGP4/SDP4 models
  }
  \item{epochRevolutionNumber}{Number of full orbital revolutions completed by
  the object
  }
  \item{objectName}{Name of the object, retrieved from the first line of the
  TLE if a Three Line Element was provided
  }
}
\references{
https://celestrak.com/columns/v04n03/#FAQ01
http://www.celestrak.com/publications/aiaa/2006-6753/AIAA-2006-6753.pdf
}
\examples{
# The file testTLE.txt provided with the package includes 29 TLE covering a
# variety of satellites, extracted from Revisiting Space Track Report #3

test_TLEs <- readTLE(paste0(path.package("asteRisk"), "/testTLE.txt"))
test_TLEs
}
