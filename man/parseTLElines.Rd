\name{parseTLElines}
\alias{parseTLElines}
\title{Parse the lines of a TLE
}
\description{
TLE (Two-/Three- Line Element) is the standard format for representing orbital
state vectors. This function parses a character vector where each element
represents a line of the TLE. The supplied character vector can have either 2
(for Two Line Elements) or 3 (for Three Line Elements) elements. The two lines
of a Two Line Element contain all the information. The additional line in a
Three Line Element is optional, and contains just the satellite name. For a
detailed description of the TLE format, see https://celestrak.org/columns/v04n03/#FAQ01.
}
\usage{
parseTLElines(lines)
}
\arguments{
  \item{lines}{Character vector where each element is a string corresponding to
  a line of the TLE. The character vector must have either 2 or 3 elements.
  }
}
\value{
A list with the following elements that define the orbital state vector of the
satellite:
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
https://celestrak.org/columns/v04n03/#FAQ01
}
\examples{
# The following lines correspond to a TLE for Italsat 2 the 26th of June, 2006
# at 00:58:29.34 UTC.

italsat2_lines <- c("ITALSAT 2",
"1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600",
"2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119")

italsat2_TLE <- parseTLElines(italsat2_lines)
italsat2_TLE
}
