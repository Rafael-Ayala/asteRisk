\name{revDay2radMin}
\alias{revDay2radMin}
\title{Converts revolutions per day to radians per minute
}
\description{
This function converts a rotation rate in revolutions per day to radians per
minute. This conversion is useful since values in TLEs are given in revolutions
per day, but the SGP4 and SDP4 propagators require the mean motion to be provided
in radians per minute.
}
\usage{
revDay2radMin(revPerDay)
}
\arguments{
  \item{revPerDay}{Value of the rotation rate in revolutions per day.
  }
}
\value{
The corresponding value of the rotation rate in radians per minute.
}
\examples{
revDay2radMin(2)
}
