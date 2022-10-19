\name{getLatestSpaceData}
\alias{getLatestSpaceData}
\title{Retrieves the latest space data
}
\description{
The \code{asteRiskData} package provides the data and coefficients required
for calculation of forces for \link{hpop} and other functions such certain
conversions between reference frames. Some of the data tables included in the
package are updated periodically with new data. These include Earth orientation 
parameters, space weather data and solar and geomagnetic storms. In order to 
perform the calculations dependent on such data for the most recent times, the 
latest available data must be retrieved.

This function automatically updates the data tables, enabling such calculations
for the most recent times.
}
\usage{
getLatestSpaceData(targets="all")
}
\arguments{
  \item{targets}{Character vector specifying the data that should be updated.
  It should be a vector containing one or more of the following strings: "all"
  (to update all data), "EOP" (Earth orientation parameters), "SW" (space weather),
  "SS" (solar storms) or "GS" (geomagnetic storms). By default, all data are 
  updated.
  }}
\value{
This function is invoked for its side effect, which is updating the data tables
used internally for calculations requiring \code{asteRiskData} package, such as 
those performed by \link{hpop}.
}
\references{
http://www.celestrak.org/SpaceData/EOP-All.txt
https://celestrak.org/SpaceData/SW-All.txt
https://sol.spacenvironment.net/jb2008/indices.html
}
\examples{
if(interactive()) {
if(requireNamespace("asteRiskData", quietly = TRUE)) {
# The table of Earth orientation parameters distributed with asteRiskData
# comprises data up to the 21st of March, 2021

asteRiskData::earthPositions[nrow(asteRiskData::earthPositions),]

# The table can be easily updated to include the most recent available data

getLatestSpaceData(targets="all")
asteRiskData::earthPositions[nrow(asteRiskData::earthPositions),]
}
}
}
