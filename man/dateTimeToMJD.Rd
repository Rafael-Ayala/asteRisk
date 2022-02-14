\name{dateTimeToMJD}
\alias{dateTimeToMJD}
\title{Calculate Modified Julian Date for a given date and time
}
\description{
The Julian Date (JD) of a given date and time is the number of days since noon of
Monday 1st of January 4713 BC, including a fractional part. Modified Julian Date
(MJD) are instead the number of days since 00:00 of November 17th, 1858. The 
difference JD and MJD for a given instant is always 2400000.5, which is the JD
of the reference time for MJD.

This function calculates the MJD of a date and time, provided as a date-time
character string in UTC time. The output refers by default to the MJD in UTC,
but different time systems can be chosen: UTC (Coordinated Universal Time),
UT1 (Universal Time), TT (Terrestrial Time) and TDB (Barycentric Dynamical Time).
}
\usage{
dateTimeToMJD(dateTime, timeSystem="UTC")
}
\arguments{
  \item{dateTime}{Date-time string with the date and time in UTC corresponding
  to the provided geodetic coordinates.
  }
  \item{timeSystem}{Time system into which the MJD should be calculated. Should
  be one from "UTC" (Coordinated Universal Time; default), "UT1" (Universal Time), 
  "TT" (Terrestrial Time) and "TDB" (Barycentric Dynamical Time).
  }
}
\value{
The MJD for the specified date and time in the chosen time system.
}
\references{
https://gssc.esa.int/navipedia/index.php/Julian_Date
https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems
}
\examples{
if(requireNamespace("asteRiskData", quietly = TRUE)) {
# Let's calculate the MJD of the 12th of June, 2000 at 10:00:00 UTC time, in UTC

MJD_UTC <- dateTimeToMJD("2000-06-12 10:00:00", timeSystem = "UTC")

# Let's now calculate the MJD for the same instant in TDB:

MJD_TDB <- dateTimeToMJD("2000-06-12 10:00:00", timeSystem = "TDB")

# We can now calculate the difference in seconds, which matches the difference 
# between UTC and TDB for that day:

(MJD_UTC - MJD_TDB) * 86400
}
}
