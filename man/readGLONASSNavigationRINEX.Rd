\name{readGLONASSNavigationRINEX}
\alias{readGLONASSNavigationRINEX}
\title{Read a RINEX navigation file for GLONASS satellites
}
\description{
RINEX (Receiver Independent Exchange Format) is one of the most widely used
formats for providing data of satellite navigation systems. The RINEX standard
defines several file types, among which navigation files are used to distribute
positional information of the satellites. The exact information provided in a
RINEX navigation file varies for each satellite navigation system. This function
reads RINEX navigation files for satellites of the GLONASS constellation,
operated by Russia.
}
\usage{
readGLONASSNavigationRINEX(filename)
}
\arguments{
  \item{filename}{Path to the GLONASS RINEX navigation file.
  }
}
\value{
A list with two elements. The first element, named \code{header}, is a list with
the following elements:
  \item{rinexVersion }{Version of the RINEX format used in the file}
  \item{rinexFileType }{Type of RINEX file}
  \item{generatorProgram }{Program used to generate the RINEX file}
  \item{generatorEntity }{Individual or organization that generated the file}
  \item{fileCreationDateString }{Date-time string indicating when the file was
  created}
  \item{refYear }{Reference year for system time correction}
  \item{refMonth }{Reference month for system time correction}
  \item{refDay }{Reference day for system time correction}
  
  \item{sysTimeCorrection }{Correction to system time scale to fine-tune GLONASS
  time to UTC in seconds. Since GLONASS time is linked to UTC, it should be a
  very small amount}
  \item{leapSeconds }{Leap seconds introduced since 1980. Useful to convert to
  GPS time}
  \item{comments }{Miscellaneous comments found in the header of the RINEX file}

The second element is named \code{messages}, and it contains one element for
each navigation message found in the RINEX file. Each of these elements is a 
list with the following elements that provide information about the position
of the GLONASS satellite:
  \item{satelliteNumber }{Slot number of the satellite within the GLONASS
  constellation. It can be converted to a PRN code by adding 37 to it}
  \item{epochYearShort }{Epoch year in 2-digit format. If lower than 80, the
  meaning should be taken as 20XX, while if larger than 80, it refers to 19XX.}
  \item{epochMonth }{Epoch month}
  \item{epochDay }{Epoch day}
  \item{epochHour }{Epoch hour}
  \item{epochMinute }{Epoch minute}
  \item{epochSecond }{Epoch second}
  \item{ephemerisUTCTime }{A \link{nanotime} object indicating the time
  corresponding to the reported position (ephemeris) in the present message. 
  The time is in UTC, obtained by applying the individual clock bias of the 
  particular satellite (\code{clockBias} field of each message) and the latest 
  global GLONASS time bias with respect to UTC (\code{sysTimeCorrection} field
  of the header) to the uncorrected ephemeris time, given by the previous
  time fields. Corrections are performed as described in the GLONASS system
  specifications (https://russianspacesystems.ru/wp-content/uploads/2016/08/ICD_GLONASS_eng_v5.1.pdf)}
  \item{clockBias }{Clock bias (i.e., constant offset) that should be applied to
  the satellite time in order to obtain an even more accurate UTC time. In 
  seconds}
  \item{relativeFreqBias }{Clock drift of the satellite clock that should be
  applied in combination with the time difference to the reference time in order
  to obtain an even more accurate UTC time. In seconds per second}
  \item{messageFrameTime }{Second of the UTC day when the message was 
  transmitted}
  \item{positionX }{X coordinate of the position of the satellite in km, in the
  ECEF system of coordinates}
  \item{positionY }{Y coordinate of the position of the satellite in km, in the
  ECEF system of coordinates}
  \item{positionZ }{Z coordinate of the position of the satellite in km, in the
  ECEF system of coordinates}
  \item{velocityX }{X component of the velocity of the satellite in km/s, in the
  ECEF system of coordinates}
  \item{velocityY }{Y component of the velocity of the satellite in km/s, in the
  ECEF system of coordinates}
  \item{velocityZ }{Z component of the velocity of the satellite in km/s, in the
  ECEF system of coordinates}
  \item{accelX }{X component of the accel of the satellite in km/s, in the
  ECEF system of coordinates}
  \item{accelY }{Y component of the accel of the satellite in km/s, in the
  ECEF system of coordinates}
  \item{accelZ }{Z component of the accel of the satellite in km/s, in the
  ECEF system of coordinates}
  \item{satelliteHealthCode }{Code indicating the health of the satellite. 0 if
  healthy}
  \item{freqNumber }{Frequency number (k) of the GLONASS satellite. The two 
  frequencies in MHz, f1 and f2, used by the satellite to transmit data can be 
  calculated as follows: f1 = 1602 + k*9/16 and f2 = 1246 + k*7/16}
  \item{informationAge }{Age in days of the observation data used to generate
  the provided ephemeris}
}
\references{
https://gage.upc.edu/gFD/
https://www.navcen.uscg.gov/pubs/gps/rinex/rinex.txt
ftp://www.ngs.noaa.gov/cors/RINEX211.txt
http://acc.igs.org/misc/rinex304.pdf
https://russianspacesystems.ru/wp-content/uploads/2016/08/ICD_GLONASS_eng_v5.1.pdf
}
\examples{
# The file testGLONASSRINEX.txt provided with the package includes 5 navigation
# messages from 4 GLONASS satellites

testGLONASSnav <- readGLONASSNavigationRINEX(paste0(path.package("asteRisk"), 
"/testGLONASSRINEX.txt"))
testGLONASSnav$header
testGLONASSnav$messages
}
