\name{readGLONASSNavigationRINEX}
\alias{readGLONASSNavigationRINEX}
\title{Read a RINEX navigation file for GLONASS satellites
}
\description{
RINEX (Receiver Independent Exchange Format) is one of the most widely used
formats for providing data of satellite navigation systems. The RINEX standard
defines several structured text file types, among which navigation files are used 
to distribute positional information of the satellites. The exact information 
provided in a RINEX navigation file varies for each satellite navigation system. 
This function reads RINEX navigation files for satellites of the GLONASS 
constellation, operated by Russia.
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
the information contained in the header of the RINEX file. For files using RINEX
major version 2, it contains the following elements:
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
  very small amount. This is the negative of the parameter typically referred to
  as tauC. }
  \item{leapSeconds }{Leap seconds introduced since 1980. Useful to convert to
  GPS time}
  \item{comments }{Miscellaneous comments found in the header of the RINEX file}

For files using RINEX major version 3, it contains the following elements:
  \item{rinexVersion }{Version of the RINEX format used in the file}
  \item{rinexFileType }{Type of RINEX file}
  \item{satelliteSystem }{Character indicating the satellite system. For GLONASS,
  it should be "R"}
  \item{generatorProgram }{Program used to generate the RINEX file}
  \item{generatorEntity }{Individual or organization that generated the file}
  \item{fileCreationDateString }{Date-time string indicating when the file was
  created}
  \item{systemTimeCorrectionType }{String indicating the type of system time
  correction, which defines the exact meaning of the system time correction
  parameters A0 and A1 and between what time systems these allow conversion. 
  Possible values are (the definition of A0/A1 are given between brackets): 
  GAUT, for GAL to UTC (a0/a1); GPUT, for GPS to UTC (a0/a1); SBUT, for SBAS to UTC 
  (a0/a1); GLUT, for GLO to UTC (TauC/0); GAGP or GPGA, for GPS to GAL (A0G/A1G);
  GLGP, for GLO to GPS (TauGPS/0); QZGP, for QZS to GPS (a0/a1); QZUT, for QZS to
  UTC (a0/a1); BDUT, for BDS to UTC (A0UTC/A1UTC); IRUT, for IRN to UTC (A0UTC/A1UTC);
  IRGP, for IRN to GPS (A0/A1). Note that GLONASS RINEX navigation files will
  typically contain system time corrections of types either GLUT or GLGP}
  \item{timeCorrectionA0 }{A0 parameter (bias) for system time correction in seconds.
  For corrections of type GLUT, this is equivalent to the field \code{sysTimeCorrection}
  found in headers of GLONASS RINEX navigation files of major version 2, i.e.,
  the negative of tauC. For corrections of type GLGP, this is equivalent to the
  parameter tauGPS}
  \item{timeCorrectionA1 }{A1 parameter (drift) for system time correction in seconds/seconds.
  Note that this is set to 0 in GLUT and GLGP correction types}
  \item{timeCorrectionReferenceTime }{Reference time for system time corrections,
  in seconds into GPS/GAL/BDS/QZS/IRN/SBAS week. The correction to be applied to 
  the system time should be calculated as:
  Correction=A0 + A1*(ephemerisTime - timeCorrectionReferenceTime).
  Note that GLONASS time is aligned with UTC + 3 hours at a precision of 1 ms, and
  timeCorrectionReferenceTime is set to 0.}
  \item{timeCorrectionReferenceWeekNumber }{Reference week number for the system 
  time correction reference time. For GPS/GAL/QZS/IRN/SBAS, it is a continuous
  week scale since 6th of January 1980. For BDS, it is a continuous week scale 
  since 1st of January 2006. For GLONASS, it is set to 0}
  \item{timeCorrectionSatelliteNumber }{String indicating the GNSS satellite that
  is broadcasting the time system difference for system time correction }
  \item{UTCType }{String indicating the exact UTC type (specific UTC realization) 
  to which system time correction parameters refers to }
  \item{leapSeconds }{Current number of leap seconds }
  \item{deltaTimeLeapSeconds }{Number of leap seconds at a reference date for leap
  seconds specified by the following week number and week day. Note that the date 
  specified by said week number and day can be either in the past or future }
  \item{deltaTimeLeapSecondsWeekNumber }{Week number for the reference date to
  which deltaTimeLeapSeconds refers. Given as continuous week number from either
  6th of January 1980 (GPS week number) or from 1st of January 2006 (BDS week
  number). Field leapSecondsTimeSystemIdentifier specifies which of the 2 systems
  is used.}
  \item{deltaTimeLeapSecondsDayNumber }{Day of week number for the reference date 
  to which deltaTimeLeapSeconds refers. Ranges from 1-7 for reference dates in GPS,
  or 0-6 for reference dates in BDS. The first day of each week (1 for GPS, 0 for
  BDS) is considered to be Sunday }
  \item{leapSecondsTimeSystemIdentifier }{String indicating the time system used
  for the week number and day number for the future/past leap seconds reference 
  date. Can be either "GPS" or "BDS" (BeiDou)}
  \item{leapSecondsTimeSystemIdentifier }{String indicating the time system used
  for the week number and day number for the future/past leap seconds reference 
  date. Can be either "GPS" or "BDS" (BeiDou)}
  \item{ionosphericCorrections }{List where each element contains the fields required
  to perform a given type of ionospheric corrections. The contents of each element
  of this list are detailed below }
  \item{comments }{Miscellaneous comments found in the header of the RINEX file}
  
As mentioned above, for RINEX files of version 3.00 and above, the header element
contains a list named \code{ionosphericCorrections}. Each element of this list
is in turn a list itself, with the following elements:
  \item{ionosphericCorrectionType }{A string indicating the type of ionospheric
  correction for which this element provides coefficients. Can be GAL, GPSA, GPSB,
  QZSA, QZSB, BDSA, BDSB, IRNA or IRNB. The specific coefficients given in the
  following elements vary depending on the type of ionospheric correction. For GAL,
  3 coefficients are given (a0-a2). For GPSA, QZSA, BDSA and IRNA, 4 coefficients
  are given (alpha0-alpha3). For GPSB, QZSB, BDSB and IRNB, 4 coefficients are
  given (beta0-beta1)}
  \item{coefficients }{3 (for GAL correction) or 4 (for all other corrections)
  elements providing the ionospheric correction coefficients. The names of the
  elements vary depending on the type of correction, as stated above}
  \item{ionosphericCorrectionTimeMark }{Character indicating the transmission time.
  The field is only mandatory for corrections of type BDSA and BDSB. It can be a letter
  from "A" to "X", with each letter indicating a 1-hour interval: A=00h-01h, B=01h-02hB,
  ..., X=23h-24h }
  \item{ionosphericCorrectionTimeMark }{Character indicating the transmission time.
  The field is only mandatory for corrections of type BDSA and BDSB. It can be a letter
  from "A" to "X", with each letter indicating a 1-hour interval: A=00h-01h, B=01h-02hB,
  ..., X=23h-24h }
  \item{ionosphericCorrectionSV }{String indicating the satellite that provided
  the ionospheric correction parameters. The field is only mandatory for BDSA and
  BDSB corrections. If multiple sources are available for the same type of BDSA/BDSB
  correction, these should be given priority according to the satellite that provided
  them as follows: BDS GEO satellites first, followed by BDS IGSO satellites,
  and finally BDS MEO satellites}

The second element is named \code{messages}, and it contains one element for
each navigation message found in the RINEX file. Each of these elements is a 
list with the following elements that provide information about the position
of the GLONASS satellite:
  \item{satelliteNumber }{Slot number of the satellite within the GLONASS
  constellation. It can be converted to a PRN code by adding 37 to it}
  \item{epochYear }{Epoch year in 4-digit format.}
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
  specifications (http://gauss.gge.unb.ca/GLONASS.ICD.pdf)}
  \item{clockBias }{Clock bias (i.e., constant offset) that should be applied to
  the satellite time in order to obtain an even more accurate UTC time. In 
  seconds}
  \item{relativeFreqBias }{Clock drift of the satellite clock that should be
  applied in combination with the time difference to the reference time in order
  to obtain an even more accurate UTC time. In seconds per second}
  \item{messageFrameTime }{Second of the UTC day when the message was 
  transmitted}
  \item{positionX }{X coordinate of the position of the satellite in km, in the
  ITRF system of coordinates}
  \item{positionY }{Y coordinate of the position of the satellite in km, in the
  ITRF system of coordinates}
  \item{positionZ }{Z coordinate of the position of the satellite in km, in the
  ITRF system of coordinates}
  \item{velocityX }{X component of the velocity of the satellite in km/s, in the
  ITRF system of coordinates}
  \item{velocityY }{Y component of the velocity of the satellite in km/s, in the
  ITRF system of coordinates}
  \item{velocityZ }{Z component of the velocity of the satellite in km/s, in the
  ITRF system of coordinates}
  \item{accelX }{X component of the accel of the satellite in km/s, in the
  ITRF system of coordinates}
  \item{accelY }{Y component of the accel of the satellite in km/s, in the
  ITRF system of coordinates}
  \item{accelZ }{Z component of the accel of the satellite in km/s, in the
  ITRF system of coordinates}
  \item{satelliteHealthCode }{Code indicating the health of the satellite. 0 if
  healthy}
  \item{freqNumber }{Frequency number (k) of the GLONASS satellite. The two 
  frequencies in MHz, f1 and f2, used by the satellite to transmit data can be 
  calculated as follows: f1 = 1602 + k*9/16 and f2 = 1246 + k*7/16}
  \item{informationAge }{Age in days of the observation data used to generate
  the provided ephemeris}

For GLONASS RINEX navigation files of version 3.05 and above, each message element
contains the following additional elements:
  \item{GLONASSType }{String indicating the type of GLONASS satellite. Can be either 
  "GLO" for fist-generation GLONASS satellites, or "GLO-M/K" for second or third
  generation satellites}
  \item{updatedDataFlag }{Logical indicating if the provided ephemeris data is
  up to date}
  \item{numberSatellitesAlmanac }{Number of satellites in the almanac for the current
  transmitted frame. Can be either 4 or 5. Note that almanac data provide coarser
  information about the location of multiple satellites in a GNSS constellation, and
  these are not actually included in RINEX navigation files.}
  \item{ephemerisValidityTimeInterval }{Length of the time interval for which the
  ephemeris is valid, in minutes. Can be 0, 30, 45 or 60.}
  \item{parityEphemerisValidityTimeInterval }{String indicating the evenness or 
  oddity of the time interval for which ephemeris is valid. "Odd" for intervals of
  45 minutes, and "Even" for the rest.}
  \item{tauCSource }{String indicating the source providing the value of the tauC
  parameter, given in the header in field \code{sysTimeCorrection} for RINEX navigation
  messages of version 2, or field \code{timeCorrectionA0} for RINEX navigation
  messages of version 3 where a time system correction of type GLUT is given.
  Can be either "On-board" (if it was computed by the on-board satellite processor)
  or "Ground" (computed and uploaded to the satellite by control segment)}
  \item{tauGPSSource }{String indicating the source providing the value of the tauGPS
  parameter, given in the header in field \code{timeCorrectionA0} for RINEX navigation
  messages of version 3 where a time system correction of type GLGP is given.
  Can be either "On-board" (if it was computed by the on-board satellite processor)
  or "Ground" (computed and uploaded to the satellite by control segment)}
  \item{totalGroupDelay }{Bias difference between codes broadcasted on L1 and
  the ionospheric-free combination of the codes broadcasted at L1 and L2, in
  seconds. This parameter, also known as timing group delay (TGD), should be 
  considered when calculating satellite clock error.}
  \item{URAI }{Value of User Range Accuracy Index (URAI). This is an index giving 
  a measurement of the accuracy of the GNSS ranging accuracy. It is an integer
  ranging from 0 to 15. 0 indicates the highest accuracy, corresponding to an
  accuracy of 1 m, while 14 indicates the lowest accuracy, of 512 m. A value
  of 15 indicates unknown accuracy. For a complete table of equivalence between
  URAI values and ranging accuracy in meters, see Table 4.4 of GLONASS Interface
  Control Document (http://gauss.gge.unb.ca/GLONASS.ICD.pdf)}
  \item{almanacHealthStatus }{String indicating the health status of the satellite
  provided in the almanac (as previously mentioned, full almanac data are not
  present in RINEX navigation files))}

}
\references{
https://gage.upc.edu/gFD/
https://www.navcen.uscg.gov/pubs/gps/rinex/rinex.txt
ftp://www.ngs.noaa.gov/cors/RINEX211.txt
http://acc.igs.org/misc/rinex304.pdf
http://gauss.gge.unb.ca/GLONASS.ICD.pdf
}
\examples{
# The file testGLONASSRINEXv2.txt provided with the package includes 5 navigation
# messages from 4 GLONASS satellites

testGLONASSnav <- readGLONASSNavigationRINEX(paste0(path.package("asteRisk"), 
"/testGLONASSRINEXv2.txt"))
testGLONASSnav$header
testGLONASSnav$messages
}
