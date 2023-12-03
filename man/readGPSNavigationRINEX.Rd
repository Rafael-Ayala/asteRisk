\name{readGPSNavigationRINEX}
\alias{readGPSNavigationRINEX}
\title{Read a RINEX navigation file for GPS satellites
}
\description{
RINEX (Receiver Independent Exchange Format) is one of the most widely used
formats for providing data of satellite navigation systems. The RINEX standard
defines several structured text file types file types, among which navigation 
files are used to distribute positional information of the satellites. The exact 
information provided in a RINEX navigation file varies for each satellite navigation 
system. This function reads RINEX navigation files for satellites of the GPS 
constellation, operated by the USA.
}
\usage{
readGPSNavigationRINEX(filename)
}
\arguments{
  \item{filename}{Path to the GPS RINEX navigation file.
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
  \item{ionAlphaA0 }{Coefficient for ionospheric correction A0}
  \item{ionAlphaA1 }{Coefficient for ionospheric correction A1}
  \item{ionAlphaA2 }{Coefficient for ionospheric correction A2}
  \item{ionAlphaA3 }{Coefficient for ionospheric correction A3}
  \item{ionBetaB0 }{Coefficient for ionospheric correction B0}
  \item{ionBetaB1 }{Coefficient for ionospheric correction B1}
  \item{ionBetaB2 }{Coefficient for ionospheric correction B2}
  \item{ionBetaB3 }{Coefficient for ionospheric correction B3}
  \item{deltaUTCA0 }{A0 parameter, corresponding to bias between GPST and UTC
  time at the reference time (Tot) given by fields \code{referenceTimeUTC} and
  \code{referenceWeekUTC}. Should be used to compute accurate time in UTC}
  \item{deltaUTCA1 }{A1 parameter, corresponding to the clock drift between
  GPST and UTC at the reference time (Tot) given by fields \code{referenceTimeUTC} 
  and \code{referenceWeekUTC}. Should be used to compute accurate time in UTC}
  \item{referenceTimeUTC }{Time in seconds of current UTC week of Tot, which is
  the reference time to correct GPST time to UTC}
  \item{referenceWeekUTC }{UTC reference week number (continuous scale, not
  modulo 1024) of Tot.}
  \item{leapSeconds }{Leap seconds introduced since the 6th of January, 1980. 
  Useful to convert to UTC time (UTC time = GPS time - leap seconds)}
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
  IRGP, for IRN to GPS (A0/A1). Note that GPS RINEX navigation files will
  typically contain system time corrections of type GPUT}
  \item{timeCorrectionA0 }{A0 parameter (bias) for system time correction in seconds}
  \item{timeCorrectionA1 }{A1 parameter (drift) for system time correction in seconds/seconds}
  \item{timeCorrectionReferenceTime }{Reference time for system time corrections,
  in seconds into GPS/GAL/BDS/QZS/IRN/SBAS week. The correction to be applied to 
  the system time should be calculated as:
  Correction=A0 + A1*(ephemerisTime - timeCorrectionReferenceTime) }
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
of the GPS satellite:
  \item{satellitePRNCode }{PRN code of the satellite. Unique PRN codes are
  assigned to all satellites in global navigation satellite systems, and 
  therefore provide an identifier for each of them}
  \item{tocYear }{Toc year in 4-digit format. Toc is the GPS time of the specific 
  satellite that should be used as the time reference to apply clock bias, clock 
  drift and possibly even clock drift rate, as well as a relativistic correction, 
  as described in the GPS system specification (https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf) 
  to obtain the corrected GPST system time. The GPST system time can be converted 
  to UTC time by subtracting leap seconds since the 6th of January 1980 and performing
  another polynomial correction to account for bias and drift between GPST and
  UTC times.}
  \item{tocMonth }{Toc month}
  \item{tocDay }{Toc day}
  \item{tocHour }{Toc hour}
  \item{tocMinute }{Toc minute}
  \item{tocSecond }{Toc second}
  \item{clockBias }{Clock bias (i.e., constant offset) at Toc that should be 
  applied to the satellite time in order to calculate accurate GPST. In seconds.
  Often referred to as af0.}
  \item{clockDrift }{Clock drift of the satellite clock at Toc that should be 
  applied to the satellite time in order to calculate accurate GPST. In seconds.
  Often referred to as af1.}
  \item{clockDriftRate }{Rate of change for the clock drift of the satellite 
  clock at Toc. It is frequently 0, but if not, it should be applied in 
  combination with clock bias and clock drift in order to correct to GPST as
  accurately as possible. In seconds per square second. Often referred to as af2.}  
  \item{IODE }{Issue-of-data ephemeris. It acts as a time-stamp or unique 
  identifier for the provided navigation data. In particular, the IODE of a 
  given navigation message should never be the same as the IODE for any other
  navigation message broadcasted by the same satellite in the past 6 days,
  although violations of this rule have been observed. Most frequently, IODE
  are not reused in a period of 7 days, so that they match exactly the IODC.}
  \item{radiusCorrectionSine }{Amplitude of the sine harmonic component for the 
  correction of orbital radius. In meters}
  \item{deltaN }{Mean motion difference from computed value. In radians per
  second. In order to obtain the real (perturbed) mean motion, first the 
  Keplerian mean motion should be calculated from the semi-major axis. Then, 
  deltaN should be added to it.}
  \item{correctedMeanMotion }{Corrected mean motion calculated by adding deltaN
  to the value computed from the semi-major axis. In radians per second}
  \item{meanAnomaly }{Mean anomaly of the satellite at epoch. In radians. This 
  indicates where the satellite is along its orbital path. It is provided as the 
  angle between the direction of the perigee and the hypothetical point where 
  the object would be if it was moving in a circular orbit with the same period 
  as its true orbit after the same amount of time since it last crossed the 
  perigee had ellapsed. Therefore, 0 denotes that the object is at the perigee. 
  This is a Keplerian orbital element.}
  \item{latitudeCorrectionCosine }{Amplitude of the cosine harmonic component
  for the correction of latitude argument. In radians}
  \item{eccentricity }{Eccentricity of the orbit of the satellite at epoch. This
  is a Keplerian orbital element.}
  \item{latitudeCorrectionSine }{Amplitude of the sine harmonic component for 
  the correction of latitude argument. In radians}
  \item{semiMajorAxis }{Semi-major axis of the orbit of the satellite at epoch.
  In meters. This is a Keplerian orbital element}
  \item{toeSecondsOfGPSWeek }{Time of the GPS week (in seconds) for the ephemeris.
  Together with the \code{toeGPSWeek}, it can be used to calculate the ephemeris
  time in GPS time of the specific satellite, to which sequential corrections
  first to GPST and then to UTC should be applied.}
  \item{inclinationCorrectionCosine }{Amplitude of the cosine harmonic component
  for the correction of inclination. In radians}
  \item{ascension }{Longitude of the ascending node of the satellite at epoch, 
  also known as right ascension of the ascending node, in radians. This is the 
  angle between the direction of the ascending node (the point where the 
  satellite crosses the equatorial plane moving north) and the direction of the 
  First Point of Aries (which indicates the location of the vernal equinox).
  This is a Keplerian orbital element.}
  \item{inclinationCorrectionSine }{Amplitude of the sine harmonic component
  for the correction of inclination. In radians}
  \item{inclination }{Mean orbital inclination of the satellite in radians. This
  is the angle between the orbital plane of the satellite and the equatorial 
  plane. This is a Keplerian orbital element.}
  \item{radiusCorrectionCosine }{Amplitude of the cosine harmonic component for 
  the correction of orbital radius. In meters.}
  \item{perigeeArgument }{Mean argument of the perigee of the object in radians. 
  This is the angle between the direction of the ascending node and the 
  direction of the perigee (the point of the orbit at which the object is 
  closest to the Earth). This is a Keplerian orbital element.}
  \item{OMEGADot }{Angular velocity of the satellite with respect to the vernal
  equinox. In radians/second.}
  \item{codesL2Channel }{Flag indicating if coarse/acquisition (C/A) code is
  being transmitted on the L2 channel (value of 1) or not (value of 0)}
  \item{toeGPSWeek }{GPS week number at epoch}
  \item{dataFlagL2P }{Flag indicating if precise (P) code is being transmitted 
  on the L2 channel (value of 1) or not (value of 0)}
  \item{satelliteAccuracy }{Accuracy of the position of the satellite, in 
  meters.}
  \item{satelliteHealthCode }{Code indicating the health of the satellite. 0 if
  healthy.}
  \item{totalGroupDelay }{Bias difference between codes broadcasted on L1 and
  the ionospheric-free combination of the codes broadcasted at L1 and L2, in
  seconds. This parameter, also known as timing group delay (TGD), should be 
  considered when calculating satellite clock error.}
  \item{IODC }{Issue-of-data clock. It acts as a time-stamp or unique identifier 
  for the provided navigation data. In particular, the IODC of a given 
  navigation message should never be the same as the IODC for any other
  navigation message broadcasted by the same satellite in the past 7 days,
  although violations of this rule have been observed. Most frequently, IODE
  are not reused in a period of 7 days instead of the mandatory 6 days, so that 
  they match exactly the IODC.}
  \item{transmissionTime }{Transmission time for the navigation message, in
  seconds of GPS week.}
  \item{fitInterval }{Flag indicating for how long the broadcasted ephemeris are
  valid since the last time the data was updated. It should be noted that in
  order to obtain positional values/orbital elements at times other than epoch,
  the corrections for perturbed orbital elements should be applied and 
  propagated. If 0, the ephemeris data are valid for up to 4 hours. If 1, they
  are valid for more than 4 hours.}
  \item{ephemerisUTCTime }{A \link{nanotime} object indicating the time
  corresponding to the reported position (ephemeris) in the present message. 
  The time is in UTC, obtained by first applying the individual clock bias, clock 
  drift and clock drift rate of the particular satellite (fields \code{clockBias},
  \code{clockDrift} and \code{clockDriftRate} of each message) and a relativistic
  correction to obtain corrected GPST time (system-wide GPS time), and then
  subtraction of leap seconds since the 6th of January 1980 and a second polynomial
  correction (with fields \code{deltaUTCA0} and \code{deltaUTCA1} from the 
  header) to obtain UTC time. Corrections are performed as described in the GPS
  system specifications (https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf).}
  \item{position_ITRF}{Position of the satellite in the ITRF frame, calculated
  from the provided orbital ephemeris following the procedure described in the
  GPS system specifications. In meters.}
  \item{velocity_ITRF}{Velocity of the satellite in the ITRF frame, calculated
  from the provided orbital ephemeris following the procedure described in the
  GPS system specifications. In meters/second.}
  \item{acceleration_ITRF}{Acceleration of the satellite in the ITRF frame, calculated
  from the provided orbital ephemeris following the procedure described in the
  GPS system specifications. In meters/squared second.}
}
\references{
https://gage.upc.edu/gFD/
https://www.navcen.uscg.gov/pubs/gps/rinex/rinex.txt
ftp://www.ngs.noaa.gov/cors/RINEX211.txt
http://acc.igs.org/misc/rinex304.pdf
https://www.icao.int/Meetings/AMC/MA/2004/GNSS/icd.pdf
https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf
}
\examples{
# The file testGPSRINEXv2.txt provided with the package includes 3 navigation
# messages from 3 GPS satellites

testGPSnav <- readGPSNavigationRINEX(paste0(path.package("asteRisk"), 
"/testGPSRINEXv2.txt"))
testGPSnav$header
testGPSnav$messages
}
