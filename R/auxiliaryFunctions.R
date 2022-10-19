deg2rad <- function(degrees) {
    return(degrees*pi/180)
}

rad2deg <- function(radians) {
    return(radians*180/pi)
}

revDay2radMin <- function(revPerDay) {
    return(revPerDay*((2*pi)/(1440)))
}

UTCdateTimeToGMST <- function(dateTime, convertUTCtoUT1 = FALSE) {
    # Formula is actually for converting from UT1 Julian date, but difference 
    # between UTC and UT1 JD will be below 0.9 seconds
    if(convertUTCtoUT1) {
        hasData()
        date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
        year <- date$year + 1900
        month <- date$mon + 1
        day <- date$mday
        hour <- date$hour
        minute <- date$min
        second <- date$sec
        MJD_UTC <- iauCal2jd(year, month, day, hour, minute, second)$DATE
        IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "n")
        deltaUT1_UTC <- IERS_results$UT1_UTC
    } else {
        deltaUT1_UTC <- 0
    }
    daysToJ2000_0 <- as.numeric(julian(as.POSIXct(dateTime, tz="UTC") + deltaUT1_UTC,
                                       origin=as.POSIXct("2000-01-01 12:00:00", tz="UTC")))
    centuriesFromJ2000 <- daysToJ2000_0/36525
    GMST <- 67310.54841 + (876600.0*3600 + 8640184.812866)*centuriesFromJ2000+ 0.093104*centuriesFromJ2000^2 - 6.2e-6*centuriesFromJ2000^3
    GMST <- GMST %% 86400
    return(GMST*pi/43200)
}

MJDToGMST <- function(MJD, convertUTCtoUT1 = FALSE) {
    if(convertUTCtoUT1) {
        MJD <- MJDUTCtoMJDUT1(MJD)
    }
    daysToJ2000_0 <- MJD - MJD_J2000
    centuriesFromJ2000 <- daysToJ2000_0/36525
    GMST <- 67310.54841 + (876600.0*3600 + 8640184.812866)*centuriesFromJ2000+ 0.093104*centuriesFromJ2000^2 - 6.2e-6*centuriesFromJ2000^3
    GMST <- GMST %% 86400
    return(GMST*pi/43200)
}

rem <- function(x, y) {
    remainder <- x - trunc(x/y)*y
    return(remainder)
}

meanMotionToSemiMajorAxis <- function(meanMotion) {
    # output will be in meters
    a <- (GM_Earth_TDB)^(1/3) / ((2*pi*meanMotion/86400)^(2/3))
    return(a)
}

semiMajorAxisToMeanMotion <- function(semiMajorAxis, outputRevsPerDay=TRUE) {
    # provide input in meters
    n <- sqrt(GM_Earth_TDB/semiMajorAxis^3)
    if(outputRevsPerDay) n <- n * (86400/(2*pi))
    return(n)
}

vectorCrossProduct3D <- function(u, v) {
    w <- c(u[2]*v[3] - u[3]*v[2],
           u[3]*v[1] - u[1]*v[3],
           u[1]*v[2] - u[2]*v[1])
    return(w)
}

getLatestSpaceData <- function(targets="all") { # TODO: MOVE TO THIS PACKAGE?
    hasData()
    asteRiskData::getLatestSpaceData_(targets = targets)
}
