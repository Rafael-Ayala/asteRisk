deg2rad <- function(degrees) {
    return(degrees*pi/180)
}

rad2deg <- function(radians) {
    return(radians*180/pi)
}

UTCdateTimeToGMST <- function(dateTime) {
    daysToJ2000_0 <- as.numeric(julian(as.POSIXct(dateTime, tz="UTC"),
                                       origin=as.POSIXct("2000-01-01 12:00:00", tz="UTC")))
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
    a <- (earth_mu)^(1/3) / ((2*pi*meanMotion/86400)^(2/3))
    return(a)
}

semiMajorAxisToMeanMotion <- function(semiMajorAxis) {
    n <- sqrt(earth_mu/semiMajorAxis^3) * (86400/(2*pi))
    return(n)
}

vectorCrossProduct3D <- function(u, v) {
    w <- c(u[2]*v[3] - u[3]*v[2],
           u[3]*v[1] - u[1]*v[3],
           u[1]*v[2] - u[2]*v[1])
    return(w)
}

getLatestSpaceData <- function() { # TODO: MOVE TO THIS PACKAGE
    hasData()
    asteRiskData::getLatestSpaceData()
}
