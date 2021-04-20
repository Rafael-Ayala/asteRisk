calculatePolarMotionMatrix <- function(julianDate) {
    MJD <- julianDate - 2400000.5
    A <- 2 * pi * (MJD - 57226) / 365.25
    C <- 2 * pi * (MJD - 57226) / 435
    xp <- (0.1033 + 0.0494*cos(A) + 0.0482*sin(A) + 0.0297*cos(C) + 0.0307*sin(C)) * 4.84813681e-6
    yp <- (0.3498 + 0.0441*cos(A) - 0.0393*sin(A) + 0.0307*cos(C) - 0.0297*sin(C)) * 4.84813681e-6
    polarMotionMatrix <- matrix(c(cos(xp), 0, -sin(xp),
                                  sin(xp)*sin(yp), cos(yp), cos(xp)*sin(yp),
                                  sin(xp)*cos(yp), -sin(yp), cos(xp)*cos(yp)),
                                nrow=3, ncol=3, byrow=TRUE)
    return(polarMotionMatrix)
}

TEMEtoECEF <- function(position_TEME, velocity_TEME=c(0,0,0), dateTime) {
    gmst <- UTCdateTimeToGMST(dateTime)
    daysToJ2000_0 <- as.numeric(julian(as.POSIXct(dateTime, tz="UTC"),
                                       origin=as.POSIXct("2000-01-01 12:00:00", tz="UTC")))
    julianDate <- daysToJ2000_0 + JD_J2000_0
    PEF_TOD_matrix <- matrix(c(cos(gmst), -sin(gmst), 0,
                               sin(gmst), cos(gmst), 0,
                               0, 0, 1),
                             nrow=3, ncol=3, byrow=TRUE)
    position_PEF <- as.vector(t(PEF_TOD_matrix) %*% position_TEME)
    polarMotionMatrix <- calculatePolarMotionMatrix(julianDate)
    position_ECEF <- as.vector(t(polarMotionMatrix) %*% position_PEF)
    omegaEarth <- c(0, 0, 7.29211514670698e-05 * (1.0  - 0.002/86400.0))
    velocity_PEF <- as.vector(t(PEF_TOD_matrix) %*% velocity_TEME) -
        c(omegaEarth[2] * position_PEF[3] - omegaEarth[3] * position_PEF[2],
          omegaEarth[3] * position_PEF[1] - omegaEarth[1] * position_PEF[3],
          omegaEarth[1] * position_PEF[2] - omegaEarth[2] * position_PEF[1])
    velocity_ECEF <- as.vector(t(polarMotionMatrix) %*% velocity_PEF)
    return(list(
        position=position_ECEF,
        velocity=velocity_ECEF
    ))
}

ECEFtoLATLON <- function(position_ECEF, degreesOutput=TRUE) {
    a <- 6378137.0
    a2 <- a^2
    f <- 1/298.257223563
    b <- a*(1-f)
    b2 <- b^2
    e <- sqrt((a2 - b2)/a2)
    eprime <- sqrt((a2 - b2)/b2)
    X <- position_ECEF[1]
    Y <- position_ECEF[2]
    Z <- position_ECEF[3]
    p <- sqrt(X^2 + Y^2)
    theta <- atan2(a*Z, b*p)
    sintheta <- sin(theta)
    costheta <- cos(theta)
    num <- Z + eprime^2 * b * sintheta^3
    denom <- p - e^2 * a * costheta^3
    lat <- atan2(num, denom)
    lon <- atan2(Y, X)
    N <- a/(sqrt(1 - e^2 * (sin(lat))^2))
    alt <- p / cos(lat) - N
    if(abs(lon) >= pi) {
        if(lon < 0) {
            lon <- lon + 2*pi
        } else {
            lon <- lon - 2*pi
        }
    }
    LATLONALT <- if(degreesOutput)  c(rad2deg(lat), rad2deg(lon), alt) else c(lat, lon, alt)
    names(LATLONALT) <- c("latitude", "longitude", "altitude")
    return(LATLONALT)
}

TEMEtoLATLON <- function(position_TEME, dateTime) {
    ECEFcoords <- TEMEtoECEF(position_TEME=position_TEME, dateTime=dateTime)
    return(ECEFtoLATLON(ECEFcoords$position))
}

ECEFtoICRF <- function(position_ECEF, velocity_ECEF=c(0, 0, 0), dateTime) {
    date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
    year <- date$year + 1900
    month <- date$mon + 1
    day <- date$mday
    hour <- date$hour
    minute <- date$min
    second <- date$sec
    Mjd_UTC <- MJday(year, month, day, hour, minute, second)
    results <- ECEFtoECI(Mjd_UTC, c(position_ECEF, velocity_ECEF))
    return(results)
}

ICRFtoECEF <- function(position_ICRF, velocity_ICRF=c(0, 0, 0), dateTime) {
    date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
    year <- date$year + 1900
    month <- date$mon + 1
    day <- date$mday
    hour <- date$hour
    minute <- date$min
    second <- date$sec
    Mjd_UTC <- MJday(year, month, day, hour, minute, second)
    results <- ECItoECEF(Mjd_UTC, c(position_ICRF, velocity_ICRF))
    return(results)
}
