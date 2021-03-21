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

ECEFtoLATLON <- function(position_ECEF) {
    # based on the 2015 C++ port by Grady Hillhouse of the 
    # 2007 Matlab implementation by David Vallado
    tolerance <- 1e-8
    position_mag <- sqrt(sum(position_ECEF^2))
    tmp <- sqrt(position_ECEF[1]^2 + position_ECEF[2]^2)
    if (abs(tmp) < tolerance) {
        rtasc <- sign(position_ECEF[3]) * pi * 0.5
    } else {
        rtasc <- atan2(position_ECEF[2], position_ECEF[1])
    }
    lon <- rtasc
    if(abs(lon) >= pi) {
        if(lon < 0) {
            lon <- lon + 2*pi
        } else {
            lon <- lon - 2*pi
        }
    }
    lat <- asin(position_ECEF[3]/position_mag)
    i <- 1
    olddelta <- lat + 10
    sintemp <- 0
    c <- 0
    while ((abs(olddelta - lat) >=  tolerance) & (i < 10)) {
        olddelta <- lat
        sintemp <- sin(lat)
        c <- XKMPER / sqrt(1 - earthEccentricity^2 * sintemp^2)
        lat <- atan( (position_ECEF[3] + c*earthEccentricity^2*sintemp) / tmp )
        i <- i+1
    }
    if ((0.5 * pi - abs(lat)) > pi/180) {
        alt <- (tmp/cos(lat)) - c
    } else {
        alt <- position_ECEF[3]/sin(lat) - c * (1 - earthEccentricity^2)
    }
    LATLONALT <- c(rad2deg(lat), rad2deg(lon), alt)
    names(LATLONALT) <- c("latitude", "longitude", "altitude")
    return(LATLONALT)
}

TEMEtoLATLON <- function(position_TEME, dateTime) {
    ECEFcoords <- TEMEtoECEF(position_TEME=position_TEME, dateTime=dateTime)
    return(ECEFtoLATLON(ECEFcoords$position))
}
