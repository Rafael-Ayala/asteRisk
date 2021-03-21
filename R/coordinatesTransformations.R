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

TEMEtoECEF <- function(position_TEME, velocity_TEME, dateTime) {
    gmst <- UTCdateTimeToGMST(dateTime)
    daysToJ2000_0 <- as.numeric(julian(as.POSIXct(dateTime, tz="UTC"),
                                       origin=as.POSIXct("2000-01-01 12:00:00", tz="UTC")))
    julianDate <- daysToJ2000_0 + JD_J2000_0
    PEF_TOD_matrix <- matrix(c(cos(gmst), -sin(gmst), 0,
                               sin(gmst), cos(gmst), 0,
                               0, 0, 1),
                             nrow=3, ncol=3, byrow=TRUE)
    position_PEF <- as.vector(PEF_TOD_matrix %*% position_TEME)
    polarMotionMatrix <- calculatePolarMotionMatrix(julianDate)
    position_ECEF <- as.vector(polarMotionMatrix %*% position_PEF)
    omegaEarth <- c(0, 0, 7.29211514670698e-05 * (1.0  - 0.002/86400.0))
    velocity_PEF <- as.vector(PEF_TOD_matrix %*% velocity_TEME) -
        c(omegaEarth[2] * position_PEF[3] - omegaEarth[3] * position_PEF[2],
          omegaEarth[3] * position_PEF[1] - omegaEarth[1] * position_PEF[3],
          omegaEarth[1] * position_PEF[2] - omegaEarth[2] * position_PEF[1])
    velocity_ECEF <- as.vector(polarMotionMatrix %*% velocity_PEF)
    return(list(
        position=position_ECEF,
        velocity=velocity_ECEF
    ))
}


