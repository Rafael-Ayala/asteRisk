calculatePolarMotionMatrix <- function(julianDate) {
    hasData()
    # Old code where polar coefficients xp and yp are calculated instead of 
    # retrieved from data, based on Grady Hillhouse C++ implementation of
    # Vallado's MatLab implementation
    # MJD <- julianDate - 2400000.5
    # A <- 2 * pi * (MJD - 57226) / 365.25
    # C <- 2 * pi * (MJD - 57226) / 435
    # xp <- (0.1033 + 0.0494*cos(A) + 0.0482*sin(A) + 0.0297*cos(C) + 0.0307*sin(C)) * 4.84813681e-6
    # yp <- (0.3498 + 0.0441*cos(A) - 0.0393*sin(A) + 0.0307*cos(C) - 0.0297*sin(C)) * 4.84813681e-6
    # polarMotionMatrix <- matrix(c(cos(xp), 0, -sin(xp),
    #                               sin(xp)*sin(yp), cos(yp), cos(xp)*sin(yp),
    #                               sin(xp)*cos(yp), -sin(yp), cos(xp)*cos(yp)),
    #                             nrow=3, ncol=3, byrow=TRUE)
    # New code using observed/predicted values retrieved from Celestrak, still
    # using 80's nutation theory
    MJD <- trunc(julianDate - 2400000.5)
    earthPositionsRow <- asteRiskData::earthPositions[asteRiskData::earthPositions[,4] == MJD, ]
    # Multiplication factor to convert from arcseconds to radians
    xp <- earthPositionsRow[5] * 4.84813681e-6
    yp <- earthPositionsRow[6] * 4.84813681e-6
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
    # Old version uses 0.002 as constant value for Length of Day (LOD)
    # omegaEarth <- c(0, 0, 7.29211514670698e-05 * (1.0  - 0.002/86400.0))
    # Now changed to get exact value from EOP tables
    MJD <- trunc(julianDate - 2400000.5)
    LOD <- asteRiskData::earthPositions[asteRiskData::earthPositions[,4] == MJD, 8]
    omegaEarth <- c(0, 0, 7.29211514670698e-05 * (1.0  - LOD/86400.0))
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
    a <- earthRadius_WGS84
    a2 <- a^2
    f <- earthFlatteningFactor_WGS84
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

LATLONtoECEF <- function(position_LATLON, degreesInput=TRUE) {
    lat <- position_LATLON[1]
    lon <- position_LATLON[2]
    alt <- position_LATLON[3]
    if(degreesInput) {
        lat <- deg2rad(lat)
        lon <- deg2rad(lon)
    }
    # Prime-vertical radius of curvature
    N <- earthRadius_WGS84/(sqrt(1 - earthEccentricity_WGS84^2 * sin(lat)^2))
    position_ECEF <- c((N + alt) * cos(lat) * cos(lon),
                       (N + alt) * cos(lat) * sin(lon),
                       ((1 - earthEccentricity_WGS84^2) * N + alt) * sin(lat))
    return(unname(position_ECEF))
}

TEMEtoLATLON <- function(position_TEME, dateTime, degreesOutput=TRUE) {
    ECEFcoords <- TEMEtoECEF(position_TEME=position_TEME, dateTime=dateTime)
    return(ECEFtoLATLON(ECEFcoords$position, degreesOutput=degreesOutput))
}

ECEFtoGCRF <- function(position_ECEF, velocity_ECEF=c(0, 0, 0), dateTime) {
    hasData()
    date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
    year <- date$year + 1900
    month <- date$mon + 1
    day <- date$mday
    hour <- date$hour
    minute <- date$min
    second <- date$sec
    Mjd_UTC <- iauCal2jd(year, month, day, hour, minute, second)$DATE
    results <- ECEFtoECI(Mjd_UTC, c(position_ECEF, velocity_ECEF))
    return(list(
        position=as.numeric(results$position),
        velocity=as.numeric(results$velocity)
    ))
}

GCRFtoECEF <- function(position_GCRF, velocity_GCRF=c(0, 0, 0), dateTime) {
    hasData()
    date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
    year <- date$year + 1900
    month <- date$mon + 1
    day <- date$mday
    hour <- date$hour
    minute <- date$min
    second <- date$sec
    Mjd_UTC <- iauCal2jd(year, month, day, hour, minute, second)$DATE
    results <- ECItoECEF(Mjd_UTC, c(position_GCRF, velocity_GCRF))
    return(list(
        position=as.numeric(results$position),
        velocity=as.numeric(results$velocity)
    ))
}

TEMEtoGCRF <- function(position_TEME, velocity_TEME=c(0,0,0), dateTime) {
    hasData()
    ecef_results <- TEMEtoECEF(position_TEME, velocity_TEME, dateTime)
    GCRF_results <- ECEFtoGCRF(ecef_results$position, ecef_results$velocity, dateTime)
    return(GCRF_results)
}

TEMEtoGCRF2 <- function(position_TEME, velocity_TEME=c(0,0,0), dateTime) {
    hasData()
    ecef_results <- TEMEtoECEF(position_TEME, velocity_TEME, dateTime)
    GCRF_results <- ECEFtoGCRF(ecef_results$position, ecef_results$velocity, dateTime)
    return(GCRF_results)
}


GCRFtoLATLON <- function(position_GCRF, dateTime, degreesOutput=TRUE) {
    hasData()
    ECEFcoords <- GCRFtoECEF(position_GCRF=position_GCRF, dateTime=dateTime)
    return(ECEFtoLATLON(ECEFcoords$position, degreesOutput=degreesOutput))
}

LATLONtoGCRF <- function(position_LATLON, dateTime, degreesInput=TRUE) {
    hasData()
    ECEFcoords <- LATLONtoECEF(position_LATLON, degreesInput=degreesInput)
    return(ECEFtoGCRF(ECEFcoords, dateTime=dateTime))
}

KOEtoECI <- function(a, e, i, M, omega, OMEGA, keplerAccuracy=10e-8, maxKeplerIterations=100) {
    # calculate true anomaly from mean anomaly
    convergence <- FALSE
    iterations <- 0
    Eomega <- M
    while(!convergence) {
        iterations <- iterations + 1
        delta_kepler_sol <- - ( (Eomega - e*sin(Eomega) - M) / (1 - e*cos(Eomega)) )
        Eomega <- Eomega + delta_kepler_sol
        if(iterations > maxKeplerIterations | delta_kepler_sol < keplerAccuracy) {
            convergence <- TRUE
        }
    }
    nu <- 2 * atan2(sqrt(1 + e) * sin(Eomega/2), sqrt(1 - e) * cos(Eomega/2))
    R <- a * (1 - e*cos(Eomega))
    orbital_position <- R * c(cos(nu), sin(nu), 0)
    orbital_velocity <- (sqrt(GM_Earth_TCB * a)/R) * c(-sin(Eomega), sqrt(1-e^2) * cos(Eomega), 0)
    eci_position <- c(orbital_position[1] * (cos(omega) * cos(OMEGA) - sin(omega) * cos(i) * sin(OMEGA)) -
                          orbital_position[2] * (sin(omega) * cos(OMEGA) + cos(omega) * cos(i) * sin(OMEGA)),
                      orbital_position[1] * (cos(omega) * sin(OMEGA) + sin(omega) * cos(i) * cos(OMEGA)) +
                          orbital_position[2] * (cos(omega) * cos(i) * cos(OMEGA) - sin(omega) * sin(OMEGA)),
                      orbital_position[1] * sin(omega) * sin(i) + orbital_position[2] * cos(omega) * sin(i))
    eci_speed <- c(orbital_velocity[1] * (cos(omega) * cos(OMEGA) - sin(omega) * cos(i) * sin(OMEGA)) -
                       orbital_velocity[2] * (sin(omega) * cos(OMEGA) + cos(omega) * cos(i) * sin(OMEGA)),
                   orbital_velocity[1] * (cos(omega) * sin(OMEGA) + sin(omega) * cos(i) * cos(OMEGA)) +
                       orbital_velocity[2] * (cos(omega) * cos(i) * cos(OMEGA) - sin(omega) * sin(OMEGA)),
                   orbital_velocity[1] * sin(omega) * sin(i) + orbital_velocity[2] * cos(omega) * sin(i))
    return(list(
        position=eci_position,
        velocity=eci_speed
    ))
}


ECItoKOE <- function(position_ECI, velocity_ECI) {
    # calculate orbital momentum
    eps <- .Machine$double.eps
    h <- vectorCrossProduct3D(position_ECI, velocity_ECI)
    e_vector <- vectorCrossProduct3D(velocity_ECI, h)/GM_Earth_TCB - position_ECI/sqrt(sum(position_ECI^2))
    # This is equivalent to the following expression by Vallado 2007
    # e_vector2 <- ((sum(velocity_ECI^2) - GM_Earth_TCB/sqrt(sum(position_ECI^2)))*position_ECI - 
    # (position_ECI%*%velocity_ECI)*velocity_ECI )/GM_Earth_TCB
    e <- sqrt(sum(e_vector^2))
    node_vector <- c(-h[2], h[1], 0)
    E <- sum(velocity_ECI^2)/2 - GM_Earth_TCB/sqrt(sum(position_ECI^2))
    if(abs(E) > eps) {
        a <- -GM_Earth_TCB/(2*E)
        # p <- a*(1-e^2)
    } else {
        # p <- sum(h^2)/GM_Earth_TCB
        a <- Inf
    }
    # Vallado´s implementation defines p always as follows
    p <- sum(h^2)/GM_Earth_TCB
    i <- acos(h[3]/sqrt(sum(h^2)))
    # determine special orbit cases
    orbitType <- "ei" # general case: non-circular (elliptical) orbit with inclination
    if(e < eps) { # almost 0 eccentricity: circular orbits
        if(i < eps | abs(i - pi) < eps) { # no inclination: equatorial orbits
            orbitType <- "ce"
        } else {
            orbitType <- "ci" # circular inclined
        }
    } else if(i < eps | abs(i - pi) < eps) { # elliptical equatorial orbit
        orbitType <- "ee"
    }
    if(sqrt(sum(node_vector^2)) > eps) {
        cosOMEGA <- node_vector[1]/sqrt(sum(node_vector^2))
        if(abs(cosOMEGA) > 1) {
            cosOMEGA <- sign(cosOMEGA)
        }
        OMEGA <- acos(cosOMEGA)
        if(node_vector[2] < 0) {
            OMEGA <- 2*pi - OMEGA
        }
    } else {
        OMEGA <- NaN
    }
    if(orbitType == "ei") {
        omega <- acos((node_vector%*%e_vector)/(sqrt(sum(node_vector^2))*e))
        if(e_vector[3] < 0) {
            omega <- 2*pi - omega
        }
    } else {
        omega <- NaN
    }
    if(orbitType %in% c("ei", "ee")) {
        nu <- acos((e_vector%*%position_ECI)/(e*sqrt(sum(position_ECI^2))))
        if(position_ECI %*% velocity_ECI < 0) {
            nu <- 2*pi - nu
        }
    } else {
        nu <- NaN
    }
    # non-standard orbital parameters
    # argument of latitude - for non-equatorial orbits
    if(orbitType %in% c("ci", "ei")) {
        arglat <- acos((node_vector%*%position_ECI)/(sqrt(sum(node_vector^2))*sqrt(sum(position_ECI^2))))
        if(position_ECI[3] < 0) {
            arglat <- 2*pi - arglat
        }
    } else {
        arglat <- NaN
    }
    # longitude of perigee - elliptical equatorial orbits
    if(orbitType == "ee") {
        cosLonPer <- e_vector[1]/e
        if(abs(cosLonPer) > 1) {
            cosLonPer <- sign(cosLonPer)
        }
        lonPer <- acos(cosLonPer)
        if(e_vector[2] < 0) {
            lonPer <- 2*pi - lonPer
        }
        if(i > 0.5*pi) {
            lonPer <- 2*pi - lonPer
        }
    } else {
        lonPer <- NaN
    }
    # true longitude - circular equatorial orbits
    if((sqrt(sum(position_ECI^2)) > eps) & orbitType == "ce") {
        cosTrueLon <- position_ECI[1] / sqrt(sum(position_ECI))
        if(abs(cosTrueLon) > 1) {
            cosTrueLon <- sign(cosTrueLon)
        }
        trueLon <- acos(cosTrueLon)
        if(position_ECI[2] < 0) {
            trueLon <- 2*pi - trueLon
        }
        if(i > 0.5*pi) {
            trueLon <- 2*pi - trueLon
        }
    } else {
        trueLon <- NaN
    }
    # calculate mean anomaly for non-circular orbits
    if(orbitType %in% c("ei", "ee")) {
        if(e < 1 - eps) { 
            # truly elliptical orbit
            sine <- (sqrt(1-e^2)*sin(nu)) / (1 + e*cos(nu))
            cose <- (e + cos(nu)) / (1 + e*cos(nu))
            e0 <- atan2(sine, cose)
            m <- e0 - e*sin(e0)
        } else if((e > 1 + eps) & (abs(nu) + 1e-5 < pi - acos(1/e))) { 
            # hyperbolic orbit
            sine <- (sqrt(-1+e^2)*sin(nu)) / (1 + e*cos(nu))
            e0 <- asinh(sine)
            m <- e*sinh(e0) - e0
        } else if(abs(nu) < 168*pi/180) {
            # parabolic orbit
            e0 <- tan(nu/2)
            m <- e0 + e0^3/3
        }
        if(e < 1) {
            m <- rem(m, 2*pi)
            if(m < 0) {
                m <- m + 2*pi
            }
            e0 <- rem(e0, 2*pi)
        }
    } else {
        m <- NaN
    }
    return(list(
        semiMajorAxis = a,
        eccentricity = e,
        inclination = i,
        meanAnomaly = as.vector(m),
        argumentPerigee = as.vector(omega),
        longitudeAscendingNode = OMEGA,
        trueAnomaly = as.vector(nu),
        argumentLatitude = as.vector(arglat),
        longitudePerigee = lonPer,
        trueLongitude = trueLon
    ))
}

### New coordinate transformations with quaternion rotations

rotationMODtoGCRF <- function(MJD_UTC) {
    MJD_TT <- MJDUTCtoMJDTT(MJD_UTC)
    precession <- IAU76_precession(MJD_TT)
    rotation <- anglesToQuaternion(c(precession$z,
                                     -precession$theta,
                                     precession$zeta),
                                   "ZYZ")
    return(rotation)
}

rotationGCRFtoMOD <- function(MJD_UTC) {
    inverseRotation <- rotationMODtoGCRF(MJD_UTC)
    return(Conj(inverseRotation))
}

rotationTEMEtoMOD <- function(MJD_UTC, deltaDeltaPsi = 0, deltaDeltaEps = 0) {
    MJD_TT <- MJDUTCtoMJDTT(MJD_UTC)
    nutation <- IAU76_nutation(MJD_TT)
    nutation$deltaPsi <- nutation$deltaPsi + deltaDeltaPsi
    nutation$deltaEps <- nutation$deltaEps + deltaDeltaEps
    obliquity <-  nutation$meanEclipticObliquity + nutation$deltaEps
    T_TT <- (MJD_TT - MJD_J2000)/36525 # Julian centuries in TT
    meanLongitudeAscendingNodeMoon <- (125.04452222 + T_TT * (-(5*360 + 134.1362608) + T_TT * (0.0020708 + T_TT * 2.2e-6))) * pi/180
    meanLongitudeAscendingNodeMoon <- meanLongitudeAscendingNodeMoon %% (2*pi)
    equationEquinoxes1982 <- nutation$deltaPsi * cos(nutation$meanEclipticObliquity) +
        (0.00264 * sin(meanLongitudeAscendingNodeMoon) + 0.000063 * sin(2 * meanLongitudeAscendingNodeMoon)) * pi/648000
    quatTEMEtoTOD <- anglesToQuaternion(-equationEquinoxes1982, "Z")
    quatTODtoMOD <- anglesToQuaternion(c(obliquity, nutation$deltaPsi, -nutation$meanEclipticObliquity), "XZX")
    return(quatTEMEtoTOD * quatTODtoMOD)
}

rotationMODtoTEME <- function(MJD_UTC, deltaDeltaPsi = 0, deltaDeltaEps = 0) {
    inverseRotation <- rotationTEMEtoMOD(MJD_UTC, deltaDeltaPsi, deltaDeltaEps)
    return(Conj(inverseRotation))
}

rotationGCRFtoTEME <- function(MJD_UTC) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "l")
    quatGCRFtoMOD <- rotationGCRFtoMOD(MJD_UTC)
    quatMODtoTEME <- rotationMODtoTEME(MJD_UTC, IERS_results$dpsi, IERS_results$deps)
    return(quatGCRFtoMOD * quatMODtoTEME)
}

rotationTEMEtoGCRF <- function(MJD_UTC) {
    inverseRotation <- rotationGCRFtoTEME(MJD_UTC, deltaDeltaPsi, deltaDeltaEps)
    return(Conj(inverseRotation))
}

rotationPEFtoMOD <- function(MJD_UTC, deltaDeltaPsi = 0, deltaDeltaEps = 0) {
    MJD_TT <- MJDUTCtoMJDTT(MJD_UTC)
    MJD_UT1 <- MJDUTCtoMJDUT1(MJD_UTC)
    nutation <- IAU76_nutation(MJD_TT)
    nutation$deltaPsi <- nutation$deltaPsi + deltaDeltaPsi
    nutation$deltaEps <- nutation$deltaEps + deltaDeltaEps
    obliquity <-  nutation$meanEclipticObliquity + nutation$deltaEps
    T_TT <- (MJD_TT - MJD_J2000)/36525 # Julian centuries in TT
    meanLongitudeAscendingNodeMoon <- (125.04452222 + T_TT * (-(5*360 + 134.1362608) + T_TT * (0.0020708 + T_TT * 2.2e-6))) * pi/180
    meanLongitudeAscendingNodeMoon <- meanLongitudeAscendingNodeMoon %% (2*pi)
    equationEquinoxes1982 <- nutation$deltaPsi * cos(nutation$meanEclipticObliquity) +
        (0.00264 * sin(meanLongitudeAscendingNodeMoon) + 0.000063 * sin(2 * meanLongitudeAscendingNodeMoon)) * pi/648000
    thetaGMST <- MJDToGMST(MJD_UT1)
    thetaGAST <- thetaGMST + equationEquinoxes1982
    quatPEFtoTOD <- anglesToQuaternion(-thetaGAST, "Z")
    quatTODtoMOD <- anglesToQuaternion(c(obliquity, nutation$deltaPsi, -nutation$meanEclipticObliquity), "XZX")
    return(quatPEFtoTOD * quatTODtoMOD)
}

rotationMODtoPEF <- function(MJD_UTC, deltaDeltaPsi = 0, deltaDeltaEps = 0) {
    inverseRotation <- rotationPEFtoMOD(MJD_UTC, deltaDeltaPsi, deltaDeltaEps)
    return(Conj(inverseRotation))
}

rotationGCRFtoJ2000 <- function(MJD_UTC) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "l")
    quatGCRFtoMOD <- rotationGCRFtoMOD(MJD_UTC)
    quatMODtoPEF <- rotationMODtoPEF(MJD_UTC, IERS_results$dpsi, IERS_results$deps)
    quatPEFtoMOD <- rotationPEFtoMOD(MJD_UTC, 0, 0)
    quatMODtoJ2000 <- rotationMODtoGCRF(MJD_UTC)
    return(quatGCRFtoMOD * quatMODtoPEF * quatPEFtoMOD * quatMODtoJ2000)
}

rotationJ2000toGCRF <- function(MJD_UTC) {
    inverseRotation <- rotationGCRFtoJ2000(MJD_UTC)
    return(Conj(inverseRotation))
}