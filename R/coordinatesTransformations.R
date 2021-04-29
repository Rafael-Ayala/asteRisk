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

TEMEtoLATLON <- function(position_TEME, dateTime, degreesOutput=TRUE) {
    ECEFcoords <- TEMEtoECEF(position_TEME=position_TEME, dateTime=dateTime)
    return(ECEFtoLATLON(ECEFcoords$position, degreesOutput=degreesOutput))
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
    orbital_velocity <- (sqrt(earth_mu * a)/R) * c(-sin(Eomega), sqrt(1-e^2) * cos(Eomega), 0)
    eci_position <- c(orbital_position[1] * (cos(omega) * cos(OMEGA) - sin(omega) * cos(i) * sin(OMEGA)) -
                          orbital_position[2] * (sin(omega) * cos(OMEGA) + cos(omega) * cos(i) * sin(OMEGA)),
                      orbital_position[1] * (cos(omega) * sin(OMEGA) + sin(omega) * cos(i) * cos(OMEGA)) +
                          orbital_position[2] * (cos(omega) * cos(i) * cos(OMEGA) - sin(omega) * sin(OMEGA)),
                      orbital_position[1] * sin(omega) * sin(i) + orbital_position[2] * cos(omega) * sin(i))
    eci_speed <- c(orbital_velocity[1] * (cos(omega) * cos(OMEGA) - sin(omega) * cos(i) * sin(OMEGA)) -
                       orbital_velocity[2] * (sin(omega) * cos(OMEGA) + cos(omega) * cos(i) * sin(OMEGA)),
                   orbital_velocity[1] * (cos(omega) * sin(OMEGA) + sin(omega) * cos(i) * cos(OMEGA)) +
                       orbital_velocity[2] * (cos(omega) * cos(i) * cos(OMEGA) - sin(omega) * sin(OMEGA)),
                   orbital_velocity[1] * sin(omega) * sin(i) + orbital_position[2] * cos(omega) * sin(i))
    return(list(
        position=eci_position,
        velocity=eci_speed
    ))
}


ECItoKOE <- function(position_ECI, velocity_ECI) {
    # calculate orbital momentum
    eps <- .Machine$double.eps
    h <- vectorCrossProduct3D(position_ECI, velocity_ECI)
    e_vector <- vectorCrossProduct3D(velocity_ECI, h)/earth_mu - position_ECI/sqrt(sum(position_ECI^2))
    # This is equivalent to the following expression by Vallado 2007
    # e_vector2 <- ((sum(velocity_ECI^2) - earth_mu/sqrt(sum(position_ECI^2)))*position_ECI - 
    # (position_ECI%*%velocity_ECI)*velocity_ECI )/earth_mu
    e <- sqrt(sum(e_vector^2))
    node_vector <- c(-h[2], h[1], 0)
    E <- sum(velocity_ECI^2)/2 - earth_mu/sqrt(sum(position_ECI^2))
    if(abs(E) > eps) {
        a <- -earth_mu/(2*E)
        # p <- a*(1-e^2)
    } else {
        # p <- sum(h^2)/earth_mu
        a <- Inf
    }
    # Vallado´s implementation defines p always as follows
    p <- sum(h^2)/earth_mu
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
    if((sqrt(sum(position_ECI)) > eps) & orbitType == "ce") {
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
            m <- e0 - c*sin(e0)
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
        meanAnomaly = m,
        argumentPerigee = omega,
        longitudeAscendingNode = OMEGA,
        trueAnomaly = nu,
        argumentLatitude = arglat,
        longitudePerigee = lonPer,
        trueLongitude = trueLon
    ))
}