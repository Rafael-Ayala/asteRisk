MJDUTCtoMJDTT <- function(MJD_UTC) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "n")
    leapSeconds <- IERS_results$TAI_UTC
    MJD_TT <- MJD_UTC + (leapSeconds + 32.184)/86400
    return(MJD_TT)
}

MJDTTtoMJDUTC <- function(MJD_TT) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_TT, interp = "n")
    leapSeconds <- IERS_results$TAI_UTC
    MJD_UTC <- MJD_TT - (leapSeconds + 32.184)/86400
    return(MJD_TT)
}

MJDUTCtoMJDUT1 <- function(MJD_UTC) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "n")
    MJD_UT1 <- MJD_UTC + IERS_results$UT1_UTC/86400
    return(MJD_UT1)
}

MJDUT1toMJDUTC <- function(MJD_UT1) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UT1, interp = "n")
    MJD_UTC <- MJD_UT1 - IERS_results$UT1_UTC/86400
    return(MJD_UTC)
}

MJDUTCtoMJDTDB <- function(MJD_UTC) {
    MJD_TT <- MJDUTCtoMJDTT(MJD_UTC)
    MJD_TDB <- Mjday_TDB(MJD_TT)
    return(MJD_TT)
}

IAU76_precession <- function(targetTime, inputInSeconds = FALSE) {
    if(inputInSeconds) {
        T_TT <- targetTime/3155760000
    } else {
        T_TT <- (targetTime - MJD_J2000)/36525 # Julian centuries in TT
    }
    zeta <- T_TT * (2306.2181 + T_TT * (0.30188 + T_TT * 0.017998)) * pi/648000
    theta <- T_TT * (2004.3109 + T_TT * (-0.42665 + T_TT * -0.041833)) * pi/648000
    z <- T_TT * (2306.2181 + T_TT * (1.09468 + T_TT * 0.018203)) * pi/648000
    ts <- 1/3155760000
    # added expressions for derivatives from SPICE's zzeprc76
    dzeta <- ts * (2306.2181 + T_TT * (0.60375999999999996 + T_TT * 3 * 0.017998))*pi/(180*3600)
    dtheta <- ts * (2004.3109 + T_TT * (-0.85329999999999995 + T_TT * 3 * -0.041833))*pi/(180*3600)
    dz <- ts * (2306.2181 + T_TT * (2.1893600000000002 + T_TT * 3 * 0.018203))*pi/(180*3600)
    zeta <- zeta %% (2*pi)
    theta <- theta %% (2*pi)
    z <- z %% (2*pi)
    return(list(zeta=zeta,
                theta=theta,
                z=z,
                dZeta=dzeta,
                dTheta=dtheta,
                dZ=dz))
}

IAU76_nutation <- function(targetTime, numberCoefficients = 106, inputInSeconds = FALSE) {
    ## target time can be either MJD_TT (i.e. modified julian date days) or directly 
    ## seconds since J2000, for exact matching with SPICE's routines
    ## this is matching SPICE's 1980 nutation routines
    if(!(is.numeric(numberCoefficients) & numberCoefficients >= 1 & numberCoefficients <= 106)) {
        numberCoefficients <- 106
        warning(strwrap("Invalid number of nutation coefficients provided. All
                         106 coefficients will be used."))
    }
    if(inputInSeconds) {
        T_TT <- targetTime/3155760000
    } else {
        T_TT <- (targetTime - MJD_J2000)/36525 # Julian centuries in TT
    }
    ts <- 1/31557600
    meanEclipticObliquity <- (23.439291 + T_TT * (-0.0130042 + T_TT * (-1.64e-7 + T_TT * 5.04e-7))) * pi/180
    dMeanEclipticObliquity <- ts * (-46.815 + T_TT * (-0.0011800000000000001 + T_TT * 3 * 0.001813)) * pi/64800000 
    meanEclipticObliquity <- meanEclipticObliquity %% (2*pi)
    ## Delaunay orbital elements of Sun and Moon, in radians
    meanAnomalyMoon <- (134.96298139 + T_TT * ((1325*360 + 198.8673981) + T_TT * (0.0086972 + T_TT * 1.78e-5))) * pi/180
    meanAnomalySun <- (357.52772333 + T_TT * ((99*360 + 359.0503400) + T_TT * (-0.0001603 + T_TT * -3.3e-6))) * pi/180
    # Mu is the difference between mean longitude of the Moon and mean longitude 
    # of the ascending node of the Moon. Also usually named F
    muMoon <- (93.27191028 + T_TT * ((1342*360 + 82.0175381) + T_TT * (-0.0036825 + T_TT * 3.1e-6))) * pi/180
    meanElongationMoonSun <- (297.85036306 + T_TT * ((1236*360 + 307.1114800) + T_TT * (-0.0019142 + T_TT * 5.3e-6))) * pi/180
    meanLongitudeAscendingNodeMoon <- (125.04452222 + T_TT * (-(5*360 + 134.1362608) + T_TT * (0.0020708 + T_TT * 2.2e-6))) * pi/180
    meanAnomalyMoon <- meanAnomalyMoon %% (2*pi)
    meanAnomalySun <- meanAnomalySun %% (2*pi)
    muMoon <- muMoon %% (2*pi)
    # we also need derivatives to match all the output given by SPICE's zzwahr
    # just take derivatives of the polynomials in Horner's form
    dMeanAnomalyMoon <- ts * ((1325*360 + 198.8673981) + T_TT * (2* 0.0086972 + 3* T_TT * 1.78e-5)) * pi/180
    dMeanAnomalySun <- ts * ((99*360 + 359.0503400) + T_TT * (2*-0.0001603 + 3* T_TT * -3.3e-6)) * pi/180
    dMuMoon <- ts * ((1342*360 + 82.0175381) + T_TT * (2 * -0.0036825 + 3* T_TT * 3.1e-6)) * pi/180
    dMeanElongationMoonSun <- ts * ((1236*360 + 307.1114800) + T_TT * (2*-0.0019142 + 3* T_TT * 5.3e-6)) * pi/180
    dMeanLongitudeAscendingNodeMoon <- ts * (-(5*360 + 134.1362608) + T_TT * (2* 0.0020708 + 3* T_TT * 2.2e-6)) * pi/180
    # end of derivatives of Delaunay orbital elements
    meanElongationMoonSun <- meanElongationMoonSun %% (2*pi)
    meanLongitudeAscendingNodeMoon <- meanLongitudeAscendingNodeMoon %% (2*pi)
    api <- asteRiskData::nut_IAU1980[1:numberCoefficients ,1] * meanAnomalyMoon +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,2] * meanAnomalySun +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,3] * muMoon +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,4] * meanElongationMoonSun +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,5] * meanLongitudeAscendingNodeMoon
    dApi <- asteRiskData::nut_IAU1980[1:numberCoefficients ,1] * dMeanAnomalyMoon +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,2] * dMeanAnomalySun +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,3] * dMuMoon +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,4] * dMeanElongationMoonSun +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,5] * dMeanLongitudeAscendingNodeMoon
    deltaPsi <- sum((asteRiskData::nut_IAU1980[1:numberCoefficients ,7] +
                         asteRiskData::nut_IAU1980[1:numberCoefficients ,8] * T_TT) *
                        sin(api))
    deltaEps <- sum((asteRiskData::nut_IAU1980[1:numberCoefficients ,9] +
                         asteRiskData::nut_IAU1980[1:numberCoefficients ,10] * T_TT) *
                        cos(api))
    dDeltaPsi <- sum((asteRiskData::nut_IAU1980[1:numberCoefficients ,7] +
                          asteRiskData::nut_IAU1980[1:numberCoefficients ,8] * T_TT) *
                         cos(api) * dApi)
    dDeltaEps <- sum((asteRiskData::nut_IAU1980[1:numberCoefficients ,9] +
                          asteRiskData::nut_IAU1980[1:numberCoefficients ,10] * T_TT) *
                         -sin(api) * dApi)
    # Convert from tenths of milliarcseconds to radians
    deltaPsi <- deltaPsi * pi /(10000 * 3600 * 180)
    deltaEps <- deltaEps * pi /(10000 * 3600 * 180)
    dDeltaPsi <- dDeltaPsi * pi /(10000 * 3600 * 180)/100
    dDeltaEps <- dDeltaEps * pi /(10000 * 3600 * 180)/100
    return(list(meanEclipticObliquity = meanEclipticObliquity,
                deltaPsi = deltaPsi,
                deltaEps = deltaEps,
                dMeanEclipticObliquity = dMeanEclipticObliquity,
                dDeltaPsi = dDeltaPsi,
                dDeltaEps = dDeltaEps
    ))
}

anglesToQuaternion <- function(angles, axes) {
    if(!(length(angles) >= 1 & length(angles) <= 3)) {
        stop("Please provide between 1 and 3 Euler angles describing the rotation")
    }
    if(length(angles) != nchar(axes)) {
        stop("The number of angles should match the number of rotation axes.")
    }
    axes <- tolower(axes)
    if(length(angles) == 1) {
        if(!(axes) %in% c("x", "y", "z")) {
            stop(strwrap("Please provide a valid rotation X, Y and Z."))
        }
        sin1 <- sin(angles/2)
        cos1 <- cos(angles/2)
        if(cos1 < 0) {
            sin1 <- -sin1
            cos1 <- -cos1
        }
        if(axes == "x") {
            rotationQuat <- quaternion(Re = cos1, i = sin1, j = 0, k = 0)
        } else if (axes == "y") {
            rotationQuat <- quaternion(Re = cos1, i = 0, j = sin1, k = 0)
        } else if (axes == "z") {
            rotationQuat <- quaternion(Re = cos1, i = 0, j = 0, k = sin1)
        }
    } else if(length(angles) == 2) {
        if(!(axes) %in% c("zy", "xy", "xz", "yx", "yz", "zx", "zy")) {
            stop(strwrap("Please provide a valid rotation sequence. Possible 
                          values for 2-angle rotations are ZY, XY, XZ, YX, YZ, 
                          ZX and ZY."))
        }
        sin1 <- sin(angles[1]/2)
        cos1 <- cos(angles[1]/2)
        sin2 <- sin(angles[2]/2)
        cos2 <- cos(angles[2]/2)
        q0 <- cos1 * cos2
        q0_sign <- sign(q0)
        if(axes == "zy") {
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * -sin1 * sin2, 
                                       j = q0_sign * cos1 * sin2, 
                                       k = q0_sign * sin1 * cos2)
        } else if(axes == "xy") {
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * sin1 * cos2, 
                                       j = q0_sign * cos1 * sin2, 
                                       k = q0_sign * sin1 * sin2)
        } else if(axes == "xz") {
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * sin1 * cos2, 
                                       j = q0_sign * -sin1 * sin2, 
                                       k = q0_sign * cos1 * sin2)
        } else if(axes == "yx") {
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * cos1 * sin2, 
                                       j = q0_sign * sin1 * cos2, 
                                       k = q0_sign * -sin1 * sin2)
        } else if(axes == "yz") {
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * sin1 * sin2, 
                                       j = q0_sign * sin1 * cos2, 
                                       k = q0_sign * cos1 * sin2)
        } else if(axes == "zx") {
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * cos1 * sin2, 
                                       j = q0_sign * sin1 * sin2, 
                                       k = q0_sign * sin1 * cos2)
        } else if(axes == "zy") {
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * -sin1 * sin2, 
                                       j = q0_sign * cos1 * sin2, 
                                       k = q0_sign * sin1 * cos2)
        }
    } else if(length(angles) == 3) {
        if(!(axes) %in% c("zyx", "xyx", "xyz", "xzx", "xzy", "yxy", "yxz", "yzx", 
                          "yzy", "zxy", "zxz", "zyz")) {
            stop(strwrap("Please provide a valid rotation sequence. Possible 
                          values for 3-angle rotations are ZYX, XYX, XYZ, XZX,
                          XZY, YXY, YXZ, YZX, YZY, ZXY, ZXZ and ZYZ."))
        }
        sin1 <- sin(angles[1]/2)
        cos1 <- cos(angles[1]/2)
        sin2 <- sin(angles[2]/2)
        cos2 <- cos(angles[2]/2)
        sin3 <- sin(angles[3]/2)
        cos3 <- cos(angles[3]/2)
        if(axes == "zyx") {
            q0 <- cos1 * cos2 * cos3 + sin1 * sin2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * cos2 * sin3 - sin1 * sin2 * cos3), 
                                       j = q0_sign * (cos1 * sin2 * cos3 + sin1 * cos2 * sin3), 
                                       k = q0_sign * (sin1 * cos2 * cos3 - cos1 * sin2 * sin3))
        } else if(axes == "xyx") {
            q0 <- cos1 * cos2 * cos3 - sin1 * cos2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * cos2 * sin3 + sin1 * cos2 * cos3), 
                                       j = q0_sign * (cos1 * sin2 * cos3 + sin1 * sin2 * sin3), 
                                       k = q0_sign * (sin1 * sin2 * cos3 - cos1 * sin2 * sin3))
        } else if(axes == "xyz") {
            q0 <- cos1 * cos2 * cos3 - sin1 * sin2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (sin1 * cos2 * cos3 + cos1 * sin2 * sin3), 
                                       j = q0_sign * (cos1 * sin2 * cos3 - sin1 * cos2 * sin3), 
                                       k = q0_sign * (cos1 * cos2 * sin3 + sin1 * sin2 * cos3))
        } else if(axes == "xzx") {
            q0 <- cos1 * cos2 * cos3 - sin1 * cos2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * cos2 * sin3 + sin1 * cos2 * cos3), 
                                       j = q0_sign * (cos1 * sin2 * sin3 - sin1 * sin2 * cos3), 
                                       k = q0_sign * (cos1 * sin2 * cos3 + sin1 * sin2 * sin3))
        } else if(axes == "xzy") {
            q0 <- cos1 * cos2 * cos3 + sin1 * sin2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (sin1 * cos2 * cos3 - cos1 * sin2 * sin3), 
                                       j = q0_sign * (cos1 * cos2 * sin3 - sin1 * sin2 * cos3), 
                                       k = q0_sign * (cos1 * sin2 * cos3 + sin1 * cos2 * sin3))
        } else if(axes == "yxy") {
            q0 <- cos1 * cos2 * cos3 - sin1 * cos2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * sin2 * cos3 + sin1 * sin2 * sin3), 
                                       j = q0_sign * (cos1 * cos2 * sin3 + sin1 * cos2 * cos3), 
                                       k = q0_sign * (cos1 * sin2 * sin3 - sin1 * sin2 * cos3))
        } else if(axes == "yxz") {
            q0 <- cos1 * cos2 * cos3 + sin1 * sin2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * sin2 * cos3 + sin1 * cos2 * sin3), 
                                       j = q0_sign * (sin1 * cos2 * cos3 - cos1 * sin2 * sin3), 
                                       k = q0_sign * (cos1 * cos2 * sin3 - sin1 * sin2 * cos3))
        } else if(axes == "yzx") {
            q0 <- cos1 * cos2 * cos3 - sin1 * sin2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * cos2 * sin3 + sin1 * sin2 * cos3), 
                                       j = q0_sign * (sin1 * cos2 * cos3 + cos1 * sin2 * sin3), 
                                       k = q0_sign * (cos1 * sin2 * cos3 - sin1 * cos2 * sin3))
        } else if(axes == "yzy") {
            q0 <- cos1 * cos2 * cos3 - sin1 * cos2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (sin1 * sin2 * cos3 - cos1 * sin2 * sin3), 
                                       j = q0_sign * (cos1 * cos2 * sin3 + sin1 * cos2 * cos3), 
                                       k = q0_sign * (cos1 * sin2 * cos3 + sin1 * sin2 * sin3))
        } else if(axes == "zxy") {
            q0 <- cos1 * cos2 * cos3 - sin1 * sin2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * sin2 * cos3 - sin1 * cos2 * sin3), 
                                       j = q0_sign * (cos1 * cos2 * sin3 + sin1 * sin2 * cos3), 
                                       k = q0_sign * (sin1 * cos2 * cos3 + cos1 * sin2 * sin3))
        } else if(axes == "zxz") {
            q0 <- cos1 * cos2 * cos3 - sin1 * cos2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * sin2 * cos3 + sin1 * sin2 * sin3), 
                                       j = q0_sign * (sin1 * sin2 * cos3 - cos1 * sin2 * sin3), 
                                       k = q0_sign * (cos1 * cos2 * sin3 + sin1 * cos2 * cos3))
        } else if(axes == "zyz") {
            q0 <- cos1 * cos2 * cos3 - sin1 * cos2 * sin3
            q0_sign <- sign(q0)
            rotationQuat <- quaternion(Re = abs(q0), 
                                       i = q0_sign * (cos1 * sin2 * sin3 - sin1 * sin2 * cos3), 
                                       j = q0_sign * (cos1 * sin2 * cos3 + sin1 * sin2 * sin3), 
                                       k = q0_sign * (cos1 * cos2 * sin3 + sin1 * cos2 * cos3))
        }
    }
    return(rotationQuat)
}

quaternionToDCM <- function(quaternion) {
    q0 <- Re(quaternion)
    q1 <- i(quaternion)
    q2 <- j(quaternion)
    q3 <- k(quaternion)
    return(matrix(c(q0^2 + q1^2 - q2^2 -q3^2, 2*(q1*q2 + q0 * q3), 2*(q1*q3 - q0*q2),
                    2*(q1*q2 - q0*q3), q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 + q0 * q1),
                    2*(q1*q3 + q0 * q2), 2*(q2*q3 - q0*q1), q0^2 - q1^2 - q2^2 + q3^2),
                  nrow=3, byrow = TRUE))
}

latitudeToGeocentric <- function(geodeticLatitude, degreesInput=TRUE,
                                 degreesOutput=TRUE, flatteningFactor=earthFlatteningFactor_WGS84) {
    if(degreesInput) {
        geodeticLatitude <- deg2rad(geodeticLatitude)
    }
    geocentricLatitude <- atan((1 - flatteningFactor)^2*tan(geodeticLatitude))
    if(degreesOutput) {
        geocentricLatitude <- rad2deg(geocentricLatitude)
    }
    return(geocentricLatitude)
}


eulerAngleRatesToAngularVelocity <- function(angles, angleRates, axes) {
    if(length(angles) != 3 || length(angleRates) != 3) {
        stop("Please provide 3 Euler angles describing the rotation
             and their rates (derivatives)")
    }
    if(nchar(axes) != 3) {
        stop("The rotation sequence should comprise 3 axes")
    }
    axes <- tolower(axes)
    if(!(axes) %in% c("zyx", "xyx", "xyz", "xzx", "xzy", "yxy", "yxz", "yzx", 
                      "yzy", "zxy", "zxz", "zyz")) {
        stop(strwrap("Please provide a valid 3-angle rotation sequence. Possible 
                          values for 3-angle rotations are ZYX, XYX, XYZ, XZX,
                          XZY, YXY, YXZ, YZX, YZY, ZXY, ZXZ and ZYZ."))
    }
    s1 <- sin(angles[1])
    c1 <- cos(angles[1])
    s2 <- sin(angles[2])
    c2 <- cos(angles[2])
    s3 <- sin(angles[3])
    c3 <- cos(angles[3])
    if(axes == "xyx") {
        Bmatrix <- matrix(c(c2, 0, 1,
                            s2*s3, c3, 0,
                            s2*c3, -s3, 0),
                          nrow = 3, byrow = TRUE)
    } else if(axes == "xyz") {
        Bmatrix <- matrix(c(c2*c3, s3, 0,
                            -c2*s3, c3, 0,
                            s2, 0, 1),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "xzx") {
        Bmatrix <- matrix(c(c2, 0, 1,
                            -s2*c3, s3, 0,
                            s2*s3, c3, 0),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "xzy") {
        Bmatrix <- matrix(c(c2*c3, -s3, 0,
                            -s2, 0, 1,
                            c2*s3, c3, 0),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "yxy") {
        Bmatrix <- matrix(c(s2*s3, c3, 0,
                            c2, 0, 1,
                            -s2*c3, s3, 0),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "yxz") {
        Bmatrix <- matrix(c(c2*s3, c3, 0,
                            c2*c3, -s3, 0,
                            -s2, 0, 1),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "yzx") {
        Bmatrix <- matrix(c(s2, 0, 1,
                            c2*c3, s3, 0,
                            -c2*s3, c3, 0),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "yzy") {
        Bmatrix <- matrix(c(s2*c3, -s3, 0,
                            c2, 0, 1,
                            s2*s3, c3, 0),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "zxy") {
        Bmatrix <- matrix(c(-c2*s3, c3, 0,
                            s2, 0, 1,
                            c2*c3, s3, 0),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "zxz") {
        Bmatrix <- matrix(c(s2*s3, c3, 0,
                            s2*c3, -s3, 0,
                            c2, 0, 1),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "zyx") {
        Bmatrix <- matrix(c(-s2, 0, 1,
                            c2*s3, c3, 0,
                            c2*c3, -s3, 0),
                          nrow = 3, byrow = TRUE)    
    } else if(axes == "zyz") {
        Bmatrix <- matrix(c(-s2*c3, s3, 0,
                            s2*s3, c3, 0,
                            c2, 0, 1),
                          nrow = 3, byrow = TRUE)    
    }
    return(drop(
        Bmatrix %*% c(angleRates)
    ))
}

eulerAnglesToDCM <- function(angles, axes) {
    if(length(angles) != 3) {
        stop("Please provide 3 Euler angles describing the rotation")
    }
    if(nchar(axes) != 3) {
        stop("The rotation sequence should comprise 3 axes")
    }
    axes <- tolower(axes)
    if(!(axes) %in% c("zyx", "xyx", "xyz", "xzx", "xzy", "yxy", "yxz", "yzx", 
                      "yzy", "zxy", "zxz", "zyz")) {
        stop(strwrap("Please provide a valid 3-angle rotation sequence. Possible 
                          values for 3-angle rotations are ZYX, XYX, XYZ, XZX,
                          XZY, YXY, YXZ, YZX, YZY, ZXY, ZXZ and ZYZ."))
    }
    s1 <- sin(angles[1])
    c1 <- cos(angles[1])
    s2 <- sin(angles[2])
    c2 <- cos(angles[2])
    s3 <- sin(angles[3])
    c3 <- cos(angles[3])
    if(axes == "xyz") {
        DCM <- matrix(c(c2*c3, -c2*s3, s2,
                        s1*s2*c3+c1*s3, -s1*s2*s3+c1*c3, -s1*c2,
                        -c1*s2*c3+s1*s3, c1*s2*s3+s1*c3, c1*c2),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "xzy") {
        DCM <- matrix(c(c2*c3, -s2, c2*s3,
                        c1*s2*c3+s1*s3, c1*c2, c1*s2*s3-s1*c3,
                        s1*s2*c3-c1*s3, s1*c2, s1*s2*s3+c1*c3),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "xyx") {
        DCM <- matrix(c(c2, s2*s3, s2*c3,
                        s1*s2, c1*c3-s1*c2*s3, -c1*s3-s1*c2*c3,
                        -c1*s2, s1*c3+c1*c2*s3, -s1*s3+c1*c2*c3),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "xzx") {
        DCM <- matrix(c(c2, -s2*c3, s2*s3,
                        c1*s2, c1*c2*c3-s1*s3, -c1*c2*s3-s1*c3,
                        s1*s2, s1*c2*c3+c1*s3, -s1*c2*s3+c1*c3),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "yxz") {
        DCM <- matrix(c(s1*s2*s3+c1*c3, s1*s2*c3-c1*s3, s1*c2,
                        c2*s3, c2*c3, -s2,
                        c1*s2*s3-s1*c3, c1*s2*c3+s1*s3, c1*c2),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "yzx") {
        DCM <- matrix(c(c1*c2, -c1*s2*c3+s1*s3, c1*s2*s3+s1*c3,
                        s2, c2*c3, -c2*s3,
                        -s1*c2, s1*s2*c3+c1*s3, -s1*s2*s3+c1*c3),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "yxy") {
        DCM <- matrix(c(-s1*c2*s3+c1*c3, s1*s2, s1*c2*c3+c1*s3,
                        s2*s3, c2, -s2*c3,
                        -c1*c2*s3-s1*c3, c1*s2, c1*c2*c3-s1*s3),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "yzy") {
        DCM <- matrix(c(c1*c2*c3-s1*s3, -c1*s2, c1*c2*s3+s1*c3,
                        s2*c3, c2, s2*s3,
                        -s1*c2*c3-c1*s3, s1*s2, -s1*c2*s3+c1*c3),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "zxy") {
        DCM <- matrix(c(-s1*s2*s3+c1*c3, -s1*c2, s1*s2*c3+c1*s3,
                        c1*s2*s3+s1*c3, c1*c2, -c1*s2*c3+s1*s3,
                        -c2*s3, s2, c2*c3),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "zyx") {
        DCM <- matrix(c(c1*c2, c1*s2*s3-s1*c3, c1*s2*c3+s1*s3,
                        s1*c2, s1*s2*s3+c1*c3, s1*s2*c3-c1*s3,
                        -s2, c2*s3, c2*c3),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "zxz") {
        DCM <- matrix(c(-s1*c2*s3+c1*c3, -s1*c2*c3-c1*c3, s1*s2,
                        c1*c2*s3+s1*c3, c1*c2*c3-s1*s3, -c1*s2,
                        s2*s3, s2*c3, c2),
                      nrow = 3, byrow = TRUE)
    } else if(axes == "zyz") {
        DCM <- matrix(c(c1*c2*c3 - s1*s3, -c1*c2*s3-s1*c3, c1*s2,
                        s1*c2*c3+c1*s3, -s1*c2*s3+c1*c3, s1*s2,
                        -s2*c3, s2*s3, c2),
                      nrow = 3, byrow = TRUE)
    } 
    return(DCM)
}

quaternionToDCM <- function(quaternion) {
    q0 <- Re(quaternion)
    q1 <- i(quaternion)
    q2 <- j(quaternion)
    q3 <- k(quaternion)
    return(matrix(c(q0^2 + q1^2 - q2^2 -q3^2, 2*(q1*q2 + q0 * q3), 2*(q1*q3 - q0*q2),
                    2*(q1*q2 - q0*q3), q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 + q0 * q1),
                    2*(q1*q3 + q0 * q2), 2*(q2*q3 - q0*q1), q0^2 - q1^2 - q2^2 + q3^2),
                  nrow=3, byrow = TRUE))
}

latitudeToGeocentric <- function(geodeticLatitude, degreesInput=TRUE,
                                 degreesOutput=TRUE, flatteningFactor=earthFlatteningFactor_WGS84) {
    if(degreesInput) {
        geodeticLatitude <- deg2rad(geodeticLatitude)
    }
    geocentricLatitude <- atan((1 - flatteningFactor)^2*tan(geodeticLatitude))
    if(degreesOutput) {
        geocentricLatitude <- rad2deg(geocentricLatitude)
    }
    return(geocentricLatitude)
}

zztwovxf <- function(axdef, indexa, plndef, indexp) {
    # adaptation of CSPICE's zztwovxf_
    i1 <- indexa
    if(i1 == 1) {
        i2 <- 2
        i3 <- 3
    } else if (i1 == 2) {
        i2 <- 3
        i1 <- 1
    } else if (i1 == 3) {
        i2 <- 1
        i3 <- 2
    }
    axdefNorm <- axdef/max(axdef)
    plndefNorm <- plndef/max(plndef)
    xformi1 <- dvhat(axdef[1:3], axdef[4:6])
    if(indexp == i2) {
        xformi3 <- unlist(vectorCrossProductDerivative3D_unitNorm(axdef[1:3], axdef[4:6], plndef[1:3], plndef[4:6]))
        xformi2 <- unlist(vectorCrossProductDerivative3D_unitNorm(xformi3[1:3], xformi3[4:6], axdef[1:3], axdef[4:6]))
    } else {
        xformi2 <- unlist(vectorCrossProductDerivative3D_unitNorm(plndef[1:3], plndef[4:6], axdef[1:3], axdef[4:6]))
        xformi3 <- unlist(vectorCrossProductDerivative3D_unitNorm(axdef[1:3], axdef[4:6], xformi2[1:3], xformi2[4:6]))
    }
    stateTransformationMatrix <- matrix(0, nrow=6, ncol=6)
    stateTransformationMatrix[,i1] <- xformi1
    stateTransformationMatrix[,i1+3] <- c(0,0,0,xformi1[1:3])
    stateTransformationMatrix[,i2] <- xformi2
    stateTransformationMatrix[,i2+3] <- c(0,0,0,xformi2[1:3])
    stateTransformationMatrix[,i3] <- xformi3
    stateTransformationMatrix[,i3+3] <- c(0,0,0,xformi3[1:3])
    return(stateTransformationMatrix)
}

