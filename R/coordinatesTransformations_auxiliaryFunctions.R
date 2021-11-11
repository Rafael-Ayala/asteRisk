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
    MJD_UT1 <- MJD_UT1 - IERS_results$UT1_UTC/86400
    return(MJD_UT1)
}

IAU76_precession <- function(MJD_TT) {
    T_TT <- (MJD_TT - MJD_J2000)/36525 # Julian centuries in TT
    zeta <- T_TT * (2306.2181 + T_TT * (0.30188 + T_TT * 0.017998)) * pi/648000
    theta <- T_TT * (2004.3109 + T_TT * (-0.42665 + T_TT * -0.041833)) * pi/648000
    z <- T_TT * (2306.2181 + T_TT * (1.09468 + T_TT * 0.018203)) * pi/648000
    zeta <- zeta %% (2*pi)
    theta <- theta %% (2*pi)
    z <- z %% (2*pi)
    return(list(zeta=zeta,
                theta=theta,
                z=z))
}

IAU76_nutation <- function(MJD_TT, numberCoefficients = 106) {
    if(!(is.numeric(numberCoefficients) & numberCoefficients >= 1 & numberCoefficients <= 106)) {
        numberCoefficients <- 106
        warning(strwrap("Invalid number of nutation coefficients provided. All
                         106 coefficients will be used."))
    }
    T_TT <- (MJD_TT - MJD_J2000)/36525 # Julian centuries in TT
    meanEclipticObliquity <- (23.439291 + T_TT * (-0.0130042 + T_TT * (-1.64e-7 + T_TT * 5.04e-7))) * pi/180
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
    meanElongationMoonSun <- meanElongationMoonSun %% (2*pi)
    meanLongitudeAscendingNodeMoon <- meanLongitudeAscendingNodeMoon %% (2*pi)
    api <- asteRiskData::nut_IAU1980[1:numberCoefficients ,1] * meanAnomalyMoon +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,2] * meanAnomalySun +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,3] * muMoon +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,4] * meanElongationMoonSun +
        asteRiskData::nut_IAU1980[1:numberCoefficients ,5] * meanLongitudeAscendingNodeMoon
    deltaPsi <- sum((asteRiskData::nut_IAU1980[1:numberCoefficients ,7] +
                     asteRiskData::nut_IAU1980[1:numberCoefficients ,8] * T_TT) *
                    sin(api))
    deltaEps <- sum((asteRiskData::nut_IAU1980[1:numberCoefficients ,9] +
                         asteRiskData::nut_IAU1980[1:numberCoefficients ,10] * T_TT) *
                        cos(api))
    # Convert from tenths of milliarcseconds to radians
    deltaPsi <- deltaPsi * pi /(10000 * 3600 * 180)
    deltaEps <- deltaEps * pi /(10000 * 3600 * 180)
    return(list(meanEclipticObliquity = meanEclipticObliquity,
                deltaPsi = deltaPsi,
                deltaEps = deltaEps))
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