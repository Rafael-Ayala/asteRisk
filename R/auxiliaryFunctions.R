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

rawToInt <- function(x) {
    readBin(x, integer())
}

rawToDouble <- function(x) { 
    readBin(x, numeric())
}

interpolateSPKTypes8_9_12_13_18_19 <- function(sortedEpochs, segmentData,
                                               requiredRecords, requiredRecordsWindows,
                                               targetEpochsRecords, interpolationType) {
    results <- matrix(nrow=length(sortedEpochs), ncol=10)
    counter <- 1
    if(interpolationType == "Lagrange") {
        for(i in 1:length(requiredRecords)) {
            epochsCurrentRecord <- sortedEpochs[targetEpochsRecords==requiredRecords[i]]
            requiredRecordsIndexes <- requiredRecordsWindows[i, 1]:requiredRecordsWindows[i, 2]
            allRequiredRecords <- segmentData[requiredRecordsIndexes, ]
            refEpochs <- allRequiredRecords[, 1]
            xPositions0 <- allRequiredRecords[, 2]
            yPositions0 <- allRequiredRecords[, 3]
            zPositions0 <- allRequiredRecords[, 4]
            xVelocities0 <- allRequiredRecords[, 5]
            yVelocities0 <- allRequiredRecords[, 6]
            zVelocities0 <- allRequiredRecords[, 7]
            lagrangeXPos <- lagrange_(refEpochs, xPositions0)
            lagrangeYPos <- lagrange_(refEpochs, yPositions0)
            lagrangeZPos <- lagrange_(refEpochs, zPositions0)
            lagrangeXVel <- lagrange_(refEpochs, xVelocities0)
            lagrangeYVel <- lagrange_(refEpochs, yVelocities0)
            lagrangeZVel <- lagrange_(refEpochs, zVelocities0)
            xPositions <- lagrangeXPos(epochsCurrentRecord, 0)
            yPositions <- lagrangeYPos(epochsCurrentRecord, 0)
            zPositions <- lagrangeZPos(epochsCurrentRecord, 0)
            xVelAccel <- lagrangeXVel(epochsCurrentRecord, 1)
            yVelAccel <- lagrangeYVel(epochsCurrentRecord, 1)
            zVelAccel <- lagrangeZVel(epochsCurrentRecord, 1)
            xVelocities <- xVelAccel[1,]
            yVelocities <- yVelAccel[1,]
            zVelocities <- zVelAccel[1,]
            xAccelerations <- xVelAccel[2,]
            yAccelerations <- yVelAccel[2,]
            zAccelerations <- zVelAccel[2,]
            # newResultsBlock <- cbind(epochsCurrentRecord, xPositions, yPositions, zPositions,
            #                          xVelocities, yVelocities, zVelocities, xAcceleration,
            #                          yAcceleration, zAcceleration)
            if(length(epochsCurrentRecord) == 1) {
                results[counter, ] <- c(
                    epochsCurrentRecord, xPositions, yPositions, zPositions, 
                    xVelocities, yVelocities, zVelocities, xAccelerations,
                    yAccelerations, zAccelerations
                )
                counter <- counter + 1
            } else {
                results[counter:(counter + length(epochsCurrentRecord) - 1), ] <- cbind(
                    epochsCurrentRecord, xPositions, yPositions, zPositions,
                    xVelocities, yVelocities, zVelocities, xAccelerations,
                    yAccelerations, zAccelerations
                )
                counter <- counter + length(epochsCurrentRecord)
            }
        }
    } else if(interpolationType == "Hermite-joint") {
        for(i in 1:length(requiredRecords)) {
            epochsCurrentRecord <- sortedEpochs[targetEpochsRecords==requiredRecords[i]]
            allRequiredRecords <- segmentData[requiredRecordsWindows[i, ], ]
            refEpochs <- allRequiredRecords[, 1]
            xPositions0 <- allRequiredRecords[, 2]
            yPositions0 <- allRequiredRecords[, 3]
            zPositions0 <- allRequiredRecords[, 4]
            xVelocities0 <- allRequiredRecords[, 5]
            yVelocities0 <- allRequiredRecords[, 6]
            zVelocities0 <- allRequiredRecords[, 7]
            hermiteX <- hermite_(refEpochs, xPositions0, xVelocities0)
            hermiteY <- hermite_(refEpochs, yPositions0, yVelocities0)
            hermiteZ <- hermite_(refEpochs, zPositions0, zVelocities0)
            x <- hermiteX(epochsCurrentRecord, 2)
            y <- hermiteY(epochsCurrentRecord, 2)
            z <- hermiteZ(epochsCurrentRecord, 2)
            xPositions <- x[1, ]
            yPositions <- y[1, ]
            zPositions <- z[1, ]
            xVelocities <- x[2,]
            yVelocities <- y[2,]
            zVelocities <- z[2,]
            xAccelerations <- x[3,]
            yAccelerations <- y[3,]
            zAccelerations <- z[3,]
            # newResultsBlock <- cbind(epochsCurrentRecord, xPositions, yPositions, zPositions,
            #                          xVelocities, yVelocities, zVelocities, xAcceleration,
            #                          yAcceleration, zAcceleration)
            if(length(epochsCurrentRecord) == 1) {
                results[counter, ] <- c(
                    epochsCurrentRecord, xPositions, yPositions, zPositions, 
                    xVelocities, yVelocities, zVelocities, xAccelerations,
                    yAccelerations, zAccelerations
                )
                counter <- counter + 1
            } else {
                results[counter:(counter + length(epochsCurrentRecord) - 1), ] <- cbind(
                    epochsCurrentRecord, xPositions, yPositions, zPositions,
                    xVelocities, yVelocities, zVelocities, xAccelerations,
                    yAccelerations, zAccelerations
                )
                counter <- counter + length(epochsCurrentRecord)
            }
        }
    } else if(interpolationType == "Hermite") {
        for(i in 1:length(requiredRecords)) {
            epochsCurrentRecord <- sortedEpochs[targetEpochsRecords==requiredRecords[i]]
            allRequiredRecords <- segmentData[requiredRecordsWindows[i, ], ]
            refEpochs <- allRequiredRecords[, 1]
            xPositions0 <- allRequiredRecords[, 2]
            yPositions0 <- allRequiredRecords[, 3]
            zPositions0 <- allRequiredRecords[, 4]
            xVelocities0_1 <- allRequiredRecords[, 5]
            yVelocities0_1 <- allRequiredRecords[, 6]
            zVelocities0_1 <- allRequiredRecords[, 7]
            xVelocities0_2 <- allRequiredRecords[, 8]
            yVelocities0_2 <- allRequiredRecords[, 9]
            zVelocities0_2 <- allRequiredRecords[, 10]
            xAccelerations0 <- allRequiredRecords[, 11]
            yAccelerations0 <- allRequiredRecords[, 12]
            zAccelerations0 <- allRequiredRecords[, 13]
            hermitePosVelX <- hermite_(refEpochs, xPositions0, xVelocities0_1)
            hermitePosVelY <- hermite_(refEpochs, yPositions0, yVelocities0_1)
            hermitePosVelZ <- hermite_(refEpochs, zPositions0, zVelocities0_1)
            hermiteVelAccelX <- hermite_(refEpochs, xVelocities0_2, xAccelerations0)
            hermiteVelAccelY <- hermite_(refEpochs, yVelocities0_2, yAccelerations0)
            hermiteVelAccelZ <- hermite_(refEpochs, zVelocities0_2, zAccelerations0)
            xPositions <- hermitePosVelX(epochsCurrentRecord, 0)
            yPositions <- hermitePosVelX(epochsCurrentRecord, 0)
            zPositions <- hermitePosVelX(epochsCurrentRecord, 0)
            xVelAccel <- hermiteVelAccelX(epochsCurrentRecord, 1)
            yVelAccel <- hermiteVelAccelY(epochsCurrentRecord, 1)
            zVelAccel <- hermiteVelAccelZ(epochsCurrentRecord, 1)
            xVelocities <- xVelAccel[1,]
            yVelocities <- yVelAccel[1,]
            zVelocities <- zVelAccel[1,]
            xAccelerations <- xVelAccel[2,]
            yAccelerations <- yVelAccel[2,]
            zAccelerations <- zVelAccel[2,]
            if(length(epochsCurrentRecord) == 1) {
                results[counter, ] <- c(
                    epochsCurrentRecord, xPositions, yPositions, zPositions, 
                    xVelocities, yVelocities, zVelocities, xAccelerations,
                    yAccelerations, zAccelerations
                )
                counter <- counter + 1
            } else {
                results[counter:(counter + length(epochsCurrentRecord) - 1), ] <- cbind(
                    epochsCurrentRecord, xPositions, yPositions, zPositions,
                    xVelocities, yVelocities, zVelocities, xAccelerations,
                    yAccelerations, zAccelerations
                )
                counter <- counter + length(epochsCurrentRecord)
            }
        }
    }
    return(results)
}

vectorCrossProductDerivative3D <- function(position1, velocity1, position2, velocity2) {
    # equivalent to SPICE's dvcross
    crossProd <- vectorCrossProduct3D(position1, position2)
    crossProdRate <- vectorCrossProduct3D(velocity1, position2) + vectorCrossProduct3D(position1, velocity2)
    return(list(
        crossProduct = crossProd,
        crossProductDerivative = crossProdRate
    ))
}

vectorCrossProductDerivative3D_unitNorm <- function(position1, velocity1, position2, velocity2) {
    # equivalent to SPICE's ducross
    scale1 <- max(position1)
    scale2 <- max(position2)
    position1 <- position1/scale1
    velocity1 <- velocity1/scale1
    position2 <- position2/scale2
    velocity2 <- velocity2/scale2
    unnormResults <- vectorCrossProductDerivative3D(position1, velocity1, position2, velocity2)
    results <- dvhat(unnormResults$crossProduct, unnormResults$crossProductDerivative)
    return(list(
        normCrossProduct = results[1:3],
        normCrossProductDerivative = results[4:6]
    ))
}


