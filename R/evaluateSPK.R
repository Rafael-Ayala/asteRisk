evaluateType1SPKSegment <- function(segment, targetEpochs) {
    # For types 1 and 21
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    segmentData <- segment$segmentData
    originalEpochsOrder <- order(order(targetEpochs))
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    recordFinalEpochs <- sapply(segmentData, `[[`, "finalEpoch")
    checkEpochsInInterval2 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= recordFinalEpochs[length(recordFinalEpochs)]
    if(!checkEpochsInInterval1) {
        if(checkEpochsInInterval2) {
            message(strwrap("At least one target epoch out of interval covered by 
                            segment metadata start and final epoch, but actually 
                            included in one of the records. Proceeding with evaluation 
                            for all target epochs"))
        } else {
            stop("At least one target epoch out of interval covered by any record in the segment.")
        }
    }
    targetEpochsRecords <- sapply(sortedEpochs, function(x) {which.max(recordFinalEpochs > x)})
    targetEpochsRecords[targetEpochs == segmentEndEpoch] <- length(recordFinalEpochs)
    requiredRecords <- unique(targetEpochsRecords)
    #results <- vector(mode="list", length=length(sortedEpochs))
    results <- matrix(nrow=length(sortedEpochs), ncol=10)
    counter <- 1
    for(record in requiredRecords) {
        recordData <- segmentData[[record]]
        epochsCurrentRecord <- sortedEpochs[targetEpochsRecords==record]
        referenceEpoch <- recordData$referenceEpoch
        referencePosition <- recordData$referencePosition
        referenceVelocity <- recordData$referenceVelocity
        stepsizeVector <- recordData$stepsizeFunctionVector
        MDA <- recordData$MDA
        KQMAX1 <- recordData$maxIntegrationOrderP1
        MQ2 <- KQMAX1 - 2
        KS <- KQMAX1 - 1
        for(targetEpoch in epochsCurrentRecord) {
            h <- targetEpoch - referenceEpoch
            # etaValues <- numeric(MQ2)
            # gammaValues <- numeric(MQ2)
            etaValues <- numeric(KQMAX1)
            gammaValues <- numeric(KQMAX1)
            etaValues[1] <- h/stepsizeVector[1]
            gammaValues[1] <- etaValues[1]
            # for(i in 2:(MQ2)) {
            for(i in 2:(KQMAX1)) {
                etaValues[i] <- h/stepsizeVector[i]
                gammaValues[i] <- (stepsizeVector[i-1] + h)/stepsizeVector[i]
            }
            S <- c(0,0,0)
            for(i in (MQ2+1):2) {
                S <- S + MDA[i,]*gammaValues[i-1]
            }
            S <- S + MDA[1,]
            acceleration <- S
            wValues <- 1/(1:KQMAX1)
            for(k in KQMAX1:3) {
                j <- 0
                for(i in k:KQMAX1) {
                    j <- j+1
                    wValues[i] <- gammaValues[j]*wValues[i-1] - etaValues[j]*wValues[i]
                }
            }
            S <- colSums(MDA[(MQ2+1):1,]*wValues[(MQ2+2):2])
            position <- referencePosition + h*(referenceVelocity + h*S)
            wValuesVelocity <- wValues
            for(i in 2:(MQ2+1)) {
                wValuesVelocity[i] <- gammaValues[i-1] * wValuesVelocity[i-1] - etaValues[i-1] * wValuesVelocity[i]
            }
            S <- colSums(MDA[(MQ2+1):1,]*wValuesVelocity[(MQ2+1):1])
            velocity <- referenceVelocity + h * S
            # results[[counter]] <- list(position=position,
            #                            velocity=velocity,
            #                            acceleration=acceleration)
            results[counter,] <- c(targetEpoch, position, velocity, acceleration)
            counter <- counter +1
        }
    }
    results <- results[originalEpochsOrder, ]
    if(length(targetEpochs) == 1) {
        results <- matrix(results, ncol=length(results))
    }
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ",
                           "accelerationX", "accelerationY", "accelerationZ")
    return(results)
}

evaluateType2SPKSegment <- function(segment, targetEpochs) {
    # For types 2, 3 and 14
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    chebyshevDegree <- segment$segmentData$polynomialDegree
    nCoeffs <- chebyshevDegree + 1
    segmentData <- segment$segmentData$chebyshevCoefficients
    colnamesSegmentData <- colnames(segmentData)
    originalEpochsOrder <- order(order(targetEpochs))
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    recordInitialEpochs <- segmentData[,"initialEpoch"]
    recordFinalEpochs <- segmentData[,"midPoint"] + segmentData[,"intervalRadius"]
    # checkEpochsInInterval2 <- recordFinalEpochs[1] <= sortedEpochs[1] &&
    #     sortedEpochs[length(sortedEpochs)] <= recordFinalEpochs[length(recordFinalEpochs)]
    if(!checkEpochsInInterval1) {
        stop("At least one target epoch out of interval covered by any record in the segment.")
    }
    targetEpochsRecords <- sapply(sortedEpochs, function(x) {which.max(recordFinalEpochs > x)})
    targetEpochsRecords[targetEpochs == segmentEndEpoch] <- length(recordFinalEpochs)
    requiredRecords <- unique(targetEpochsRecords)
    results <- matrix(nrow=length(sortedEpochs), ncol=10)
    counter <- 1
    positionXCoeffsCols <- grep("positionXCoeff", colnamesSegmentData)
    positionYCoeffsCols <- grep("positionYCoeff", colnamesSegmentData)
    positionZCoeffsCols <- grep("positionZCoeff", colnamesSegmentData)
    if(segmentTypeCode != 2) {
        velocityXCoeffsCols <- grep("velocityXCoeff", colnamesSegmentData)
        velocityYCoeffsCols <- grep("velocityYCoeff", colnamesSegmentData)
        velocityZCoeffsCols <- grep("velocityZCoeff", colnamesSegmentData)
        for(record in requiredRecords) {
            recordData <- segmentData[record,]
            epochsCurrentRecord <- sortedEpochs[targetEpochsRecords==record]
            positionXCoeffs <- recordData[positionXCoeffsCols]
            positionYCoeffs <- recordData[positionYCoeffsCols]
            positionZCoeffs <- recordData[positionZCoeffsCols]
            velocityXCoeffs <- recordData[velocityXCoeffsCols]
            velocityYCoeffs <- recordData[velocityYCoeffsCols]
            velocityZCoeffs <- recordData[velocityZCoeffsCols]
            # intervalStart <- recordData[1]
            intervalMidPoint <- recordData[2]
            intervalRadius <- recordData[3]
            #intervalEnd <- recordData[2] + recordData[3]
            for(targetEpoch in epochsCurrentRecord) {
                positionX <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                    intervalRadius, positionXCoeffs, 0, 1, FALSE)
                positionY <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                    intervalRadius, positionYCoeffs, 0, 1, FALSE)
                positionZ <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                    intervalRadius, positionZCoeffs, 0, 1, FALSE)
                velAccelX <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                    intervalRadius, velocityXCoeffs, 1, 1, FALSE)
                velAccelY <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                    intervalRadius, velocityYCoeffs, 1, 1, FALSE)
                velAccelZ <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                    intervalRadius, velocityZCoeffs, 1, 1, FALSE)
                position <- c(positionX, positionY, positionZ)
                velAccel <- rbind(velAccelX, velAccelY, velAccelZ)
                velocity <- velAccel[,1]
                acceleration <- velAccel[,2]
                results[counter,] <- c(targetEpoch, position, velocity, acceleration)
                counter <- counter + 1
            }
        }
    } else {
        for(record in requiredRecords) {
            recordData <- segmentData[record,]
            epochsCurrentRecord <- sortedEpochs[targetEpochsRecords==record]
            positionXCoeffs <- recordData[positionXCoeffsCols]
            positionYCoeffs <- recordData[positionYCoeffsCols]
            positionZCoeffs <- recordData[positionZCoeffsCols]
            # intervalStart <- recordData[1]
            intervalMidPoint <- recordData[2]
            intervalRadius <- recordData[3]
            # intervalEnd <- recordData[2] + recordData[3]
            for(targetEpoch in epochsCurrentRecord) {
                x <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                            intervalRadius, positionXCoeffs, 2, 1, FALSE)
                y <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                            intervalRadius, positionYCoeffs, 2, 1, FALSE)
                z <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                            intervalRadius, positionZCoeffs, 2, 1, FALSE)
                xyz <- rbind(x, y, z)
                position <- xyz[,1]
                velocity <- xyz[,2]
                acceleration <- xyz[,3]
                results[counter,] <- c(targetEpoch, position, velocity, acceleration)
                counter <- counter + 1
            }
        }
    }
    results <- results[originalEpochsOrder, ]
    if(length(targetEpochs) == 1) {
        results <- matrix(results, ncol=length(results))
    }
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ",
                           "accelerationX", "accelerationY", "accelerationZ")
    return(results)
}

evaluateType5SPKSegment <- function(segment, targetEpochs) {
    # for type 5 only
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    segmentData <- segment$segmentData$stateVectors
    centralBodyGM <- segment$segmentData$centralBodyGM
    originalEpochsOrder <- order(order(targetEpochs))
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    recordEpochs <- segmentData[,"epoch"]
    # checkEpochsInInterval2 <- recordFinalEpochs[1] <= sortedEpochs[1] &&
    #     sortedEpochs[length(sortedEpochs)] <= recordFinalEpochs[length(recordFinalEpochs)]
    if(!checkEpochsInInterval1) {
        stop("At least one target epoch out of interval covered the segment.")
    }
    targetEpochsRecords <- sapply(sortedEpochs, function(x) {which.max(recordEpochs > x)})
    targetEpochsRecords[targetEpochs == segmentEndEpoch] <- length(recordEpochs)
    requiredRecords <- unique(targetEpochsRecords)
    results <- matrix(nrow=length(sortedEpochs), ncol=7)
    counter <- 1
    for(record in requiredRecords) {
        epochsCurrentRecord <- sortedEpochs[record==targetEpochsRecords]
        refState1 <- segmentData[record - 1, ]
        refState2 <- segmentData[record, ]
        propStates1 <- twoBody(centralBodyGM, refState1[2:4], refState1[5:7], epochsCurrentRecord - refState1[1])
        propStates2 <- twoBody(centralBodyGM, refState2[2:4], refState2[5:7], epochsCurrentRecord - refState2[1])
        weightArg <- (pi*(epochsCurrentRecord - refState1[1])) / (refState2[1] - refState1[1])
        weights <- 0.5 + 0.5 * cos(weightArg)
        weightDers <- -sin(weightArg)*0.5*pi/(refState2[1] - refState1[1])
        position <- weights * propStates1[, 2:4] + (1 - weights) * propStates2[, 2:4]
        velocity <- weights * propStates1[, 5:7] + (1 - weights) * propStates2[, 5:7] +
            weightDers * (propStates1[, 2:4] - propStates2[, 2:4])
        if(length(epochsCurrentRecord) == 1) {
            results[counter, ] <- c(
                epochsCurrentRecord, position, velocity 
            )
            counter <- counter + 1
        } else {
            results[counter:(counter + length(epochsCurrentRecord) - 1), ] <- cbind(
                epochsCurrentRecord, position, velocity 
            )
            counter <- counter + length(epochsCurrentRecord)
        }
    }
    results <- results[originalEpochsOrder, ]
    if(length(targetEpochs) == 1) {
        results <- matrix(results, ncol=length(results))
    }
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ")
    return(results)
}

evaluateType8SPKSegment <- function(segment, targetEpochs) {
    # For types 8, 9, 12 and 13
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    if(segmentTypeCode == 8 || segmentTypeCode == 9) {
        interpolationType <- "Lagrange"
        lagrangeDegree <- segment$segmentData$polynomialDegree
        windowSize <- lagrangeDegree + 1
    } else if(segmentTypeCode == 12 || segmentTypeCode == 13) {
        interpolationType <- "Hermite-joint"
        windowSize <- segment$segmentData$windowSize
    }
    segmentData <- segment$segmentData$stateVectors
    originalEpochsOrder <- order(order(targetEpochs))
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    # checkEpochsInInterval2 <- recordFinalEpochs[1] <= sortedEpochs[1] &&
    #     sortedEpochs[length(sortedEpochs)] <= recordFinalEpochs[length(recordFinalEpochs)]
    if(!checkEpochsInInterval1) {
        stop("At least one target epoch out of interval covered by any record in the segment.")
    }
    recordEpochs <- segmentData[,"epoch"]
    numberRecords <- length(recordEpochs)
    targetEpochsFirstAfterRecords <- sapply(sortedEpochs, function(x) {which.max(recordEpochs > x)})
    targetEpochsFirstAfterRecords[targetEpochs == segmentEndEpoch] <- length(recordEpochs)
    intervalFirstAfterRecords <- recordEpochs[targetEpochsFirstAfterRecords] - sortedEpochs
    if(windowSize %% 2 == 1) {
        targetEpochsLastBeforeRecords <- targetEpochsFirstAfterRecords - 1
        intervalLastBeforeRecords <- sortedEpochs - recordEpochs[targetEpochsLastBeforeRecords]
        whichCloser <- max.col(-cbind(intervalLastBeforeRecords, intervalFirstAfterRecords), ties.method = "first")
        targetEpochsClosestRecords <- cbind(targetEpochsLastBeforeRecords, targetEpochsFirstAfterRecords)[
            cbind(1:length(sortedEpochs), whichCloser)
        ]
        targetEpochsRecords <- targetEpochsClosestRecords
        requiredRecords <- unique(targetEpochsClosestRecords)
        requiredRecordsWindows <- matrix(nrow=length(requiredRecords), ncol=2)
        halfWindowSizeMinus1 <- (windowSize-1)/2
        for(i in 1:length(requiredRecords)) {
            centralRecord <- requiredRecords[i]
            if(centralRecord - halfWindowSizeMinus1 <= 1) {
                requiredRecordsWindows[i, ] <- c(1, windowSize)
            } else if (centralRecord + halfWindowSizeMinus1 >= numberRecords) {
                requiredRecordsWindows[i, ] <- c((numberRecords - windowSize + 1), numberRecords)
            } else {
                requiredRecordsWindows[i, ] <- c((centralRecord - halfWindowSizeMinus1), (centralRecord + halfWindowSizeMinus1))
            }
        }
    } else {
        targetEpochsRecords <- targetEpochsFirstAfterRecords
        requiredRecords <- unique(targetEpochsFirstAfterRecords)
        requiredRecordsWindows <- matrix(nrow=length(requiredRecords), ncol=2)
        halfWindowSize <- windowSize/2
        for(i in 1:length(requiredRecords)) {
            centralRecord <- requiredRecords[i]
            if(centralRecord - halfWindowSize <= 1) {
                requiredRecordsWindows[i, ] <- c(1, windowSize)
            } else if (centralRecord + halfWindowSize - 1 >= numberRecords) {
                requiredRecordsWindows[i, ] <- c((numberRecords - windowSize + 1), numberRecords)
            } else {
                requiredRecordsWindows[i, ] <- c((centralRecord - halfWindowSize), (centralRecord + halfWindowSize - 1))
            }
        }
    }
    results <- interpolateSPKTypes8_9_12_13_18_19(sortedEpochs, segmentData,
                                                  requiredRecords, requiredRecordsWindows,
                                                  targetEpochsRecords, interpolationType)
    results <- results[originalEpochsOrder, ]
    if(length(targetEpochs) == 1) {
        results <- matrix(results, ncol=length(results))
    }
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ",
                           "accelerationX", "accelerationY", "accelerationZ")
    return(results)
}

evaluateType10SPKSegment <- function(segment, targetEpochs) {
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    segmentData <- segment$segmentData$TLEs
    constants <- segment$segmentData$constants
    J2 <- constants$J2
    J3 <- constants$J3
    J4 <- constants$J4
    sqrtGM <- constants$sqrtGM
    earthRadius <- constants$earthRadius
    ae <- constants$distUnitsPerRadius
    highAltBound <- constants$highAltBound
    lowAltBound <- constants$lowAltBound
    k2 <- 0.5 * J2 * ae ^ 2
    A30 <- -J3 * ae ^ 3
    k4 <- (-3 / 8) * J4 * ae ^ 4
    ke <- sqrtGM # already in correct units
    q0 <- highAltBound / earthRadius + ae
    s0 <- lowAltBound / earthRadius + ae
    qzms2t <- ((highAltBound - lowAltBound)/earthRadius)^4
    originalEpochsOrder <- order(order(targetEpochs))
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    recordEpochs <- segmentData[,"epoch"]
    # checkEpochsInInterval2 <- recordFinalEpochs[1] <= sortedEpochs[1] &&
    #     sortedEpochs[length(sortedEpochs)] <= recordFinalEpochs[length(recordFinalEpochs)]
    targetEpochsBeforeFirstRec <- recordEpochs[1] >= sortedEpochs
    targetEpochsAfterLastRec <- recordEpochs[length(recordEpochs)] <= sortedEpochs
    targetsOutsideRangeFlag <- any(c(targetEpochsBeforeFirstRec, targetEpochsAfterLastRec))
    if(targetsOutsideRangeFlag) {
        warning(strwrap("At least one target epoch before epoch of first TLE or 
        after epoch of last TLE. The closest TLE will be evaluated for these 
        target epochs, but results might be unreliable for distant epochs"))
        firstRecord <- segmentData[1, ]
        lastRecord <- segmentData[nrow(segmentData), ]
        deepFlagFirstRecord <- 2*pi/firstRecord[["meanMotion"]] >= 225
        deepFlagLastRecord <- 2*pi/lastRecord[["meanMotion"]] >= 225
        epochsBeforeFirstRec <- sortedEpochs[targetEpochsBeforeFirstRec]
        epochsAfterLastRec <- sortedEpochs[targetEpochsAfterLastRec]
        resultsBeforeFirstRec <- matrix(nrow=length(epochsBeforeFirstRec), ncol=7)
        resultsAfterLastRec <- matrix(nrow=length(epochsAfterLastRec), ncol=7)
        sortedEpochsNoOut <- sortedEpochs[!(targetEpochsBeforeFirstRec | targetEpochsAfterLastRec)]
        counter <- 1
        for(targetTime in epochsBeforeFirstRec) {
            if(deepFlagFirstRecord) {
                propState <- sdp4_(n0=firstRecord[["meanMotion"]], e0=firstRecord[["eccentricity"]], 
                                   i0=firstRecord[["inclination"]], M0=firstRecord[["meanAnomaly"]], 
                                   omega0=firstRecord[["perigeeArgument"]], OMEGA0=firstRecord[["ascension"]], 
                                   Bstar=firstRecord[["Bstar"]], initialDateTime=
                                       as.POSIXct(TDBSecondsToUTCSeconds_J2000(firstRecord[["epoch"]]), origin="2000-01-01 12:00:00", tz="GMT"),
                                   targetTime=(targetTime - firstRecord[["epoch"]])/60,
                                   keplerAccuracy=10e-12, maxKeplerIterations=10,
                                   J2, k2, A30, k4, ke, q0, s0, earthRadius, ae, qzms2t)
            } else {
                propState <- sgp4_(n0=firstRecord[["meanMotion"]], e0=firstRecord[["eccentricity"]], 
                                   i0=firstRecord[["inclination"]], M0=firstRecord[["meanAnomaly"]], 
                                   omega0=firstRecord[["perigeeArgument"]], OMEGA0=firstRecord[["ascension"]], 
                                   Bstar=firstRecord[["Bstar"]],
                                   targetTime=(targetTime - firstRecord[["epoch"]])/60,
                                   keplerAccuracy=10e-12, maxKeplerIterations=10,
                                   J2=J2, k2=k2, A30=A30, k4=k4, ke=ke, q0=q0, 
                                   s0=s0, earthRadius=earthRadius, ae=ae, qzms2t=qzms2t)
            }
            gcrfState <- TEMEtoGCRF(propState$position*1000, propState$velocity*1000, SPICEAlgorithm = TRUE,
                                    ephemerisTime = targetTime)
            resultsBeforeFirstRec[counter, ] <- c(targetTime, gcrfState$position, gcrfState$velocity)
            counter <- counter + 1
        }
        counter <- 1
        for(targetTime in epochsAfterLastRec) {
            if(deepFlagLastRecord) {
                propState <- sdp4_(n0=lastRecord[["meanMotion"]], e0=lastRecord[["eccentricity"]], 
                                   i0=lastRecord[["inclination"]], M0=lastRecord[["meanAnomaly"]], 
                                   omega0=lastRecord[["perigeeArgument"]], OMEGA0=lastRecord[["ascension"]], 
                                   Bstar=lastRecord[["Bstar"]], initialDateTime=
                                       as.POSIXct(TDBSecondsToUTCSeconds_J2000(lastRecord[["epoch"]]), origin="2000-01-01 12:00:00", tz="GMT"),
                                   targetTime=(targetTime - lastRecord[["epoch"]])/60,
                                   keplerAccuracy=10e-12, maxKeplerIterations=10,
                                   J2, k2, A30, k4, ke, q0, s0, earthRadius, ae, qzms2t)
            } else {
                propState <- sgp4_(n0=lastRecord[["meanMotion"]], e0=lastRecord[["eccentricity"]], 
                                   i0=lastRecord[["inclination"]], M0=lastRecord[["meanAnomaly"]], 
                                   omega0=lastRecord[["perigeeArgument"]], OMEGA0=lastRecord[["ascension"]], 
                                   Bstar=lastRecord[["Bstar"]],
                                   targetTime=(targetTime - lastRecord[["epoch"]])/60,
                                   keplerAccuracy=10e-12, maxKeplerIterations=10,
                                   J2=J2, k2=k2, A30=A30, k4=k4, ke=ke, q0=q0, 
                                   s0=s0, earthRadius=earthRadius, ae=ae, qzms2t=qzms2t)
            }
            gcrfState <- TEMEtoGCRF(propState$position*1000, propState$velocity*1000, SPICEAlgorithm = TRUE,
                                    ephemerisTime = targetTime)
            resultsAfterLastRec[counter, ] <- c(targetTime, gcrfState$position, gcrfState$velocity)
            counter <- counter + 1
        }
    } else {
        sortedEpochsNoOut <- sortedEpochs
        resultsBeforeFirstRec <- NULL
        resultsAfterLastRec <- NULL
    }
    targetEpochsRecords <- sapply(sortedEpochsNoOut, function(x) {which.max(recordEpochs >= x)})
    requiredRecords <- unique(targetEpochsRecords)
    results <- matrix(nrow=length(sortedEpochsNoOut), ncol=7)
    counter <- 1
    for(record in requiredRecords) {
        epochsCurrentRecord <- sortedEpochsNoOut[record==targetEpochsRecords]
        refState1 <- segmentData[record - 1, ]
        refState2 <- segmentData[record, ]
        meanMotion1 <- refState1[["meanMotion"]]
        deepFlag1 <- 2*pi/meanMotion1 >= 225
        meanMotion2 <- refState2[["meanMotion"]]
        deepFlag2 <- 2*pi/meanMotion2 >= 225
        weightArg <- (pi*(epochsCurrentRecord - refState1[1])) / (refState2[1] - refState1[1])
        weights <- 0.5 + 0.5 * cos(weightArg)
        weightDers <- -sin(weightArg)*0.5*pi/(refState2[1] - refState1[1])
        for(i in 1:length(epochsCurrentRecord)) {
            targetTime <- epochsCurrentRecord[i]
            weight <- weights[i]
            weightDer <- weightDers[i]
            if(deepFlag1) {
                propStates1 <- sdp4_(n0=refState1[["meanMotion"]], e0=refState1[["eccentricity"]], 
                                     i0=refState1[["inclination"]], M0=refState1[["meanAnomaly"]], 
                                     omega0=refState1[["perigeeArgument"]], OMEGA0=refState1[["ascension"]], 
                                     Bstar=refState1[["Bstar"]], initialDateTime=
                                         as.POSIXct(TDBSecondsToUTCSeconds_J2000(refState1[["epoch"]]), origin="2000-01-01 12:00:00", tz="GMT"),
                                     targetTime=(targetTime - refState1[["epoch"]])/60,
                                     keplerAccuracy=10e-12, maxKeplerIterations=10,
                                     J2, k2, A30, k4, ke, q0, s0, earthRadius, ae, qzms2t)
            } else {
                propStates1 <- sgp4_(n0=refState1[["meanMotion"]], e0=refState1[["eccentricity"]], 
                                     i0=refState1[["inclination"]], M0=refState1[["meanAnomaly"]], 
                                     omega0=refState1[["perigeeArgument"]], OMEGA0=refState1[["ascension"]], 
                                     Bstar=refState1[["Bstar"]],
                                     targetTime=(targetTime - refState1[["epoch"]])/60,
                                     keplerAccuracy=10e-12, maxKeplerIterations=10,
                                     J2=J2, k2=k2, A30=A30, k4=k4, ke=ke, q0=q0, 
                                     s0=s0, earthRadius=earthRadius, ae=ae, qzms2t=qzms2t)
            }
            if(deepFlag2) {
                propStates2 <- sdp4_(n0=refState2[["meanMotion"]], e0=refState2[["eccentricity"]], 
                                     i0=refState2[["inclination"]], M0=refState2[["meanAnomaly"]], 
                                     omega0=refState2[["perigeeArgument"]], OMEGA0=refState2[["ascension"]], 
                                     Bstar=refState2[["Bstar"]], initialDateTime=
                                         as.POSIXct(TDBSecondsToUTCSeconds_J2000(refState2[["epoch"]]), origin="2000-01-01 12:00:00", tz="GMT"),
                                     targetTime=(targetTime - refState2[["epoch"]])/60,
                                     keplerAccuracy=10e-12, maxKeplerIterations=10,
                                     J2, k2, A30, k4, ke, q0, s0, earthRadius, ae, qzms2t)
            } else {
                propStates2 <- sgp4_(n0=refState2[["meanMotion"]], e0=refState2[["eccentricity"]], 
                                     i0=refState2[["inclination"]], M0=refState2[["meanAnomaly"]], 
                                     omega0=refState2[["perigeeArgument"]], OMEGA0=refState2[["ascension"]], 
                                     Bstar=refState2[["Bstar"]],
                                     targetTime=(targetTime - refState2[["epoch"]])/60,
                                     keplerAccuracy=10e-12, maxKeplerIterations=10,
                                     J2=J2, k2=k2, A30=A30, k4=k4, ke=ke, q0=q0, 
                                     s0=s0, earthRadius=earthRadius, ae=ae, qzms2t=qzms2t)
            }
            position <- weight * propStates1$position + (1 - weight) * propStates2$position
            velocity <- weight * propStates1$velocity + (1 - weight) * propStates2$velocity +
                weightDer * (propStates1$position - propStates2$position)
            gcrfState <- TEMEtoGCRF(position*1000, velocity*1000, SPICEAlgorithm = TRUE,
                                    ephemerisTime = targetTime)
            results[counter, ] <- c(targetTime, gcrfState$position, gcrfState$velocity)
            counter <- counter + 1
        }
    }
    results <- rbind(resultsBeforeFirstRec, results, resultsAfterLastRec)
    results <- results[originalEpochsOrder, ]
    if(length(targetEpochs) == 1) {
        results <- matrix(results, ncol=length(results))
    }
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ")
    return(results)
}

evaluateType15SPKSegment <- function(segment, targetEpochs) {
    # for type 15 only. Seems like type 15 segments have a single record
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    if(any(targetEpochs < segmentStartEpoch | targetEpochs > segmentEndEpoch)) {
        stop("At least one target epoch out of interval covered the segment.")
    }
    segmentData <- segment$segmentData
    epochPeriapsis <- segmentData$epochPeriapsis
    trajectoryPoleVector <- c(segmentData$unitVectorTrajectoryPoleX, segmentData$unitVectorTrajectoryPoleY,
                    segmentData$unitVectorTrajectoryPoleZ)
    periapsisVector <- c(segmentData$unitVectorPeriapsisX, segmentData$unitVectorPeriapsisY,
                         segmentData$unitVectorPeriapsisZ)
    p <- segmentData$semiLatusRectum
    e <- segmentData$eccentricity
    J2Flag <- segmentData$J2ProcessingFlag
    centralBodyPoleVector <- c(segmentData$unitVectorCentralBodyPoleX, segmentData$unitVectorCentralBodyPoleY,
                               segmentData$unitVectorCentralBodyPoleZ)
    GM <- segmentData$centralBodyGM
    J2 <- segmentData$centralBodyJ2
    centralBodyRadius <- segmentData$centralBodyRadius # rpl in SPICE
    ## checks block
    if(p <= 0) {
        stop("Non-positive semi latus rectum")
    }
    if(e < 0) {
        stop("Negative eccentricity")
    }
    if(GM <= 0) {
        stop("Non-positive GM")
    }
    if(centralBodyRadius < 0) {
        stop("Negative central body radius")
    }
    trajectoryPoleVectorMod <- sqrt(sum(trajectoryPoleVector^2))
    if(trajectoryPoleVectorMod == 0) {
        stop("0-length trajectory pole vector")
    }
    periapsisVectorMod <- sqrt(sum(periapsisVector^2))
    if(periapsisVectorMod == 0) {
        stop("0-length periapsis vector")
    }
    centralBodyPoleVectorMod <- sqrt(sum(centralBodyPoleVector^2))
    if(centralBodyPoleVectorMod == 0) {
        stop("0-length central body pole vector")
    }
    trajectoryPoleVector <- trajectoryPoleVector/trajectoryPoleVectorMod #tp
    periapsisVector <- periapsisVector/periapsisVectorMod #pa
    centralBodyPoleVector <- centralBodyPoleVector/centralBodyPoleVectorMod #pv
    dotProductTrajectoryPeriapsis <- sum(periapsisVector * trajectoryPoleVector)
    if(abs(dotProductTrajectoryPeriapsis) > 1e-5) {
        stop("Trajectory pole and periapsis vectors not orthogonal")
    }
    periapsisDist <- p/(e+1) # called near__ in SPICE
    periapsisSpeed <- sqrt(GM/p) * (e+1)
    periapsisPosition <- periapsisDist * periapsisVector
    periapsisVelocity <- vectorCrossProduct3D(trajectoryPoleVector,
                                              periapsisVector) * periapsisSpeed
    results <- matrix(nrow=length(targetEpochs), ncol=7)
    diffTimes <- targetEpochs - epochPeriapsis
    numberEpochs <- length(diffTimes)
    twoBodyResults <- twoBody(GM, periapsisPosition, periapsisVelocity, diffTimes)
    if(J2Flag != 3 && J2 != 0 && e < 1 && periapsisDist > centralBodyRadius) {
        oneMe2 <- 1 - e*e
        meanAnomRate <- oneMe2/p * sqrt(GM*oneMe2/p) # dmdt in SPICE
        meanAnomaly <- meanAnomRate * diffTimes
        theta <- rem(meanAnomaly, 2*pi)
        thetaValuesOverPi <- theta > pi
        thetaValuesBelowNegPi <- theta < -pi
        theta[thetaValuesOverPi] <- theta[thetaValuesOverPi] - 2*pi
        theta[thetaValuesBelowNegPi] <- theta[thetaValuesBelowNegPi] + 2*pi
        # if(theta > pi) {
        #     theta <- theta - 2*pi
        # } else if(theta < -pi) {
        #     theta <- theta + 2*pi
        # }
        k2pi <- meanAnomaly - theta
        trueAnomaly <- numeric(numberEpochs)
        for(i in 1:numberEpochs) {
            trueAnomaly[i] <- angularSeparation(periapsisVector, twoBodyResults[i, 2:4])
        }
        trueAnomaly <- trueAnomaly*sign(theta)
        trueAnomaly <- trueAnomaly + k2pi
        cosinc <- sum(centralBodyPoleVector * trajectoryPoleVector)
        ratio <- centralBodyRadius/p
        z <- trueAnomaly * 1.5 * J2 * ratio * ratio
        dnode <- - z * cosinc
        dperi <- z * (cosinc * cosinc * 2.5 - 0.5)
        if(J2Flag == 1) {
            # regression of the line of nodes only
            for(i in 1:numberEpochs) {
                results[i, 2:4] <- rotate3DVectorAngleAxis(twoBodyResults[i, 2:4], 
                                                           centralBodyPoleVector,
                                                           dnode[i])
                results[i, 5:7] <- rotate3DVectorAngleAxis(twoBodyResults[i, 5:7], 
                                                           centralBodyPoleVector,
                                                           dnode[i])
            }
        } else if(J2Flag == 2) {
            # precession of the line of apsides only
            for(i in 1:numberEpochs) {
                results[i, 2:4] <- rotate3DVectorAngleAxis(twoBodyResults[i, 2:4], 
                                                           trajectoryPoleVector,
                                                           dperi[i])
                results[i, 5:7] <- rotate3DVectorAngleAxis(twoBodyResults[i, 5:7], 
                                                           trajectoryPoleVector,
                                                           dperi[i])
            }
        } else {
            # all corrections
            for(i in 1:numberEpochs) {
                results[i, 2:4] <- rotate3DVectorAngleAxis(rotate3DVectorAngleAxis(twoBodyResults[i, 2:4], 
                                                                                   trajectoryPoleVector,
                                                                                   dperi[i]), 
                                                           centralBodyPoleVector,
                                                           dnode[i])
                results[i, 5:7] <- rotate3DVectorAngleAxis(rotate3DVectorAngleAxis(twoBodyResults[i, 5:7], 
                                                                                   trajectoryPoleVector,
                                                                                   dperi[i]), 
                                                           centralBodyPoleVector,
                                                           dnode[i])
            }
        }
    } else {
        results <- twoBodyResults
    }
    results[, 1] <- targetEpochs
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ")
    return(results)
}

evaluateType17SPKSegment <- function(segment, targetEpochs) {
    # for type 17 only. Seems like type 17 segments have a single record also
    # OUTPUT HAS BEEN VERIFIED TO MATCH CSPICE'S
    # but we dont enforce eccentricity < 0.9. Instead give a warning
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    if(any(targetEpochs < segmentStartEpoch | targetEpochs > segmentEndEpoch)) {
        stop("At least one target epoch out of interval covered the segment.")
    }
    segmentData <- segment$segmentData
    epochPeriapsis <- segmentData$epochPeriapsis
    semiMajorAxis <- segmentData$semiMajorAxis
    H0 <- segmentData$equinoctialH
    K0 <- segmentData$equinoctialK
    s <- sqrt(H0 * H0 + K0 * K0)
    meanLongitude0 <- segmentData$meanLongitude
    P0 <- segmentData$equinoctialP
    Q0 <- segmentData$equinoctialQ
    longitudePeriapsisRate <- segmentData$longitudePeriapsisDerivative
    meanLongitudeRate <- segmentData$meanLongitudeDerivative
    longitudeAscendingNodeRate <- segmentData$longitudeAscendingNodeDerivative
    poleRightAscension <- segmentData$equatorialPoleRightAscension
    poleDeclination <- segmentData$equatorialPoleDeclination
    if(semiMajorAxis <= 0) {
        stop("Non-positive semi major axis")
    }
    diffTimes <- targetEpochs - epochPeriapsis
    results <- evaluateEquinoctialElements(semiMajorAxis, H0, K0, meanLongitude0,
                                           P0, Q0, longitudePeriapsisRate, meanLongitudeRate,
                                           longitudeAscendingNodeRate, poleRightAscension,
                                           poleDeclination, diffTimes)
    results[, 1] <- targetEpochs
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ")
    return(results)
}

evaluateType18SPKSegment <- function(segment, targetEpochs) {
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    segmentSubTypeCode <- segment$segmentData$subTypeCode
    interpolationType <- segment$segmentData$interpolationType
    windowSize <- segment$segmentData$windowSize
    polynomialDegree <- segment$segmentData$polynomialDegree
    segmentData <- segment$segmentData$stateVectors
    originalEpochsOrder <- order(order(targetEpochs))
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    # checkEpochsInInterval2 <- recordFinalEpochs[1] <= sortedEpochs[1] &&
    #     sortedEpochs[length(sortedEpochs)] <= recordFinalEpochs[length(recordFinalEpochs)]
    if(!checkEpochsInInterval1) {
        stop("At least one target epoch out of interval covered by any record in the segment.")
    }
    recordEpochs <- segmentData[,"epoch"]
    numberRecords <- length(recordEpochs)
    targetEpochsFirstAfterRecords <- sapply(sortedEpochs, function(x) {which.max(recordEpochs > x)})
    targetEpochsFirstAfterRecords[targetEpochs == segmentEndEpoch] <- length(recordEpochs)
    intervalFirstAfterRecords <- recordEpochs[targetEpochsFirstAfterRecords] - sortedEpochs
    if(windowSize %% 2 == 1) {
        targetEpochsLastBeforeRecords <- targetEpochsFirstAfterRecords - 1
        intervalLastBeforeRecords <- sortedEpochs - recordEpochs[targetEpochsLastBeforeRecords]
        whichCloser <- max.col(-cbind(intervalLastBeforeRecords, intervalFirstAfterRecords), ties.method = "first")
        targetEpochsClosestRecords <- cbind(targetEpochsLastBeforeRecords, targetEpochsFirstAfterRecords)[
            cbind(1:length(sortedEpochs), whichCloser)
        ]
        targetEpochsRecords <- targetEpochsClosestRecords
        requiredRecords <- unique(targetEpochsClosestRecords)
        requiredRecordsWindows <- matrix(nrow=length(requiredRecords), ncol=2)
        halfWindowSizeMinus1 <- (windowSize-1)/2
        for(i in 1:length(requiredRecords)) {
            centralRecord <- requiredRecords[i]
            firstRecord <- centralRecord - halfWindowSizeMinus1
            lastRecord <- centralRecord + halfWindowSizeMinus1
            if(firstRecord < 1) firstRecord <- 1
            if(lastRecord > numberRecords) lastRecord <- numberRecords
            requiredRecordsWindows[i, ] <- c(firstRecord, lastRecord)
        }
    } else {
        targetEpochsRecords <- targetEpochsFirstAfterRecords
        requiredRecords <- unique(targetEpochsFirstAfterRecords)
        requiredRecordsWindows <- matrix(nrow=length(requiredRecords), ncol=2)
        halfWindowSize <- windowSize/2
        for(i in 1:length(requiredRecords)) {
            centralRecord <- requiredRecords[i]
            firstRecord <- centralRecord - halfWindowSize
            lastRecord <- centralRecord + halfWindowSize - 1
            if(firstRecord < 1) firstRecord <- 1
            if(lastRecord > numberRecords) lastRecord <- numberRecords
            requiredRecordsWindows[i, ] <- c(firstRecord, lastRecord)
        }
    }
    results <- interpolateSPKTypes8_9_12_13_18_19(sortedEpochs, segmentData,
                                                  requiredRecords, requiredRecordsWindows,
                                                  targetEpochsRecords, interpolationType)
    results <- results[originalEpochsOrder, ]
    if(length(targetEpochs) == 1) {
        results <- matrix(results, ncol=length(results))
    }
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ",
                           "accelerationX", "accelerationY", "accelerationZ")
    return(results)
}

evaluateType19SPKSegment <- function(segment, targetEpochs) {
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    originalEpochsOrder <- order(order(targetEpochs))
    boundaryChoiceFlag <- segment$segmentData$boundaryChoiceFlag
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    # checkEpochsInInterval2 <- recordFinalEpochs[1] <= sortedEpochs[1] &&
    #     sortedEpochs[length(sortedEpochs)] <= recordFinalEpochs[length(recordFinalEpochs)]
    if(!checkEpochsInInterval1) {
        stop("At least one target epoch out of interval covered by any record in the segment.")
    }
    minisegments <- segment$segmentData$minisegments
    interpolationIntervalStarts <- sapply(minisegments, `[[`, "intervalStartEpoch")
    interpolationIntervalEnds <- sapply(minisegments, `[[`, "intervalEndEpoch")
    if(boundaryChoiceFlag) {
        targetEpochsRequiredMinisegment <- sapply(sortedEpochs, function(x) {which.max(interpolationIntervalEnds > x)})
    } else {
        targetEpochsRequiredMinisegment <- sapply(sortedEpochs, function(x) {which.max(interpolationIntervalEnds >= x)})
    }
    requiredMinisegments <- unique(targetEpochsRequiredMinisegment)
    results <- matrix(nrow=length(sortedEpochs), ncol=10)
    counter <- 1
    for(minisegmentIndex in requiredMinisegments) {
        currentMinisgment <- minisegments[[minisegmentIndex]]
        epochsCurrentMinisegment <- sortedEpochs[targetEpochsRequiredMinisegment==minisegmentIndex]
        subTypeCode <- currentMinisgment$subTypeCode
        interpolationType <- currentMinisgment$interpolationType
        windowSize <- currentMinisgment$windowSize
        polynomialDegree <- currentMinisgment$polynomialDegree
        minisegmentData <- currentMinisgment$stateVectors
        recordEpochs <- minisegmentData[,"epoch"]
        numberRecords <- length(recordEpochs)
        targetEpochsFirstAfterRecords <- sapply(epochsCurrentMinisegment, function(x) {which.max(recordEpochs > x)})
        targetEpochsFirstAfterRecords[targetEpochs == segmentEndEpoch] <- length(recordEpochs)
        intervalFirstAfterRecords <- recordEpochs[targetEpochsFirstAfterRecords] - epochsCurrentMinisegment
        if(windowSize %% 2 == 1) {
            targetEpochsLastBeforeRecords <- targetEpochsFirstAfterRecords - 1
            intervalLastBeforeRecords <- epochsCurrentMinisegment - recordEpochs[targetEpochsLastBeforeRecords]
            whichCloser <- max.col(-cbind(intervalLastBeforeRecords, intervalFirstAfterRecords), ties.method = "first")
            targetEpochsClosestRecords <- cbind(targetEpochsLastBeforeRecords, targetEpochsFirstAfterRecords)[
                cbind(1:length(epochsCurrentMinisegment), whichCloser)
            ]
            targetEpochsRecords <- targetEpochsClosestRecords
            requiredRecords <- unique(targetEpochsClosestRecords)
            requiredRecordsWindows <- matrix(nrow=length(requiredRecords), ncol=2)
            halfWindowSizeMinus1 <- (windowSize-1)/2
            for(i in 1:length(requiredRecords)) {
                centralRecord <- requiredRecords[i]
                firstRecord <- centralRecord - halfWindowSizeMinus1
                lastRecord <- centralRecord + halfWindowSizeMinus1
                if(firstRecord < 1) firstRecord <- 1
                if(lastRecord > numberRecords) lastRecord <- numberRecords
                requiredRecordsWindows[i, ] <- c(firstRecord, lastRecord)
            }
        } else {
            targetEpochsRecords <- targetEpochsFirstAfterRecords
            requiredRecords <- unique(targetEpochsFirstAfterRecords)
            requiredRecordsWindows <- matrix(nrow=length(requiredRecords), ncol=2)
            halfWindowSize <- windowSize/2
            for(i in 1:length(requiredRecords)) {
                centralRecord <- requiredRecords[i]
                firstRecord <- centralRecord - halfWindowSize
                lastRecord <- centralRecord + halfWindowSize - 1
                if(firstRecord < 1) firstRecord <- 1
                if(lastRecord > numberRecords) lastRecord <- numberRecords
                requiredRecordsWindows[i, ] <- c(firstRecord, lastRecord)
            }
        }
        resultsCurrentMinisegment <- interpolateSPKTypes8_9_12_13_18_19(epochsCurrentMinisegment, minisegmentData,
                                                                        requiredRecords, requiredRecordsWindows,
                                                                        targetEpochsRecords, interpolationType)
        results[counter:(counter + length(epochsCurrentMinisegment) - 1), ] <- resultsCurrentMinisegment
        counter <- counter + length(epochsCurrentMinisegment)
    }
    results <- results[originalEpochsOrder, ]
    if(length(targetEpochs) == 1) {
        results <- matrix(results, ncol=length(results))
    }
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ",
                           "accelerationX", "accelerationY", "accelerationZ")
    return(results)
}

evaluateType20SPKSegment <- function(segment, targetEpochs) {
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    chebyshevDegree <- segment$segmentData$polynomialDegree
    dScale <- segment$segmentData$dScale
    tScale <- segment$segmentData$tScale
    nCoeffs <- chebyshevDegree + 1
    segmentData <- segment$segmentData$chebyshevCoefficients
    colnamesSegmentData <- colnames(segmentData)
    originalEpochsOrder <- order(order(targetEpochs))
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    recordInitialEpochs <- segmentData[,"initialEpoch"]
    recordFinalEpochs <- segmentData[,"midPoint"] + segmentData[,"intervalRadius"]
    if(!checkEpochsInInterval1) {
        stop("At least one target epoch out of interval covered by any record in the segment.")
    }
    targetEpochsRecords <- sapply(sortedEpochs, function(x) {which.max(recordFinalEpochs > x)})
    targetEpochsRecords[targetEpochs == segmentEndEpoch] <- length(recordFinalEpochs)
    requiredRecords <- unique(targetEpochsRecords)
    results <- matrix(nrow=length(sortedEpochs), ncol=10)
    counter <- 1
    velocityXCoeffsCols <- grep("velocityXCoeff", colnamesSegmentData)
    velocityYCoeffsCols <- grep("velocityYCoeff", colnamesSegmentData)
    velocityZCoeffsCols <- grep("velocityZCoeff", colnamesSegmentData)
    for(record in requiredRecords) {
        recordData <- segmentData[record,]
        epochsCurrentRecord <- sortedEpochs[targetEpochsRecords==record]
        velocityXCoeffs <- recordData[velocityXCoeffsCols]
        velocityYCoeffs <- recordData[velocityYCoeffsCols]
        velocityZCoeffs <- recordData[velocityZCoeffsCols]
        midPointPositionX <- recordData["midPointPositionX"]
        midPointPositionY <- recordData["midPointPositionY"]
        midPointPositionZ <- recordData["midPointPositionZ"]
        midPointPosition <- c(midPointPositionX, midPointPositionY, midPointPositionZ) * dScale
        intervalMidPoint <- recordData[2]
        intervalRadius <- recordData[3]
        for(targetEpoch in epochsCurrentRecord) {
            posVelAccelX <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                intervalRadius, velocityXCoeffs, 1, tScale, TRUE) * dScale
            posVelAccelY <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                   intervalRadius, velocityYCoeffs, 1, tScale, TRUE) * dScale
            posVelAccelZ <- clenshawAllDerivatives(targetEpoch, nCoeffs, intervalMidPoint,
                                                   intervalRadius, velocityZCoeffs, 1, tScale, TRUE) * dScale
            posVelAccel <- rbind(posVelAccelX, posVelAccelY, posVelAccelZ)
            position <- (posVelAccel[,1] + midPointPosition) 
            velocity <- posVelAccel[,2]/tScale
            acceleration <- posVelAccel[,3]
            results[counter,] <- c(targetEpoch, position, velocity, acceleration)
            counter <- counter + 1
        }
    }
    results <- results[originalEpochsOrder, ]
    if(length(targetEpochs) == 1) {
        results <- matrix(results, ncol=length(results))
    }
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ",
                           "accelerationX", "accelerationY", "accelerationZ")
    return(results)
}

evaluateSPKSegment <- function(segment, targetEpochs) {
    segmentTypeCode <- segment$segmentSummary$SPKTypeCode
    if(segmentTypeCode == 1 ||
       segmentTypeCode == 21) {
        results <- evaluateType1SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 2 ||
              segmentTypeCode == 3 ||
              segmentTypeCode == 14) {
        results <- evaluateType2SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 5) {
        results <- evaluateType5SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 8 ||
              segmentTypeCode == 9 ||
              segmentTypeCode == 12 ||
              segmentTypeCode == 13) {
        results <- evaluateType8SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 10) {
        results <- evaluateType10SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 15) {
        results <- evaluateType15SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 17) {
        results <- evaluateType17SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 18) {
        results <- evaluateType18SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 19) {
        results <- evaluateType19SPKSegment(segment, targetEpochs)
    } else if(segmentTypeCode == 20) {
        results <- evaluateType20SPKSegment(segment, targetEpochs)
    }
    return(results)
}