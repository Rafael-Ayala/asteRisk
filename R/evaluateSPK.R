evaluateType1SPKSegment <- function(segment, targetEpochs) {
    segmentStartEpoch <- segment$segmentSummary$initialEpoch
    segmentEndEpoch <- segment$segmentSummary$finalEpoch
    segmentData <- segment$segmentData
    originalEpochsOrder <- order(targetEpochs)
    sortedEpochs <- sort(targetEpochs)
    checkEpochsInInterval1 <- segmentStartEpoch <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= segmentEndEpoch
    recordFinalEpochs <- sapply(segmentData, `[[`, "finalEpoch")
    checkEpochsInInterval2 <- recordFinalEpochs[1] <= sortedEpochs[1] &&
        sortedEpochs[length(sortedEpochs)] <= recordFinalEpochs[length(recordFinalEpochs)]
    if(!checkEpochsInInterval1) {
        if(checkEpochsInInterval2) {
            warning(strwrap("At least one target epoch out of interval covered by segment metadata
                            start and final epoch, but actually included in one of the records.
                            Proceeding with evaluation for all target epochs"))
        } else {
            stop("At least one target epoch out of interval covered by any record in the segment.")
        }
    }
    recordFinalEpochs <- sapply(segmentData, `[[`, "finalEpoch")
    targetEpochsRecords <- sapply(sortedEpochs, function(x) {which.max(recordFinalEpochs > x)})
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
    # results <- results[originalEpochsOrder]
    results <- results[originalEpochsOrder, ]
    colnames(results) <- c("epoch", "positionX", "positionY", "positionZ",
                           "velocityX", "velocityY", "velocityZ",
                           "accelerationX", "accelerationY", "accelerationZ")
    return(results)
}
