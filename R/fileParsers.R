parseTLElines <- function(lines) {
    if(length(lines) == 2) {
        line1 <- lines[1]
        line2 <- lines[2]
    } else if (length(lines) == 3) {
        line1 <- lines[2]
        line2 <- lines[3]
    } else {
        stop("Please provide a character vector with two or three elements")
    }
    launchYearShort <- as.numeric(substr(line1, 10, 11))
    if (is.na(launchYearShort)) {
        launchYear <- "Unknown"
    } else if(launchYearShort >= 80) {
        launchYear <- 1900 + launchYearShort
    } else {
        launchYear <- 2000 + launchYearShort
    }
    launchNumber = substr(line1, 12, 14)
    launchPiece = gsub(" ", "", substr(line1, 15, 17), fixed=TRUE)
    epochYearShort <- as.numeric(substr(line1, 19, 20))
    if (is.numeric(epochYearShort)) if(epochYearShort >= 80) {
        epochYear <- 1900 + epochYearShort
    } else {
        epochYear <- 2000 + epochYearShort
    }
    yearFraction <- as.numeric(substr(line1, 21, 32))
    epochDays <- floor(yearFraction)
    daysFraction <- yearFraction%%1
    epochHours <- daysFraction*24
    epochWholeHours <- floor(epochHours)
    epochMinutes <- epochHours%%1*60
    epochWholeMinutes <- floor(epochMinutes)
    epochSeconds <- epochMinutes%%1*60
    date <- as.Date(paste(epochYear, "-01-01", sep=""), tz="UTC") + epochDays - 1
    epochDateTimeString <- paste(as.character(date), " ", epochWholeHours, ":",
                                 epochWholeMinutes, ":", epochSeconds, sep="")
    parsedTLE <- list(
        NORADcatalogNumber = substr(line1, 3, 7),
        classificationLevel = switch(substr(line1, 8, 8),
                                     "U" = {"unclassified"},
                                     "C" = {"classified"},
                                     "S" = {"secret"},
                                     "unknown"),
        internationalDesignator = paste(launchYear, "-", launchNumber, launchPiece, sep=""),
        launchYear = launchYear,
        launchNumber = launchNumber,
        launchPiece = launchPiece,
        dateTime = epochDateTimeString,
        elementNumber = as.numeric(substr(line1, 65, 68)),
        inclination = as.numeric(substr(line2, 9, 16)),
        ascension = as.numeric(substr(line2, 18, 25)),
        eccentricity = as.numeric(paste("0.", substr(line2, 27, 33), sep="")),
        perigeeArgument = as.numeric(substr(line2, 35, 42)),
        meanAnomaly = as.numeric(substr(line2, 44, 51)),
        meanMotion = as.numeric(substr(line2, 53, 63)),
        meanMotionDerivative = 2*as.numeric(substr(line1, 34, 43)),
        meanMotionSecondDerivative = 6*as.numeric(paste(substr(line1, 45, 45), "0.",
                                                        substr(line1, 46, 50), "e-",
                                                        substr(line1, 52, 52), sep="")),
        Bstar = as.numeric(paste(substr(line1, 54, 54), "0.",
                                 substr(line1, 55, 59), "e-",
                                 substr(line1, 61, 61), sep="")),
        ephemerisType = switch(substr(line1, 63, 63),
                               "0" = {"Distributed data (SGP4/SDP4)"},
                               "1" = {"SGP"},
                               "2" = {"SGP4"},
                               "3" = {"SDP4"},
                               "4" = {"SGP8"},
                               "5" = {"SDP8"},
                               "unknown"),
        epochRevolutionNumber = as.numeric(substr(line2, 64, 68))
    )
    if(length(lines) == 3) {
        parsedTLE <- c(parsedTLE, objectName = trimws(lines[1]))
    }
    return(parsedTLE)
}

readTLE <- function(filename, maxTLEs=NULL) {
    if(grepl("http", filename)) {
        tmpFile <- tempfile("asteRiskTLE")
        download.file(filename, tmpFile)
        filename <- tmpFile
    }
    line1 <- readLines(filename, n=1)
    if(substr(line1, 1, 1) != "1" | nchar(trimws(line1)) <= 12) {
        numberLines <- 3
    } else {
        numberLines <- 2
    }
    if(is.null(maxTLEs)) {
        linesToRead <- -1
    } else {
        linesToRead <- maxTLEs * numberLines
    }
    lines <- readLines(filename, n=linesToRead)
    lines <- lines[lines != ""]
    # Determine if TLE are actually in 2 or 3 line format
    # if(nchar(trimws(lines[1])) <= 12) {
    numberTLEs <- length(lines)/numberLines
    if(numberTLEs == 1) {
        parsedTLEs <- parseTLElines(lines)
    } else {
        startingLines <- seq(1, by=numberLines, length.out=numberTLEs)
        parsedTLEs <- vector(mode = "list", length = numberTLEs)
        for(i in 1:length(startingLines)) {
            singleTLElines <- lines[startingLines[i]:(startingLines[i]+numberLines-1)]
            parsedTLEs[[i]] <- parseTLElines(singleTLElines)
        }
    }
    return(parsedTLEs)
}

parseGLONASSNavigationRINEXlines <- function(lines, tauC=0, rinexVersion, flagV305) {
    if(!(length(lines) == 4 || length(lines) == 5)) {
        stop("Invalid GLONASS navigation RINEX file")
    }
    line1 <- gsub("D", "E", lines[1], ignore.case=TRUE)
    line2 <- gsub("D", "E", lines[2], ignore.case=TRUE)
    line3 <- gsub("D", "E", lines[3], ignore.case=TRUE)
    line4 <- gsub("D", "E", lines[4], ignore.case=TRUE)
    if(rinexVersion == 2) {
        satelliteSlotNumber <- as.numeric(substr(line1, 1, 2))
        satelliteNumber <- paste("R", formatC(satelliteSlotNumber, width=2, flag="0"), sep="")
        epochYearShort <- as.numeric(substr(line1, 4, 5))
        epochMonth <- as.numeric(substr(line1, 7, 8))
        epochDay <- as.numeric(substr(line1, 10, 11))
        epochHour <- as.numeric(substr(line1, 13, 14))
        epochMinute <- as.numeric(substr(line1, 16, 17))
        epochSecond <- as.numeric(substr(line1, 18, 22))
        clockBias <- -as.numeric(substr(line1, 23, 41))
        relativeFreqBias <- as.numeric(substr(line1, 42, 60))
        messageFrameTime <- as.numeric(substr(line1, 61, 79))
        positionX <- as.numeric(substr(line2, 4, 22))
        velocityX <- as.numeric(substr(line2, 23, 41))
        accelX <- as.numeric(substr(line2, 42, 60))
        satelliteHealthCode <- as.numeric(substr(line2, 61, 79))
        if(satelliteHealthCode == 0) {
            satelliteHealthStatus <- "Healthy"
        } else if(satelliteHealthCode == 1) {
            satelliteHealthStatus <- "Malfunction"
        }
        positionY <- as.numeric(substr(line3, 4, 22))
        velocityY <- as.numeric(substr(line3, 23, 41))
        accelY <- as.numeric(substr(line3, 42, 60))
        freqNumber <- as.numeric(substr(line3, 61, 79))
        positionZ <- as.numeric(substr(line4, 4, 22))
        velocityZ <- as.numeric(substr(line4, 23, 41))
        accelZ <- as.numeric(substr(line4, 42, 60))
        informationAge <- as.numeric(substr(line4, 61, 79))
        if (is.numeric(epochYearShort)) if(epochYearShort >= 80) {
            epochYear <- 1900 + epochYearShort
        } else {
            epochYear <- 2000 + epochYearShort
        }
    } else if(rinexVersion == 3) {
        satelliteNumber <- substr(line1, 1, 3)
        epochYear <- as.numeric(substr(line1, 5, 8))
        epochMonth <- as.numeric(substr(line1, 10, 11))
        epochDay <- as.numeric(substr(line1, 13, 14))
        epochHour <- as.numeric(substr(line1, 16, 17))
        epochMinute <- as.numeric(substr(line1, 19, 20))
        epochSecond <- as.numeric(substr(line1, 21, 23))
        clockBias <- -as.numeric(substr(line1, 24, 42))
        relativeFreqBias <- as.numeric(substr(line1, 43, 61))
        messageFrameTime <- as.numeric(substr(line1, 62, 80))
        positionX <- as.numeric(substr(line2, 5, 23))
        velocityX <- as.numeric(substr(line2, 24, 42))
        accelX <- as.numeric(substr(line2, 43, 61))
        satelliteHealthCode <- as.numeric(substr(line2, 62, 80))
        if(satelliteHealthCode == 0) {
            satelliteHealthStatus <- "Healthy"
        } else if(satelliteHealthCode == 1) {
            satelliteHealthStatus <- "Malfunction"
        }
        positionY <- as.numeric(substr(line3, 5, 23))
        velocityY <- as.numeric(substr(line3, 24, 42))
        accelY <- as.numeric(substr(line3, 43, 61))
        freqNumber <- as.numeric(substr(line3, 62, 80))
        positionZ <- as.numeric(substr(line4, 5, 23))
        velocityZ <- as.numeric(substr(line4, 24, 42))
        accelZ <- as.numeric(substr(line4, 43, 61))
        informationAge <- as.numeric(substr(line4, 62, 80))
        if(flagV305) {
            line5 <- gsub("D", "E", lines[5], ignore.case=TRUE)
            statusFlagsInt <- as.numeric(substr(line5, 5, 23))
            if(is.na(statusFlagsInt) || statusFlagsInt >= 512 || statusFlagsInt < 0) {
                GLONASSType <- NA
                dataUpdatedFlag <- NA
                numberSatellitesAlmanac <- NA
                ephemerisValidityTimeInterval <- NA
                parityEphemerisValidityTimeInterval <- NA
                tauCSource <- NA
                tauGPSSource <- NA
            } else {
                statusFlagsBits <- intToBits(statusFlagsInt)[1:9]
                GLONASSType <- if(statusFlagsBits[9]) "GLO-M/K" else "GLO"
                updatedDataFlag <- if(statusFlagsBits[7]) TRUE else FALSE
                numberSatellitesAlmanac <- if(statusFlagsBits[6]) 5 else 4
                if(statusFlagsBits[3]) {
                    if(statusFlagsBits[4]) {
                        ephemerisValidityTimeInterval <- 60
                    } else {
                        ephemerisValidityTimeInterval <- 45
                    }
                } else {
                    if(statusFlagsBits[4]) {
                        ephemerisValidityTimeInterval <- 30
                    } else {
                        ephemerisValidityTimeInterval <- 0
                    }
                }
                parityEphemerisValidityTimeInterval <- if(statusFlagsBits[5]) "Odd" else "Even"
                tauCSource <- if(statusFlagsBits[1]) "On-board" else "Ground"
                tauGPSSource <- if(statusFlagsBits[2]) "On-board" else "Ground"
            }
            totalGroupDelay <- substr(line5, 24, 42)
            if(trimws(totalGroupDelay) == ".999999999999E+09") {
                totalGroupDelay <- NA
            } else {
                totalGroupDelay <- as.numeric(totalGroupDelay)
            }
            URAI <- as.numeric(substr(line5, 43, 61))
            healthFlagsInt <- as.numeric(substr(line5, 62, 80))
            if(is.na(healthFlagsInt) || healthFlagsInt >= 8 || healthFlagsInt < 0) {
                satelliteHealthStatus2 <- NA
                almanacHealthStatus <- NA
            } else {
                healthFlagsBits <- intToBits(healthFlagsInt)[1:3]
                satelliteHealthStatus2 <- if(healthFlagsBits[3]) "Malfunction" else "Healthy"
                if(satelliteHealthStatus2 != satelliteHealthStatus) {
                    warning("Conflicting satellite health status codes")
                }
                if(healthFlagsBits[2]) {
                    almanacHealthStatus <- if(healthFlagsBits[1]) "Healthy" else "Malfunction"
                } else {
                    almanacHealthStatus <- NA
                }
            }
        }
    }
    dateTimeString <- paste(epochYear, "-", epochMonth, "-", epochDay, " ", 
                            epochHour, ":", epochMinute, ":", epochSecond, 
                            sep="")
    correctedEphemerisUTC <- as.nanotime(as.POSIXct(dateTimeString, tz="UTC")) + clockBias*1e9 + tauC*1e9 # AQUIAQUIAQUI
    resultsList <- list(
        satelliteNumber=satelliteNumber,
        epochYear=epochYear,
        epochMonth=epochMonth,
        epochDay=epochDay,
        epochHour=epochHour,
        epochMinute=epochMinute,
        epochSecond=epochSecond,
        ephemerisUTCTime=correctedEphemerisUTC,
        clockBias=clockBias,
        relativeFreqBias=relativeFreqBias,
        messageFrameTime=messageFrameTime,
        positionX=positionX,
        positionY=positionY,
        positionZ=positionZ,
        velocityX=velocityX,
        velocityY=velocityY,
        velocityZ=velocityZ,
        accelX=accelX,
        accelY=accelY,
        accelZ=accelZ,
        satelliteHealthStatus=satelliteHealthStatus,
        freqNumber=freqNumber,
        informationAge=informationAge
    )
    if(flagV305) {
        resultsList <- c(resultsList, list(
            GLONASSType=GLONASSType,
            updatedDataFlag=updatedDataFlag,
            numberSatellitesAlmanac=numberSatellitesAlmanac,
            ephemerisValidityTimeInterval=ephemerisValidityTimeInterval,
            parityEphemerisValidityTimeInterval=parityEphemerisValidityTimeInterval,
            tauCSource=tauCSource,
            tauGPSSource=tauGPSSource,
            totalGroupDelay,
            URAI=URAI,
            almanacHealthStatus=almanacHealthStatus
        ))
    }
    return(resultsList)
}

parseGLONASSNavigationRINEXheaderLines <- function(lines) {
    line1 <- lines[1]
    rinexVersion <- trimws(substr(line1, 1, 9))
    rinexMajorVersion <- strsplit(rinexVersion, "\\.")[[1]][1]
    if(rinexMajorVersion == "2") {
        parsedHeader <- parseGLONASSNavigationRINEXheaderV2Lines(lines)
    } else if(rinexMajorVersion == "3") {
        parsedHeader <- parseNavigationRINEXheaderV3Lines(lines)
    }
    return(parsedHeader)
}

parseGLONASSNavigationRINEXheaderV2Lines <- function(lines) {
    line1 <- lines[1]
    line2 <- lines[2]
    rinexVersion <- trimws(substr(line1, 1, 9))
    rinexFileType <- trimws(substr(line1, 21, 60))
    generatorProgram <- trimws(substr(line2, 1, 20))
    generatorEntity <- trimws(substr(line2, 21, 40))
    fileCreationDateString <- trimws(substr(line2, 41, 60))
    commentLinesIndexes <- grep("COMMENT", lines)
    if(length(commentLinesIndexes) > 0) {
        comments <- trimws(substr(lines[commentLinesIndexes], 1, 60))
    } else {
        comments <- NULL
    }
    systemTimeCorrectionLineIndex <- grep("CORR TO SYSTEM TIME", lines)
    if(length(systemTimeCorrectionLineIndex) > 0) {
        systemTimeCorrectionLine <- lines[systemTimeCorrectionLineIndex]
        refYear <- as.numeric(substr(systemTimeCorrectionLine, 1, 6))
        refMonth <- as.numeric(substr(systemTimeCorrectionLine, 7, 12))
        refDay <- as.numeric(substr(systemTimeCorrectionLine, 13, 18))
        # sysTimeCorrection <- -as.numeric(gsub("D", "E", substr(systemTimeCorrectionLine, 22, 40)))
        sysTimeCorrection <- as.numeric(gsub("D", "E", substr(systemTimeCorrectionLine, 22, 40)))
        #TODO CHECK IF THIS SHOULD HAVE - OR NOT
    } else {
        systemTimeCorrectionLine <- NULL
        refYear <- NULL
        refMonth <- NULL
        refDay <- NULL
        sysTimeCorrection <- NULL
    }
    leapSecondsLineIndex <- grep("LEAP SECONDS", lines)
    if(length(leapSecondsLineIndex) > 0) {
        leapSeconds <- as.numeric(substr(lines[leapSecondsLineIndex], 1, 6))
        # Note: number of leap seconds since 6th of January, 1980. 
        # Recommended for mixed GLONASS/GPS
    } else {
        leapSeconds <- NA
    }
    return(list(
        rinexVersion=rinexVersion,
        rinexFileType=rinexFileType,
        generatorProgram=generatorProgram,
        generatorEntity=generatorEntity,
        fileCreationDateString=fileCreationDateString,
        refYear=refYear,
        refMonth=refMonth,
        refDay=refDay,
        sysTimeCorrection=sysTimeCorrection,
        leapSeconds=leapSeconds,
        comments=comments
    ))
}

parseNavigationRINEXheaderV3Lines <- function(lines) {
    line1 <- lines[1]
    line2 <- lines[2]
    rinexVersion <- trimws(substr(line1, 1, 9))
    rinexFileType <- substr(line1, 21, 21)
    if(rinexFileType != "N") {
        stop("Wrong RINEX file type field")
    }
    satelliteSystem <- substr(line1, 41, 41)
    generatorProgram <- trimws(substr(line2, 1, 20))
    generatorEntity <- trimws(substr(line2, 21, 40))
    fileCreationDateString <- trimws(substr(line2, 41, 60))
    commentLinesIndexes <- grep("COMMENT", lines)
    if(length(commentLinesIndexes) > 0) {
        comments <- trimws(substr(lines[commentLinesIndexes], 1, 60))
    } else {
        comments <- NULL
    }
    ionosphericCorrectionIndexes <- grep("IONOSPHERIC CORR", lines)
    if(length(ionosphericCorrectionIndexes) > 0) {
        ionosphericCorrections <- vector(mode="list", length = length(ionosphericCorrectionIndexes))
        for(index in ionosphericCorrectionIndexes) {
            ionosphericCorrectionLine <- lines[index]
            ionosphericCorrectionType <- trimws(substr(ionosphericCorrectionLine, 1, 4))
            coefficient1 <- as.numeric(gsub("D", "E", substr(ionosphericCorrectionLine, 6, 17)))
            coefficient2 <- as.numeric(gsub("D", "E", substr(ionosphericCorrectionLine, 18, 29)))
            coefficient3 <- as.numeric(gsub("D", "E", substr(ionosphericCorrectionLine, 30, 41)))
            coefficient4 <- as.numeric(gsub("D", "E", substr(ionosphericCorrectionLine, 42, 53)))
            timeMark <- trimws(substr(ionosphericCorrectionLine, 54, 55))
            SVID <- trimws(substr(ionosphericCorrectionLine, 56, 58))
            namesCoefficients <- switch(ionosphericCorrectionType,
                                        GAL=c("ai0", "ai1", "ai2"),
                                        GPSA=c("alpha0", "alpha1", "alpha2", "alpha3"),
                                        GPSB=c("beta0", "beta1", "beta2", "beta3"),
                                        QZSA=c("alpha0", "alpha1", "alpha2", "alpha3"),
                                        QZSB=c("beta0", "beta1", "beta2", "beta3"),
                                        BDSA=c("alpha0", "alpha1", "alpha2", "alpha3"),
                                        BDSB=c("beta0", "beta1", "beta2", "beta3"),
                                        IRNA=c("alpha0", "alpha1", "alpha2", "alpha3"),
                                        IRNB=c("beta0", "beta1", "beta2", "beta3"))
            if(ionosphericCorrectionType == "GAL") {
                ionosphericCorrectionList <- list(
                    ionosphericCorrectionType, coefficient1, coefficient2, coefficient3, 
                    timeMark, SVID
                )
            } else {
                ionosphericCorrectionList <- list(
                    ionosphericCorrectionType, coefficient1, coefficient2, coefficient3, 
                    coefficient4, timeMark, SVID
                )
            }
            names(ionosphericCorrectionList) <- c("ionosphericCorrectionType", namesCoefficients,
                                                  "ionosphericCorrectionTimeMark", "ionosphericCorrectionSV")
            ionosphericCorrections[[index]] <- ionosphericCorrectionList
        }
        comments <- trimws(substr(lines[commentLinesIndexes], 1, 60))
    } else {
        ionosphericCorrections <- NULL
    }
    systemTimeCorrectionLineIndex <- grep("TIME SYSTEM CORR", lines)
    if(length(systemTimeCorrectionLineIndex) > 0) {
        systemTimeCorrectionLine <- lines[systemTimeCorrectionLineIndex]
        systemTimeCorrectionType <- trimws(substr(systemTimeCorrectionLine, 1, 4))
        timeCorrectionA0 <- as.numeric(gsub("D", "E", substr(systemTimeCorrectionLine, 6, 22)))
        timeCorrectionA1 <- as.numeric(gsub("D", "E", substr(systemTimeCorrectionLine, 23, 38)))
        timeCorrectionReferenceTime <- as.numeric(substr(systemTimeCorrectionLine, 39, 45))
        timeCorrectionReferenceWeekNumber <- as.numeric(substr(systemTimeCorrectionLine, 46, 50))
        timeCorrectionSatelliteNumber <- trimws(substr(systemTimeCorrectionLine, 51, 57))
        UTCTypeIdentifier <- as.numeric(substr(systemTimeCorrectionLine, 58, 60))
        UTCType <- switch(UTCTypeIdentifier+1,
                          "Uknown",
                          "UTC(NIST)",
                          "UTC(USNO)",
                          "UTC(SU)",
                          "UTC(BIPM)",
                          "UTC(Europe Lab)",
                          "UTC(CRL)",
                          "UTC(NTSC) (BDS)",
                          "Unassigned")
    } else {
        systemTimeCorrectionType <- NA
        timeCorrectionA0 <- NA
        timeCorrectionA1 <- NA
        timeCorrectionReferenceTime <- NA
        timeCorrectionReferenceWeekNumber <- NA
        timeCorrectionSatelliteNumber <- NA
        UTCTypeIdentifier <- NA
        UTCType <- NA
    }
    leapSecondsLineIndex <- grep("LEAP SECONDS", lines)
    if(length(leapSecondsLineIndex) > 0) {
        leapSecondsLine <- lines[leapSecondsLineIndex]
        leapSeconds <- as.numeric(substr(leapSecondsLine, 1, 6))
        deltaTimeLeapSeconds <- as.numeric(substr(leapSecondsLine, 7, 12))
        deltaTimeLeapSecondsWeekNumber <- as.numeric(substr(leapSecondsLine, 13, 18))
        deltaTimeLeapSecondsDayNumber <- as.numeric(substr(leapSecondsLine, 19, 24))
        leapSecondsTimeSystemIdentifier <- substr(leapSecondsLine, 25, 27)
        if(leapSecondsTimeSystemIdentifier != "BDS") leapSecondsTimeSystemIdentifier <- "GPS"
    } else {
        leapSeconds <- NA
        deltaTimeLeapSeconds <- NA
        deltaTimeLeapSecondsWeekNumber <- NA
        deltaTimeLeapSecondsDayNumber <- NA
        leapSecondsTimeSystemIdentifier <- NA
    }
    return(list(
        rinexVersion=rinexVersion,
        rinexFileType=rinexFileType,
        satelliteSystem=satelliteSystem,
        generatorProgram=generatorProgram,
        generatorEntity=generatorEntity,
        fileCreationDateString=fileCreationDateString,
        systemTimeCorrectionType=systemTimeCorrectionType,
        timeCorrectionA0=timeCorrectionA0,
        timeCorrectionA1=timeCorrectionA1,
        timeCorrectionReferenceTime=timeCorrectionReferenceTime,
        timeCorrectionReferenceWeekNumber=timeCorrectionReferenceWeekNumber,
        timeCorrectionSatelliteNumber=timeCorrectionSatelliteNumber,
        UTCType=UTCType,
        leapSeconds=leapSeconds,
        deltaTimeLeapSeconds=deltaTimeLeapSeconds,
        deltaTimeLeapSecondsWeekNumber=deltaTimeLeapSecondsWeekNumber,
        deltaTimeLeapSecondsDayNumber=deltaTimeLeapSecondsDayNumber,
        leapSecondsTimeSystemIdentifier=leapSecondsTimeSystemIdentifier,
        ionosphericCorrections=ionosphericCorrections,
        comments=comments
    ))
}

readGLONASSNavigationRINEX <- function(filename) {
    lines <- readLines(filename)
    lines <- lines[lines != ""]
    endOfHeader <- grep("END OF HEADER", lines)
    if(length(endOfHeader) == 0) {
        stop("Header missing END OF HEADER label")
    }
    headerLines <- lines[1:endOfHeader]
    bodyLines <- lines[(endOfHeader+1):length(lines)]
    headerFields <- parseGLONASSNavigationRINEXheaderLines(headerLines)
    rinexVersion <- headerFields$rinexVersion
    rinexMajorVersion <- substr(rinexVersion, 1, 1)
    rinexVersionCompare <- compareVersion(rinexVersion, "3.05")
    if(rinexVersionCompare >= 0) {
        flagV305 <- TRUE
    } else {
        flagV305 <- FALSE
    }
    messageNumberLines <- if(flagV305) 5 else 4
    numberMessages <- length(bodyLines)/messageNumberLines
    parsedMessages <- vector(mode = "list", length = numberMessages)
    startingLines <- seq(1, by=messageNumberLines, length.out = numberMessages)
    for(i in 1:length(startingLines)) {
        singleMessageLines <- bodyLines[startingLines[i]:(startingLines[i]+messageNumberLines-1)]
        parsedMessages[[i]] <- parseGLONASSNavigationRINEXlines(singleMessageLines, 
                                                                tauC=headerFields$sysTimeCorrection,
                                                                rinexVersion = rinexMajorVersion,
                                                                flagV305 = flagV305)
    }
    return(list(
        header=headerFields, 
        messages=parsedMessages))
}

parseGPSNavigationRINEXlines <- function(lines, leapSeconds=0, deltaUTCA0=0,
                                         deltaUTCA1=0, referenceTimeUTC,
                                         referenceWeekUTC, rinexVersion) {
    if(length(lines) != 8) {
        stop("Invalid GPS navigation RINEX file")
    }
    if(is.na(leapSeconds)) leapSeconds <- 0
    line1 <- gsub("D", "E", lines[1], ignore.case=TRUE)
    line2 <- gsub("D", "E", lines[2], ignore.case=TRUE)
    line3 <- gsub("D", "E", lines[3], ignore.case=TRUE)
    line4 <- gsub("D", "E", lines[4], ignore.case=TRUE)
    line5 <- gsub("D", "E", lines[5], ignore.case=TRUE)
    line6 <- gsub("D", "E", lines[6], ignore.case=TRUE)
    line7 <- gsub("D", "E", lines[7], ignore.case=TRUE)
    line8 <- gsub("D", "E", lines[8], ignore.case=TRUE)
    if(rinexVersion == 2) {
        satellitePRNCode <- as.numeric(substr(line1, 1, 2))
        satelliteNumber <- paste("G", formatC(satellitePRNCode, width=2, flag="0"), sep="")
        tocYearShort <- as.numeric(substr(line1, 4, 5))
        tocMonth <- as.numeric(substr(line1, 7, 8))
        tocDay <- as.numeric(substr(line1, 10, 11))
        tocHour <- as.numeric(substr(line1, 13, 14))
        tocMinute <- as.numeric(substr(line1, 16, 17))
        tocSecond <- as.numeric(substr(line1, 18, 22))
        clockBias <- as.numeric(substr(line1, 23, 41))
        clockDrift <- as.numeric(substr(line1, 42, 60))
        clockDriftRate <- as.numeric(substr(line1, 61, 79))
        IODE <- as.numeric(substr(line2, 4, 22))
        radiusCorrectionSine <- as.numeric(substr(line2, 23, 41))
        deltaN <- as.numeric(substr(line2, 42, 60))
        meanAnomaly <- as.numeric(substr(line2, 61, 79)) # radians, M0 angle
        latitudeCorrectionCosine <- as.numeric(substr(line3, 4, 22)) 
        eccentricity <- as.numeric(substr(line3, 23, 41))
        latitudeCorrectionSine <- as.numeric(substr(line3, 42, 60))
        semiMajorAxis <- as.numeric(substr(line3, 61, 79))^2 # RINEX file contains sqrt
        calculatedMeanMotion <- semiMajorAxisToMeanMotion(semiMajorAxis, outputRevsPerDay = FALSE)
        correctedMeanMotion <- calculatedMeanMotion + deltaN
        toeSecondsOfGPSWeek <- as.numeric(substr(line4, 4, 22))
        inclinationCorrectionCosine <- as.numeric(substr(line4, 23, 41))
        ascension <- as.numeric(substr(line4, 42, 60)) # radians, OMEGA angle
        inclinationCorrectionSine <- as.numeric(substr(line4, 61, 79))
        inclination <- as.numeric(substr(line5, 4, 22)) # radians, initial inclination
        radiusCorrectionCosine <- as.numeric(substr(line5, 23, 41))
        perigeeArgument <- as.numeric(substr(line5, 42, 60)) # radians, omega angle
        OMEGADot <- as.numeric(substr(line5, 61, 79)) # radians/s, derivative of OMEGA
        inclinationRate <- as.numeric(substr(line6, 4, 22)) 
        codesL2Channel <- as.numeric(substr(line6, 23, 41))
        toeGPSWeek <- as.numeric(substr(line6, 42, 60))
        dataFlagL2P <- as.numeric(substr(line6, 61, 79))
        satelliteAccuracy <- as.numeric(substr(line7, 4, 22)) 
        satelliteHealthCode <- as.numeric(substr(line7, 23, 41))
        totalGroupDelay <- as.numeric(substr(line7, 42, 60)) 
        IODC <- as.numeric(substr(line7, 61, 79))
        transmissionTime <- as.numeric(substr(line8, 4, 22)) # seconds of GPS week
        fitInterval <- as.numeric(substr(line8, 23, 41))
        if (is.numeric(tocYearShort)) if(tocYearShort >= 80) {
            tocYear <- 1900 + tocYearShort
        } else {
            tocYear <- 2000 + tocYearShort
        }
    } else if(rinexVersion == 3) {
        satelliteNumber <- substr(line1, 1, 3)
        tocYear <- as.numeric(substr(line1, 5, 8))
        tocMonth <- as.numeric(substr(line1, 10, 11))
        tocDay <- as.numeric(substr(line1, 13, 14))
        tocHour <- as.numeric(substr(line1, 16, 17))
        tocMinute <- as.numeric(substr(line1, 19, 20))
        tocSecond <- as.numeric(substr(line1, 22, 23))
        clockBias <- as.numeric(substr(line1, 24, 42))
        clockDrift <- as.numeric(substr(line1, 43, 61))
        clockDriftRate <- as.numeric(substr(line1, 62, 80))
        IODE <- as.numeric(substr(line2, 5, 23))
        radiusCorrectionSine <- as.numeric(substr(line2, 24, 42))
        deltaN <- as.numeric(substr(line2, 43, 61))
        meanAnomaly <- as.numeric(substr(line2, 62, 80)) # radians, M0 angle
        latitudeCorrectionCosine <- as.numeric(substr(line3, 5, 23)) 
        eccentricity <- as.numeric(substr(line3, 24, 42))
        latitudeCorrectionSine <- as.numeric(substr(line3, 43, 61))
        semiMajorAxis <- as.numeric(substr(line3, 62, 80))^2 # RINEX file contains sqrt
        calculatedMeanMotion <- semiMajorAxisToMeanMotion(semiMajorAxis, outputRevsPerDay = FALSE)
        correctedMeanMotion <- calculatedMeanMotion + deltaN
        toeSecondsOfGPSWeek <- as.numeric(substr(line4, 5, 23))
        inclinationCorrectionCosine <- as.numeric(substr(line4, 24, 42))
        ascension <- as.numeric(substr(line4, 43, 61)) # radians, OMEGA angle
        inclinationCorrectionSine <- as.numeric(substr(line4, 62, 80))
        inclination <- as.numeric(substr(line5, 5, 23)) # radians, initial inclination
        radiusCorrectionCosine <- as.numeric(substr(line5, 24, 42))
        perigeeArgument <- as.numeric(substr(line5, 43, 61)) # radians, omega angle
        OMEGADot <- as.numeric(substr(line5, 62, 80)) # radians/s, derivative of OMEGA
        inclinationRate <- as.numeric(substr(line6, 5, 23)) 
        codesL2Channel <- as.numeric(substr(line6, 24, 42))
        toeGPSWeek <- as.numeric(substr(line6, 43, 61))
        dataFlagL2P <- as.numeric(substr(line6, 62, 80))
        satelliteAccuracy <- as.numeric(substr(line7, 5, 23)) 
        satelliteHealthCode <- as.numeric(substr(line7, 24, 42))
        totalGroupDelay <- as.numeric(substr(line7, 43, 61)) 
        IODC <- as.numeric(substr(line7, 62, 80))
        transmissionTime <- as.numeric(substr(line8, 5, 23)) # seconds of GPS week
        fitInterval <- as.numeric(substr(line8, 24, 42))
    }
    tocDateTimeString <- paste(tocYear, "-", tocMonth, "-", tocDay, " ", 
                               tocHour, ":", tocMinute, ":", tocSecond, 
                               sep="")
    tocGPSTseconds <- as.double(difftime(tocDateTimeString, "1980-01-06 00:00:00", tz="UTC", units="secs")) 
    tocGPSWeek <- tocGPSTseconds%/%604800
    tocSecondsOfGPSWeek <- tocGPSTseconds%%604800
    differenceToeToc <- toeSecondsOfGPSWeek - tocSecondsOfGPSWeek
    if(differenceToeToc > 302400) differenceToeToc <- differenceToeToc - 604800
    if(differenceToeToc < -302400) differenceToeToc <- differenceToeToc + 604800
    F_constant <- -2*sqrt(GM_Earth_TDB)/c_light^2
    kepler_sol_1 <- meanAnomaly
    kepler_sol_current <- kepler_sol_1
    convergence <- FALSE
    iterations <- 0
    keplerAccuracy=10e-8
    maxKeplerIterations=100
    while(!convergence) {
        iterations <- iterations + 1
        delta_kepler_sol <- ( (meanAnomaly - kepler_sol_current + eccentricity*sin(kepler_sol_current)) / (1 - eccentricity*cos(kepler_sol_current)) )
        kepler_sol_current <- kepler_sol_current + delta_kepler_sol
        if(iterations > maxKeplerIterations | delta_kepler_sol < keplerAccuracy) {
            convergence <- TRUE
        }
    }
    eccentricAnomaly <- kepler_sol_current
    relativityDeltaSVtimeToGPST <- F_constant * eccentricity*(sqrt(semiMajorAxis)) * sin(eccentricAnomaly)
    deltaSVtimeToGPST <- clockBias + clockDrift*differenceToeToc + clockDriftRate*differenceToeToc^2 + relativityDeltaSVtimeToGPST
    toeSecondsOfGPSWeekGPST <- toeSecondsOfGPSWeek - deltaSVtimeToGPST
    differenceToeGPSTTot <- toeSecondsOfGPSWeekGPST - referenceTimeUTC
    deltaGPSTToUTC <- leapSeconds + deltaUTCA0 + deltaUTCA1*(differenceToeGPSTTot + 604800*(tocGPSWeek - referenceWeekUTC))
    # ephemerisUTCDateTime <- as.POSIXct(toeGPSWeek * 604800 + toeSecondsOfGPSWeek - deltaSVtimeToGPST - deltaGPSTToUTC, origin="1980-01-06 00:00:00", tz="UTC")
    # ephemerisUTCDateTimeString <- as.character(ephemerisUTCDateTime) # In UTC
    # Modified to use nanotime package to handle nanosecond accuracy
    ephemerisUTCDateTime <- as.POSIXct(toeGPSWeek * 604800 + toeSecondsOfGPSWeek, origin="1980-01-06 00:00:00", tz="UTC")
    ephemerisUTCDateTime <- as.nanotime(ephemerisUTCDateTime) - deltaSVtimeToGPST*1e9 - deltaGPSTToUTC*1e9
    # Calculation of ECEF ephemeris according to IS-GPS-200M
    # trueAnomaly <- 2 * atan2(sqrt(1 + eccentricity) * sin(eccentricAnomaly/2), sqrt(1 - eccentricity) * cos(eccentricAnomaly/2))
    trueAnomaly <- atan2(sqrt(1 - eccentricity^2) * sin(eccentricAnomaly), cos(eccentricAnomaly) - eccentricity)
    # trueAnomaly <- 2 * atan(sqrt((1 + eccentricity) / (1 - eccentricity)) * tan(eccentricAnomaly/2))
    latitudeArgument <- trueAnomaly + perigeeArgument
    deltauk <- latitudeCorrectionSine * sin(2*latitudeArgument) + latitudeCorrectionCosine * cos(2*latitudeArgument)
    deltark <- radiusCorrectionSine * sin(2*latitudeArgument) + radiusCorrectionCosine * cos(2*latitudeArgument)
    deltaik <- inclinationCorrectionSine * sin(2*latitudeArgument) + inclinationCorrectionCosine * cos(2*latitudeArgument)
    correctedLatitudeArgument <- latitudeArgument + deltauk
    correctedRadius <- semiMajorAxis * (1 - eccentricity*cos(eccentricAnomaly)) + deltark
    correctedInclination <- inclination + deltaik
    # Positions in orbital plane
    xkprime <- correctedRadius*cos(correctedLatitudeArgument)
    ykprime <- correctedRadius*sin(correctedLatitudeArgument)
    correctedAscension <- ascension - omegaEarth*toeSecondsOfGPSWeek
    position_ECEF <- c( # in m
        xkprime*cos(correctedAscension) - ykprime*cos(correctedInclination)*sin(correctedAscension),
        xkprime*sin(correctedAscension) + ykprime*cos(correctedInclination)*cos(correctedAscension),
        ykprime*sin(correctedInclination))
    # velocity
    eccentricAnomalyDot <- correctedMeanMotion/(1 - eccentricity * cos(eccentricAnomaly))
    trueAnomalyDot <- eccentricAnomalyDot*sqrt(1 - eccentricity^2)/(1 - eccentricity*cos(eccentricAnomaly))
    correctedInclinationDot <- inclinationRate + 2*trueAnomalyDot*(
        inclinationCorrectionSine*cos(2*latitudeArgument) - inclinationCorrectionCosine*sin(2*latitudeArgument))
    correctedLatitudeArgumentDot <- trueAnomalyDot + 2*trueAnomalyDot*(
        latitudeCorrectionSine*cos(2*latitudeArgument) - latitudeCorrectionCosine*sin(2*latitudeArgument))
    correctedRadiusDot <- eccentricity * semiMajorAxis * eccentricAnomalyDot *
        sin(eccentricAnomaly) + 2*trueAnomalyDot*(
            radiusCorrectionSine*cos(2*latitudeArgument) - radiusCorrectionCosine*sin(2*latitudeArgument)
        )
    ascensionDot <- OMEGADot - omegaEarth
    # In-plane velocities
    xkprimedot <- correctedRadiusDot*cos(correctedLatitudeArgument) - 
        correctedRadius*correctedLatitudeArgumentDot*sin(correctedLatitudeArgument)
    ykprimedot <- correctedRadiusDot*sin(correctedLatitudeArgument) + 
        correctedRadius*correctedLatitudeArgumentDot*cos(correctedLatitudeArgument)
    velocity_ECEF <- c( # in m/s
        -xkprime*ascensionDot*sin(correctedAscension) + xkprimedot*cos(correctedAscension) - 
            ykprimedot*sin(correctedAscension)*cos(correctedInclination) -
            ykprime*(ascensionDot*cos(correctedAscension)*cos(correctedInclination) - correctedInclinationDot*sin(correctedAscension)*sin(correctedInclination)),
        xkprime*ascensionDot*cos(correctedAscension) + xkprimedot*sin(correctedAscension) + 
            ykprimedot*cos(correctedAscension)*cos(correctedInclination) -
            ykprime*(ascensionDot*sin(correctedAscension)*cos(correctedInclination) + correctedInclinationDot*sin(correctedAscension)*sin(correctedInclination)),
        ykprimedot*sin(correctedInclination) + ykprime*correctedInclinationDot*cos(correctedInclination)
    )
    # Acceleration
    oblateEarthAccelerationFactor <- -1.5*J2_WGS84*(GM_Earth_TDB/correctedRadius^2)*(earthRadius_WGS84/correctedRadius)^2
    acceleration_ECEF <- c(
        -GM_Earth_TDB*(position_ECEF[1]/correctedRadius^3) + oblateEarthAccelerationFactor*
            ((1 - 5*(position_ECEF[3]/correctedRadius)^2)*(position_ECEF[1]/correctedRadius)) +
            2*velocity_ECEF[2]*omegaEarth + position_ECEF[1]*omegaEarth^2,
        -GM_Earth_TDB*(position_ECEF[2]/correctedRadius^3) + oblateEarthAccelerationFactor*
            ((1 - 5*(position_ECEF[3]/correctedRadius)^2)*(position_ECEF[1]/correctedRadius)) -
            2*velocity_ECEF[1]*omegaEarth + position_ECEF[2]*omegaEarth^2,
        -GM_Earth_TDB*(position_ECEF[3]/correctedRadius^3) + oblateEarthAccelerationFactor*
            ((3-5*(position_ECEF[3]/correctedRadius)^2)*(position_ECEF[3]/correctedRadius))
    )
    return(list(
        satelliteNumber=satelliteNumber,
        tocYear=tocYear,
        tocMonth=tocMonth,
        tocDay=tocDay,
        tocHour=tocHour,
        tocMinute=tocMinute,
        tocSecond=tocSecond,
        clockBias=clockBias,
        clockDrift=clockDrift,
        clockDriftRate=clockDriftRate,
        IODE=IODE,
        radiusCorrectionSine=radiusCorrectionSine,
        deltaN=deltaN,
        correctedMeanMotion=correctedMeanMotion,
        meanAnomaly=meanAnomaly,
        latitudeCorrectionCosine=latitudeCorrectionCosine,
        eccentricity=eccentricity,
        latitudeCorrectionSine=latitudeCorrectionSine,
        semiMajorAxis=semiMajorAxis,
        toeSecondsOfGPSWeek=toeSecondsOfGPSWeek,
        inclinationCorrectionCosine=inclinationCorrectionCosine,
        ascension=ascension,
        inclinationCorrectionSine=inclinationCorrectionSine,
        inclination=inclination,
        radiusCorrectionCosine=radiusCorrectionCosine,
        perigeeArgument=perigeeArgument,
        OMEGADot=OMEGADot,
        inclinationRate=inclinationRate,
        codesL2Channel=codesL2Channel,
        toeGPSWeek=toeGPSWeek,
        dataFlagL2P=dataFlagL2P,
        satelliteAccuracy=satelliteAccuracy,
        satelliteHealthCode=satelliteHealthCode,
        totalGroupDelay=totalGroupDelay,
        IODC=IODC,
        transmissionTime=transmissionTime,
        fitInterval=fitInterval,
        ephemerisUTCTime=ephemerisUTCDateTime,
        position_ITRF=position_ECEF,
        velocity_ITRF=velocity_ECEF,
        acceleration_ITRF=acceleration_ECEF
    ))
}

parseGPSNavigationRINEXheaderV2Lines <- function(lines) {
    line1 <- lines[1]
    line2 <- lines[2]
    rinexVersion <- trimws(substr(line1, 1, 9))
    rinexFileType <- trimws(substr(line1, 21, 60))
    generatorProgram <- trimws(substr(line2, 1, 20))
    generatorEntity <- trimws(substr(line2, 21, 40))
    fileCreationDateString <- trimws(substr(line2, 41, 60))
    commentLinesIndexes <- grep("COMMENT", lines)
    if(length(commentLinesIndexes) > 0) {
        comments <- trimws(substr(lines[commentLinesIndexes], 1, 60))
    } else {
        comments <- NULL
    }
    ionAlphaLineIndex <- grep("ION ALPHA", lines)
    if(length(ionAlphaLineIndex) > 0) {
        ionAlphaLine <- lines[ionAlphaLineIndex]
        ionAlphaFields <- strsplit(trimws(substr(ionAlphaLine, 1, 60)), split="\\s+")[[1]]
        ionAlphaA0 <- as.numeric(gsub("D", "E", ionAlphaFields[1]))
        ionAlphaA1 <- as.numeric(gsub("D", "E", ionAlphaFields[2]))
        ionAlphaA2 <- as.numeric(gsub("D", "E", ionAlphaFields[3]))
        ionAlphaA3 <- as.numeric(gsub("D", "E", ionAlphaFields[4]))
    } else {
        ionAlphaA0 <- ionAlphaA1 <- ionAlphaA2 <- ionAlphaA3 <- NULL
    }
    ionBetaLineIndex <- grep("ION BETA", lines)
    if(length(ionBetaLineIndex) > 0) {
        ionBetaLine <- lines[ionBetaLineIndex]
        ionBetaFields <- strsplit(trimws(substr(ionBetaLine, 1, 60)), split="\\s+")[[1]]
        ionBetaB0 <- as.numeric(gsub("D", "E", ionBetaFields[1]))
        ionBetaB1 <- as.numeric(gsub("D", "E", ionBetaFields[2]))
        ionBetaB2 <- as.numeric(gsub("D", "E", ionBetaFields[3]))
        ionBetaB3 <- as.numeric(gsub("D", "E", ionBetaFields[4]))
    } else {
        ionBetaB0 <- ionBetaB1 <- ionBetaB2 <- ionBetaB3 <- NULL
    }
    deltaUTCLineIndex <- grep("DELTA.UTC", lines)
    if(length(deltaUTCLineIndex) > 0) {
        deltaUTCLine <- lines[deltaUTCLineIndex]
        deltaUTCA0 <- as.numeric(gsub("D", "E", substr(deltaUTCLine, 1, 22)))
        deltaUTCA1 <- as.numeric(gsub("D", "E", substr(deltaUTCLine, 23, 42)))
        referenceTimeUTC <- as.numeric(substr(deltaUTCLine, 43, 51))
        referenceWeekUTC <- as.numeric(substr(deltaUTCLine, 52, 60))
    } else {
        deltaUTCLine <- deltaUTCA0 <- deltaUTCA1 <- referenceTimeUTC <- referenceWeekUTC <- NULL
    }
    leapSecondsLineIndex <- grep("LEAP SECONDS", lines)
    if(length(leapSecondsLineIndex) > 0) {
        leapSeconds <- as.numeric(substr(lines[leapSecondsLineIndex], 1, 6))
        # Note: number of leap seconds since 6th of January, 1980. 
        # Recommended for mixed GLONASS/GPS
    } else {
        leapSeconds <- NA
    }
    return(list(
        rinexVersion=rinexVersion,
        rinexFileType=rinexFileType,
        generatorProgram=generatorProgram,
        generatorEntity=generatorEntity,
        fileCreationDateString=fileCreationDateString,
        ionAlphaA0=ionAlphaA0,
        ionAlphaA1=ionAlphaA1,
        ionAlphaA2=ionAlphaA2,
        ionAlphaA3=ionAlphaA3,
        ionBetaB0=ionBetaB0,
        ionBetaB1=ionBetaB1,
        ionBetaB2=ionBetaB2,
        ionBetaB3=ionBetaB3,
        deltaUTCA0=deltaUTCA0,
        deltaUTCA1=deltaUTCA1,
        referenceTimeUTC=referenceTimeUTC,
        referenceWeekUTC=referenceWeekUTC,
        leapSeconds=leapSeconds,
        comments=comments
    ))
}

parseGPSNavigationRINEXheaderLines <- function(lines) {
    line1 <- lines[1]
    rinexVersion <- trimws(substr(line1, 1, 9))
    rinexMajorVersion <- strsplit(rinexVersion, "\\.")[[1]][1]
    if(rinexMajorVersion == "2") {
        parsedHeader <- parseGPSNavigationRINEXheaderV2Lines(lines)
    } else if(rinexMajorVersion == "3") {
        parsedHeader <- parseNavigationRINEXheaderV3Lines(lines)
    }
    return(parsedHeader)
}

readGPSNavigationRINEX <- function(filename) {
    lines <- readLines(filename)
    lines <- lines[lines != ""]
    endOfHeader <- grep("END OF HEADER", lines)
    if(length(endOfHeader) == 0) {
        stop("Header missing END OF HEADER label")
    }
    headerLines <- lines[1:endOfHeader]
    bodyLines <- lines[(endOfHeader+1):length(lines)]
    headerFields <- parseGPSNavigationRINEXheaderLines(headerLines)
    rinexVersion <- headerFields$rinexVersion
    rinexMajorVersion <- substr(rinexVersion, 1, 1)
    messageNumberLines <- 8
    numberMessages <- length(bodyLines)/messageNumberLines
    parsedMessages <- vector(mode = "list", length = numberMessages)
    startingLines <- seq(1, by=messageNumberLines, length.out = numberMessages)
    for(i in 1:length(startingLines)) {
        singleMessageLines <- bodyLines[startingLines[i]:(startingLines[i]+messageNumberLines-1)]
        parsedMessages[[i]] <- parseGPSNavigationRINEXlines(singleMessageLines, headerFields$leapSeconds, 
                                                            headerFields$deltaUTCA0,
                                                            headerFields$deltaUTCA1, 
                                                            headerFields$referenceTimeUTC,
                                                            headerFields$referenceWeekUTC,
                                                            rinexMajorVersion)
    }
    return(list(
        header=headerFields, 
        messages=parsedMessages))
}

parseOMheader <- function(lines) {
    OMVersion <- trimws(strsplit(lines[1], split="=")[[1]][2])
    creationDate <- gsub("T", " ", trimws(strsplit(lines[length(lines) - 1], split="=")[[1]][2]))
    creator <- trimws(strsplit(lines[length(lines)], split="=")[[1]][2])
    return(list(
        OMVersion=OMVersion,
        creationDate=creationDate,
        creator=creator
    ))
}

parseOEMDataBlock <- function(lines) {
    metadataStartLineIndex <- grep("META_START", lines) + 1
    metadataEndLineIndex <- grep("META_STOP", lines) - 1
    covarianceMatrixStartLineIndex <- grep("COVARIANCE_START", lines) + 1
    if(length(covarianceMatrixStartLineIndex) > 0) {
        hasCovarianceMatrix <- TRUE
        covarianceMatrixEndLineIndex <- grep("COVARIANCE_STOP", lines) - 1
        ephemeridesEndLineIndex <- covarianceMatrixStartLineIndex - 2
    } else {
        hasCovarianceMatrix <- FALSE
        covarianceMatrixEndLineIndex <- numeric(0)
        ephemeridesEndLineIndex <- length(lines)
    }
    commentLinesIndexes <- grep("COMMENT", lines)
    metadataCommentLinesIndexes <- commentLinesIndexes[(commentLinesIndexes >= metadataStartLineIndex) &
                                                           (commentLinesIndexes <= metadataEndLineIndex)]
    if(hasCovarianceMatrix) {
        ephemeridesCommentLinesIndexes <- commentLinesIndexes[(commentLinesIndexes > metadataEndLineIndex) &
                                                                  (commentLinesIndexes < covarianceMatrixStartLineIndex)]
        covarianceMatrixCommentLinesIndexes <- commentLinesIndexes[commentLinesIndexes >= covarianceMatrixStartLineIndex]
    } else {
        ephemeridesCommentLinesIndexes <- commentLinesIndexes[commentLinesIndexes > metadataEndLineIndex]
        covarianceMatrixCommentLinesIndexes <- NULL
    }
    ephemeridesStartLineIndex <- max(metadataEndLineIndex + 2, ephemeridesCommentLinesIndexes[length(ephemeridesCommentLinesIndexes)] + 1)
    objectNameLineIndex <- grep("OBJECT_NAME", lines)
    objectName <- trimws(strsplit(lines[objectNameLineIndex], split="=")[[1]][2])
    objectID <- trimws(strsplit(lines[objectNameLineIndex+1], split="=")[[1]][2])
    refFrameCentralBody <- trimws(strsplit(lines[objectNameLineIndex+2], split="=")[[1]][2])
    refFrame <- trimws(strsplit(lines[objectNameLineIndex+3], split="=")[[1]][2])
    refFrameEpochLineIndex <-  grep("REF_FRAME_EPOCH", lines)
    if(length(refFrameEpochLineIndex) > 0) {
        refFrameEpoch <- gsub("T", " ", trimws(strsplit(lines[refFrameEpochLineIndex], split="=")[[1]][2]))
    } else {
        refFrameEpoch <- NULL
    }
    timeSystemLineIndex <- grep("TIME_SYSTEM", lines)
    timeSystem <- trimws(strsplit(lines[timeSystemLineIndex], split="=")[[1]][2])
    startTimeLineIndex <- grep("START_TIME", lines)
    startTime <- gsub("T", " ", trimws(strsplit(lines[startTimeLineIndex], split="=")[[1]][2]))
    stopTimeLineIndex <- grep("STOP_TIME", lines)
    stopTime <- gsub("T", " ", trimws(strsplit(lines[stopTimeLineIndex], split="=")[[1]][2]))
    usableStartTimeLineIndex <- grep("USEABLE_START_TIME", lines)
    usableStartTime <- gsub("T", " ", trimws(strsplit(lines[usableStartTimeLineIndex], split="=")[[1]][2]))
    usableStopTimeLineIndex <- grep("USEABLE_STOP_TIME", lines)
    usableStopTime <- gsub("T", " ", trimws(strsplit(lines[usableStopTimeLineIndex], split="=")[[1]][2]))
    interpolationMethodLineIndex <- grep("INTERPOLATION", lines)
    if(length(interpolationMethodLineIndex) > 0) {
        interpolationMethod <- trimws(strsplit(lines[interpolationMethodLineIndex], split="=")[[1]][2])
        interpolationOrder <- trimws(strsplit(lines[interpolationMethodLineIndex+1], split="=")[[1]][2])
    } else {
        interpolationMethod <- interpolationOrder <- NULL
    }
    massLineIndex <- grep("MASS", lines)
    if(length(massLineIndex) > 0) {
        mass <- as.numeric(trimws(strsplit(lines[massLineIndex], split="=")[[1]][2]))
    } else {
        mass <- NULL
    }
    dragAreaLineIndex <- grep("DRAG_AREA", lines)
    if(length(dragAreaLineIndex) > 0) {
        dragArea <- as.numeric(trimws(strsplit(lines[dragAreaLineIndex], split="=")[[1]][2]))
    } else {
        dragArea <- NULL
    }
    dragCoeffLineIndex <- grep("DRAG_COEFF", lines)
    if(length(dragCoeffLineIndex) > 0) {
        dragCoeff <- as.numeric(trimws(strsplit(lines[dragCoeffLineIndex], split="=")[[1]][2]))
    } else {
        dragCoeff <- NULL
    }
    solarRadAreaLineIndex <- grep("SOLAR_RAD_AREA", lines)
    if(length(solarRadAreaLineIndex) > 0) {
        solarRadArea <- as.numeric(trimws(strsplit(lines[solarRadAreaLineIndex], split="=")[[1]][2]))
    } else {
        solarRadArea <- NULL
    }
    solarRadCoeffLineIndex <- grep("SOLAR_RAD_COEFF", lines)
    if(length(solarRadCoeffLineIndex) > 0) {
        solarRadCoeff <- as.numeric(trimws(strsplit(lines[solarRadCoeffLineIndex], split="=")[[1]][2]))
    } else {
        solarRadCoeff <- NULL
    }
    ephemeridesColNumbers <- length(strsplit(trimws(lines[ephemeridesStartLineIndex]), split="\\s+")[[1]])
    if(ephemeridesColNumbers == 7) {
        ephemeridesColNames <- c("epoch", "position_X", "position_Y", "position_Z", 
                                 "velocity_X", "velocity_Y", "velocity_Z")
    } else if (ephemeridesColNumbers == 10) {
        ephemeridesColNames <- c("epoch", "position_X", "position_Y", "position_Z", 
                                 "velocity_X", "velocity_Y", "velocity_Z",
                                 "acceleration_X", "acceleration_Y", "acceleration_Z")
    } else {
        stop("Invalid number of columns in ephemerides data block.")
    }
    ephemerides <- read.table(header=FALSE, text=paste(lines[ephemeridesStartLineIndex:ephemeridesEndLineIndex],
                                                       collapse="\n"),
                              col.names = ephemeridesColNames)
    ephemerides[, 1] <- gsub("T", " ", ephemerides[, 1])
    if(hasCovarianceMatrix) {
        covarianceMatrixLines <- lines[covarianceMatrixStartLineIndex:covarianceMatrixEndLineIndex]
        singleCovarianceMatrixStartLineIndexes <- grep("EPOCH", covarianceMatrixLines)
        singleCovarianceMatrixEndLineIndexes <- c(singleCovarianceMatrixStartLineIndexes[-1] - 1,
                                                  covarianceMatrixEndLineIndex)
        numberCovarianceMatrixes <- length(singleCovarianceMatrixStartLineIndexes)
        parsedCovarianceMatrixes <- vector(mode = "list", length = numberCovarianceMatrixes)
        for(i in 1:length(singleCovarianceMatrixStartLineIndexes)) {
            singleCovarianceMatrixLines <- covarianceMatrixLines[singleCovarianceMatrixStartLineIndexes[i]:singleCovarianceMatrixEndLineIndexes[i]]
            covarianceMatrixEpoch <- gsub("T", " ", trimws(strsplit(singleCovarianceMatrixLines[1], split="=")[[1]][2]))
            covarianceMatrixRefFrameLineIndex <- grep("COV_REF_FRAME", singleCovarianceMatrixLines)
            if(length(covarianceMatrixRefFrameLineIndex) > 0) {
                covarianceMatrixRefFrame <- trimws(strsplit(singleCovarianceMatrixLines[covarianceMatrixRefFrameLineIndex], split="=")[[1]][2])
                covarianceMatrixDataLines <- singleCovarianceMatrixLines[3:8]
            } else {
                covarianceMatrixRefFrame <- refFrame
                covarianceMatrixDataLines <- singleCovarianceMatrixLines[2:7]
            }
            covarianceMatrix <- data.matrix(read.table(text=paste(covarianceMatrixDataLines, collapse="\n"), 
                                                       fill=TRUE, header=FALSE, col.names=1:6))
            colnames(covarianceMatrix) <- NULL
            covarianceMatrix[upper.tri(covarianceMatrix)] <- t(covarianceMatrix)[upper.tri(covarianceMatrix)]
            parsedCovarianceMatrixes[[i]] <- list(
                epoch=covarianceMatrixEpoch,
                referenceFrame=covarianceMatrixRefFrame,
                covarianceMatrix=covarianceMatrix
            )
        }
        # if(length(parsedCovarianceMatrixes) == 1) {
        #     parsedCovarianceMatrixes <- parsedCovarianceMatrixes[[1]]
        # }
    } else {
        parsedCovarianceMatrixes=NULL
    }
    return(list(
        objectName=objectName,
        objectID=objectID,
        referenceFrame=refFrame,
        refFrameEpoch=refFrameEpoch,
        centerName=refFrameCentralBody,
        timeSystem=timeSystem,
        startTime=startTime,
        endTime=stopTime,
        usableStartTime=usableStartTime,
        usableEndTime=usableStopTime,
        interpolationMethod=interpolationMethod,
        interpolationOrder=interpolationOrder,
        mass=mass,
        dragArea=dragArea,
        dragCoefficient=dragCoeff,
        solarRadArea=solarRadArea,
        solarRadCoefficient=solarRadCoeff,
        ephemerides=ephemerides,
        covarianceMatrixes=parsedCovarianceMatrixes
    ))
}

readOEM <- function(filename) {
    lines <- readLines(filename)
    lines <- lines[lines != ""]
    endOfHeader <- grep("ORIGINATOR", lines)
    headerFields <- parseOMheader(lines[1:endOfHeader])
    startLines <- grep("META_START", lines)
    endLines <- c(startLines[-1] - 1, length(lines))
    numberDataBlocks <- length(startLines)
    parsedDataBlocks <- vector(mode = "list", length = numberDataBlocks)
    for(i in 1:length(startLines)) {
        singleDataBlockLines <- lines[startLines[i]:endLines[i]]
        parsedDataBlocks[[i]] <- parseOEMDataBlock(singleDataBlockLines)
    }
    # if(length(parsedDataBlocks) == 1) {
    #     parsedDataBlocks <- parsedDataBlocks[[1]]
    # }
    return(list(
        header=headerFields, 
        dataBlocks=parsedDataBlocks))
}

parseOMMDataBlock <- function(lines) {
    metadataStartLineIndex <- grep("META_START", lines) + 1
    metadataEndLineIndex <- grep("META_STOP", lines) - 1
    covarianceMatrixStartLineIndex <- grep("COVARIANCE_START", lines) + 1
    if(length(covarianceMatrixStartLineIndex) > 0) {
        hasCovarianceMatrix <- TRUE
        covarianceMatrixEndLineIndex <- grep("COVARIANCE_STOP", lines) - 1
        ephemeridesEndLineIndex <- covarianceMatrixStartLineIndex - 2
    } else {
        hasCovarianceMatrix <- FALSE
        covarianceMatrixEndLineIndex <- numeric(0)
        ephemeridesEndLineIndex <- length(lines)
    }
    commentLinesIndexes <- grep("COMMENT", lines)
    metadataCommentLinesIndexes <- commentLinesIndexes[(commentLinesIndexes >= metadataStartLineIndex) &
                                                           (commentLinesIndexes <= metadataEndLineIndex)]
    if(hasCovarianceMatrix) {
        ephemeridesCommentLinesIndexes <- commentLinesIndexes[(commentLinesIndexes > metadataEndLineIndex) &
                                                                  (commentLinesIndexes < covarianceMatrixStartLineIndex)]
        covarianceMatrixCommentLinesIndexes <- commentLinesIndexes[commentLinesIndexes >= covarianceMatrixStartLineIndex]
    } else {
        ephemeridesCommentLinesIndexes <- commentLinesIndexes[commentLinesIndexes > metadataEndLineIndex]
        covarianceMatrixCommentLinesIndexes <- NULL
    }
    ephemeridesStartLineIndex <- max(metadataEndLineIndex + 2, ephemeridesCommentLinesIndexes[length(ephemeridesCommentLinesIndexes)] + 1)
    objectNameLineIndex <- grep("OBJECT_NAME", lines)
    objectName <- trimws(strsplit(lines[objectNameLineIndex], split="=")[[1]][2])
    objectID <- trimws(strsplit(lines[objectNameLineIndex+1], split="=")[[1]][2])
    refFrameCentralBody <- trimws(strsplit(lines[objectNameLineIndex+2], split="=")[[1]][2])
    refFrame <- trimws(strsplit(lines[objectNameLineIndex+3], split="=")[[1]][2])
    refFrameEpochLineIndex <-  grep("REF_FRAME_EPOCH", lines)
    if(length(refFrameEpochLineIndex) > 0) {
        refFrameEpoch <- gsub("T", " ", trimws(strsplit(lines[refFrameEpochLineIndex], split="=")[[1]][2]))
    } else {
        refFrameEpoch <- NULL
    }
    timeSystemLineIndex <- grep("TIME_SYSTEM", lines)
    timeSystem <- trimws(strsplit(lines[timeSystemLineIndex], split="=")[[1]][2])
    startTimeLineIndex <- grep("START_TIME", lines)
    startTime <- gsub("T", " ", trimws(strsplit(lines[startTimeLineIndex], split="=")[[1]][2]))
    stopTimeLineIndex <- grep("STOP_TIME", lines)
    stopTime <- gsub("T", " ", trimws(strsplit(lines[stopTimeLineIndex], split="=")[[1]][2]))
    usableStartTimeLineIndex <- grep("USEABLE_START_TIME", lines)
    usableStartTime <- gsub("T", " ", trimws(strsplit(lines[usableStartTimeLineIndex], split="=")[[1]][2]))
    usableStopTimeLineIndex <- grep("USEABLE_STOP_TIME", lines)
    usableStopTime <- gsub("T", " ", trimws(strsplit(lines[usableStopTimeLineIndex], split="=")[[1]][2]))
    interpolationMethodLineIndex <- grep("INTERPOLATION", lines)
    if(length(interpolationMethodLineIndex) > 0) {
        interpolationMethod <- trimws(strsplit(lines[interpolationMethodLineIndex], split="=")[[1]][2])
        interpolationOrder <- trimws(strsplit(lines[interpolationMethodLineIndex+1], split="=")[[1]][2])
    } else {
        interpolationMethod <- interpolationOrder <- NULL
    }
    massLineIndex <- grep("MASS", lines)
    if(length(massLineIndex) > 0) {
        mass <- as.numeric(trimws(strsplit(lines[massLineIndex], split="=")[[1]][2]))
    } else {
        mass <- NULL
    }
    dragAreaLineIndex <- grep("DRAG_AREA", lines)
    if(length(dragAreaLineIndex) > 0) {
        dragArea <- as.numeric(trimws(strsplit(lines[dragAreaLineIndex], split="=")[[1]][2]))
    } else {
        dragArea <- NULL
    }
    dragCoeffLineIndex <- grep("DRAG_COEFF", lines)
    if(length(dragCoeffLineIndex) > 0) {
        dragCoeff <- as.numeric(trimws(strsplit(lines[dragCoeffLineIndex], split="=")[[1]][2]))
    } else {
        dragCoeff <- NULL
    }
    solarRadAreaLineIndex <- grep("SOLAR_RAD_AREA", lines)
    if(length(solarRadAreaLineIndex) > 0) {
        solarRadArea <- as.numeric(trimws(strsplit(lines[solarRadAreaLineIndex], split="=")[[1]][2]))
    } else {
        solarRadArea <- NULL
    }
    solarRadCoeffLineIndex <- grep("SOLAR_RAD_COEFF", lines)
    if(length(solarRadCoeffLineIndex) > 0) {
        solarRadCoeff <- as.numeric(trimws(strsplit(lines[solarRadCoeffLineIndex], split="=")[[1]][2]))
    } else {
        solarRadCoeff <- NULL
    }
    ephemeridesColNumbers <- length(strsplit(trimws(lines[ephemeridesStartLineIndex]), split="\\s+")[[1]])
    if(ephemeridesColNumbers == 7) {
        ephemeridesColNames <- c("epoch", "position_X", "position_Y", "position_Z", 
                                 "velocity_X", "velocity_Y", "velocity_Z")
    } else if (ephemeridesColNumbers == 10) {
        ephemeridesColNames <- c("epoch", "position_X", "position_Y", "position_Z", 
                                 "velocity_X", "velocity_Y", "velocity_Z",
                                 "acceleration_X", "acceleration_Y", "acceleration_Z")
    } else {
        stop("Invalid number of columns in ephemerides data block.")
    }
    ephemerides <- read.table(header=FALSE, text=paste(lines[ephemeridesStartLineIndex:ephemeridesEndLineIndex],
                                                       collapse="\n"),
                              col.names = ephemeridesColNames)
    ephemerides[, 1] <- gsub("T", " ", ephemerides[, 1])
    if(hasCovarianceMatrix) {
        covarianceMatrixLines <- lines[covarianceMatrixStartLineIndex:covarianceMatrixEndLineIndex]
        singleCovarianceMatrixStartLineIndexes <- grep("EPOCH", covarianceMatrixLines)
        singleCovarianceMatrixEndLineIndexes <- c(singleCovarianceMatrixStartLineIndexes[-1] - 1,
                                                  covarianceMatrixEndLineIndex)
        numberCovarianceMatrixes <- length(singleCovarianceMatrixStartLineIndexes)
        parsedCovarianceMatrixes <- vector(mode = "list", length = numberCovarianceMatrixes)
        for(i in 1:length(singleCovarianceMatrixStartLineIndexes)) {
            singleCovarianceMatrixLines <- covarianceMatrixLines[singleCovarianceMatrixStartLineIndexes[i]:singleCovarianceMatrixEndLineIndexes[i]]
            covarianceMatrixEpoch <- gsub("T", " ", trimws(strsplit(singleCovarianceMatrixLines[1], split="=")[[1]][2]))
            covarianceMatrixRefFrameLineIndex <- grep("COV_REF_FRAME", singleCovarianceMatrixLines)
            if(length(covarianceMatrixRefFrameLineIndex) > 0) {
                covarianceMatrixRefFrame <- trimws(strsplit(singleCovarianceMatrixLines[covarianceMatrixRefFrameLineIndex], split="=")[[1]][2])
                covarianceMatrixDataLines <- singleCovarianceMatrixLines[3:8]
            } else {
                covarianceMatrixRefFrame <- refFrame
                covarianceMatrixDataLines <- singleCovarianceMatrixLines[2:7]
            }
            covarianceMatrix <- data.matrix(read.table(text=paste(covarianceMatrixDataLines, collapse="\n"), 
                                                       fill=TRUE, header=FALSE, col.names=1:6))
            colnames(covarianceMatrix) <- NULL
            covarianceMatrix[upper.tri(covarianceMatrix)] <- t(covarianceMatrix)[upper.tri(covarianceMatrix)]
            parsedCovarianceMatrixes[[i]] <- list(
                epoch=covarianceMatrixEpoch,
                referenceFrame=covarianceMatrixRefFrame,
                covarianceMatrix=covarianceMatrix
            )
        }
        # if(length(parsedCovarianceMatrixes) == 1) {
        #     parsedCovarianceMatrixes <- parsedCovarianceMatrixes[[1]]
        # }
    } else {
        parsedCovarianceMatrixes=NULL
    }
    return(list(
        objectName=objectName,
        objectID=objectID,
        referenceFrame=refFrame,
        refFrameEpoch=refFrameEpoch,
        centerName=refFrameCentralBody,
        timeSystem=timeSystem,
        startTime=startTime,
        endTime=stopTime,
        usableStartTime=usableStartTime,
        usableEndTime=usableStopTime,
        interpolationMethod=interpolationMethod,
        interpolationOrder=interpolationOrder,
        mass=mass,
        dragArea=dragArea,
        dragCoefficient=dragCoeff,
        solarRadArea=solarRadArea,
        solarRadCoefficient=solarRadCoeff,
        ephemerides=ephemerides,
        covarianceMatrixes=parsedCovarianceMatrixes
    ))
}

readPlanetLabsStateVectors <- function(filename) {
    
}

.readBinDAF <- function(filename) {
    recordLengthBytes <- 1024
    fileSize <- file.info(filename)$size
    rawData <- readBin(filename, "raw", n=fileSize*2)
    fileRecord <- parseBinDAFFileRecord(rawData[1:recordLengthBytes])
    if((fileRecord$endianString == "LTL-IEEE" && .Platform$endian == "big") || 
       (fileRecord$endianString == "BIG-IEEE" && .Platform$endian == "little")) {
        rawData <- readBin(filename, "raw", n=fileSize*2, endian = "swap")
        fileRecord <- parseBinDAFFileRecord(rawData[1:recordLengthBytes])
    }
    firstSummaryRecord <- fileRecord$firstSummaryRecNum
    if(firstSummaryRecord > 2) {
        lastCommentsAddress <- recordLengthBytes*(firstSummaryRecord-1)
        comments <- parseBinDAFCommentsRecords(rawData[(recordLengthBytes+1):lastCommentsAddress])
    } else {
        comments <- NULL
    }
    numIntsSummary <- fileRecord$numIntsSummary
    numDoublesSummary <- fileRecord$numDoublesSummary
    summaryRecordSizeDoubles <- fileRecord$summaryRecordSizeDoubles
    numCharsName <- fileRecord$numCharsName
    currentSummaryRecord <- parseBinDAFSummaryRecord(rawData[(lastCommentsAddress+1):(lastCommentsAddress+1024)],
                                                     numDoublesSummary, numIntsSummary, summaryRecordSizeDoubles)
    currentNameRecord <- parseBinDAFNameRecord(rawData[(lastCommentsAddress+1025):(lastCommentsAddress+2048)],
                                               currentSummaryRecord$numberOfSummaries, numCharsName)
    initialArrayAddressesDoubles <- sapply(currentSummaryRecord$summaries, `[[`, "initialArrayAddress")
    finalArrayAddressesDoubles <- sapply(currentSummaryRecord$summaries, `[[`, "finalArrayAddress")
    arraysBytesIndexes <- Map(`:`, 8*(initialArrayAddressesDoubles - 1) + 1, 8*finalArrayAddressesDoubles) 
    currentArrays <- lapply(arraysBytesIndexes, function(i) readBin(rawData[i], what="double", n=length(i)/8))
    allSummaries <- currentSummaryRecord$summaries
    allNames <- currentNameRecord
    allArrays <- currentArrays
    while(currentSummaryRecord$nextSummaryRecord != 0) {
        newSummaryRecordInitialAddress <- recordLengthBytes*(currentSummaryRecord$nextSummaryRecord-1) + 1
        currentSummaryRecord <- parseBinDAFSummaryRecord(rawData[newSummaryRecordInitialAddress:(newSummaryRecordInitialAddress+1023)],
                                                         numDoublesSummary, numIntsSummary, summaryRecordSizeDoubles)
        currentNameRecord <- parseBinDAFNameRecord(rawData[(newSummaryRecordInitialAddress+1024):(newSummaryRecordInitialAddress+2047)],
                                                   currentSummaryRecord$numberOfSummaries, numCharsName)
        initialArrayAddressesDoubles <- sapply(currentSummaryRecord$summaries, `[[`, "initialArrayAddress")
        finalArrayAddressesDoubles <- sapply(currentSummaryRecord$summaries, `[[`, "finalArrayAddress")
        arraysBytesIndexes <- Map(`:`, 8*(initialArrayAddressesDoubles - 1) + 1, 8*finalArrayAddressesDoubles) 
        currentArrays <- lapply(arraysBytesIndexes, function(i) readBin(rawData[i], what="double", n=length(i)/8))
        allSummaries <- c(allSummaries, currentSummaryRecord$summaries)
        allNames <- c(allNames, currentNameRecord)
        allArrays <- c(allArrays, currentArrays)
    }
    return(list(
        fileRecord=fileRecord,
        comments=comments,
        summaries=allSummaries,
        names=allNames,
        arrays=allArrays
    ))
}

readBinDAF <- function(filename) {
    DAF <- .readBinDAF(filename)
    allNames <- DAF$names
    allSummaries <- DAF$summaries
    allArrays <- DAF$arrays
    arraysList <- vector(mode="list", length=length(allArrays))
    for(i in 1:length(arraysList)) {
        arraysList[[i]] <- list(
            arrayName=allNames[i],
            arraySummary=allSummaries[[i]],
            arrayElements=allArrays[[i]]
        )
    }
    return(list(
        metadata=DAF$fileRecord,
        comments=DAF$comments,
        arrays=arraysList
    ))
}

parseBinDAFFileRecord <- function(rawData) {
    fileType <- trimws(rawToChar(rawData[1:8]), "right")
    numDoublesSummary <- rawToInt(rawData[9:12])
    numIntsSummary <- rawToInt(rawData[13:16])
    # This size is in doubles, i.e., multiply by 8 to get bytes
    summaryRecordSizeDoubles <- numDoublesSummary + (numIntsSummary+1)%/%2 
    numCharsName <- 8*summaryRecordSizeDoubles
    description <- trimws(rawToChar(rawData[17:76]), "right")
    firstSummaryRecNum <- rawToInt(rawData[77:80])
    lastSummaryRecNum <- rawToInt(rawData[81:84])
    firstFreeAddress <- rawToInt(rawData[85:88])
    endianString <- rawToChar(rawData[89:96])
    ftpStrBytes <- rawData[700:728]
    ftpString <- rawToChar(ftpStrBytes[ftpStrBytes!="00"])
    if(ftpString != "FTPSTR:\r:\n:\r\n:\r:\x81:\020\xce:ENDFTP") {
        warning("Corrupt FTP string. Please verify integrity of the file.")
        #TODO add ways to check FTP string better? Only match initial FTPSTR: ?
    }
    return(list(
        fileType=fileType,
        numDoublesSummary=numDoublesSummary,
        numIntsSummary=numIntsSummary,
        summaryRecordSizeDoubles=summaryRecordSizeDoubles,
        numCharsName=numCharsName,
        description=description,
        firstSummaryRecNum=firstSummaryRecNum,
        lastSummaryRecNum=lastSummaryRecNum,
        firstFreeAddress=firstFreeAddress,
        endianString=endianString,
        ftpString=ftpString
    ))
}

parseBinDAFCommentsRecords <- function(rawData) {
    commentsBytes <- rawData[1:(which(rawData=="04")-1)]
    commentsBytes[commentsBytes=="00"] <- charToRaw("\n")
    commentsLines <- strsplit(rawToChar(commentsBytes), "\n")[[1]]
}

parseBinDAFSummaryRecord <- function(rawData, numberDoubles, numberInts, summarySizeDoubles) {
    if(numberInts < 2) {
        stop("At least 2 integers should be present in each summary")
    }
    nextSummaryRecord <- rawToDouble(rawData[1:8])
    previoustSummaryRecord <- rawToDouble(rawData[9:16])
    numberOfSummaries <- rawToDouble(rawData[17:24])
    splitBytes <- split(rawData[-(1:24)], ceiling(seq_along(rawData[-(1:24)])/(summarySizeDoubles*8)))
    summaries <- lapply(splitBytes[1:numberOfSummaries], parseBinDAFSummary, 
                        numberDoubles=numberDoubles, numberInts=numberInts)
    return(list(
        nextSummaryRecord = nextSummaryRecord,
        previoustSummaryRecord = previoustSummaryRecord,
        numberOfSummaries = numberOfSummaries,
        summaries = summaries
    ))
}

parseBinDAFSummary <- function(rawData, numberDoubles, numberInts) {
    doubles <- readBin(rawData[1:(numberDoubles*8)], what="double", n=numberDoubles, size=8)
    integers <- readBin(rawData[(numberDoubles*8+1):(numberDoubles*8 + 1 + numberInts*4)], what="integer", n=numberInts, size=4)
    summary <- c(as.list(doubles), as.list(integers))
    names(summary) <- c(
        paste(rep("Double", times=numberDoubles), seq(1, by=1, length.out=numberDoubles), sep=""),
        #paste(rep("Integer", times=min(0, numberInts-2)), seq(1, by=1, length.out=min(0, numberInts-2)), sep=""),
        paste(rep("Integer", times=numberInts-2), seq(1, by=1, length.out=numberInts-2), sep=""),
        "initialArrayAddress", "finalArrayAddress"
    )
    return(summary)
}

parseBinDAFNameRecord <- function(rawData, numberSummaries, numberChars) {
    splitBytes <- split(rawData, ceiling(seq_along(rawData)/(numberChars)))
    arrayNames <- trimws(lapply(splitBytes[1:numberSummaries], rawToChar))
    return(arrayNames)
}

readBinSPK <- function(filename) {
    DAF <- .readBinDAF(filename)
    formattedArraySummaries <- lapply(DAF$summaries, formatSPKSummary)
    formattedArrays <- mapply(formatSPKArray, DAF$arrays, formattedArraySummaries, 
                              SIMPLIFY = FALSE, USE.NAMES = FALSE)
    segments <- mapply(list, DAF$names, formattedArraySummaries, formattedArrays, 
                       SIMPLIFY = FALSE, USE.NAMES = FALSE)
    for(i in 1:length(segments)){
        names(segments[[i]]) <- c("segmentName", "segmentSummary", "segmentData")
    }
    return(list(
        comments=DAF$comments,
        segments=segments
    ))
}

formatSPKSummary <- function(SPKArraySummary) {
    SPKType <- switch(SPKArraySummary[[6]],
                      "Modified Difference Arrays", # 1 - Supported
                      "Chebyshev position (equal time steps)", # 2 - Supported
                      "Chebyshev position and velocity (equal time steps)", # 3 - Supported
                      "Hubble special form", # 4 - NOT Supported
                      "Discrete states (two body propagation)", # 5 - Supported
                      "Phobos/Deimos special form", # 6 - NOT Supported
                      "Precessing classical elements", # 7 - NOT Supported
                      "Lagrange interpolation (equal time steps)", # 8 - Supported
                      "Lagrange interpolation (unequal time steps)", # 9 - Supported
                      "Two-Line Elements", # 10 - Supported
                      "Unknown", # 11 - NOT Supported (not defined yet)
                      "Hermite interpolation (equal time steps)", # 12 - Supported
                      "Hermite interpolation (unequal time steps)", # 13 - Supported
                      "Chebyshev position and velocity (unequal time steps)", # 14 - Supported
                      "Precessing conic elements", # 15 - Supported
                      "ESA Infrared Space Observatory special form", # 16 - NOT Supported
                      "Equinoctial elements", # 17 - Supported
                      "ESOC/DDID interpolation", # 18 - Supported
                      "ESOC/DDID piecewise interpolation", # 19 - Supported
                      "Chebyshev velocity (equal time steps)", # 20 - Supported
                      "Extended Modified Difference Arrays" # 21 - Supported
    )
    formattedSPKArraySummary <- c(SPKType, SPKArraySummary)
    names(formattedSPKArraySummary) <- c("SPKType", "initialEpoch", "finalEpoch", 
                                         "targetNAIFCode", "centerNAIFCode", 
                                         "frameNAIFCode", "SPKTypeCode",
                                         "initialArrayAddress", "finalArrayAddress")
    return(formattedSPKArraySummary)
}

formatSPKArray <- function(SPKArray, SPKArraySummary) {
    SPKTypeCode <- SPKArraySummary$SPKTypeCode
    if(SPKTypeCode %in% c(4, 6, 7, 11, 16)) {
        stop("Unsupported SPK type")
    } else if(SPKTypeCode == 1) {
        numberRecords <- SPKArray[length(SPKArray)]
        formattedArray <- vector(mode="list", length=numberRecords)
        finalEpochs <- SPKArray[(71*numberRecords + 1):(72*numberRecords)]
        for(i in 1:numberRecords) {
            startingIndex <- 71*(i - 1) + 1
            formattedArray[[i]] <- list(
                referenceEpoch = SPKArray[startingIndex],
                finalEpoch =finalEpochs[i],
                stepsizeFunctionVector = SPKArray[(startingIndex+1):(startingIndex+15)],
                referencePosition = SPKArray[startingIndex+c(16, 18, 20)],
                referenceVelocity = SPKArray[startingIndex+c(17, 19, 21)],
                MDA = matrix(SPKArray[(startingIndex+22):(startingIndex+66)], ncol=3),
                maxIntegrationOrderP1 = SPKArray[startingIndex+67],
                integrationOrderArray = SPKArray[(startingIndex+68):(startingIndex+70)]
            )
        }
    } else if(SPKTypeCode == 2) {
        numberRecords <- SPKArray[length(SPKArray)]
        elementsPerRecord <- SPKArray[length(SPKArray)-1]
        intervalLength <- SPKArray[length(SPKArray)-2]
        firstInitialEpoch <- SPKArray[length(SPKArray)-3]
        numberCoefficients <- (elementsPerRecord - 2)/3
        records <- SPKArray[1:(length(SPKArray)-4)]
        recordStartIndexes <- seq(from=0, by=elementsPerRecord, length.out=numberRecords)
        midPointIndexes <- 1 + recordStartIndexes
        intervalRadiiIndexes <- 2 + recordStartIndexes
        midPoint <- records[midPointIndexes]
        intervalRadius <- records[intervalRadiiIndexes]
        allCoefficients <- records[-c(midPointIndexes, intervalRadiiIndexes)]
        chebyshevCoefficients <- matrix(allCoefficients, ncol=3*numberCoefficients, nrow=numberRecords,
                                        byrow = TRUE)
        colnames(chebyshevCoefficients) <- paste(rep(c("positionXCoeff", "positionYCoeff", "positionZCoeff"), 
                                                     each=numberCoefficients), 1:numberCoefficients, sep="")
        initialEpoch <- seq(from=firstInitialEpoch, by=intervalLength, length.out=numberRecords)
        # formattedArray <- cbind(initialEpochs, midPoints, intervalRadii, chebyshevCoefficients)
        formattedArray <- list(
            polynomialDegree=numberCoefficients - 1,
            chebyshevCoefficients=cbind(initialEpoch, midPoint, intervalRadius, chebyshevCoefficients)
        )
    } else if(SPKTypeCode == 3) {
        numberRecords <- SPKArray[length(SPKArray)]
        elementsPerRecord <- SPKArray[length(SPKArray)-1]
        intervalLength <- SPKArray[length(SPKArray)-2]
        firstInitialEpoch <- SPKArray[length(SPKArray)-3]
        numberCoefficients <- (elementsPerRecord - 2)/6
        records <- SPKArray[1:(length(SPKArray)-4)]
        recordStartIndexes <- seq(from=0, by=elementsPerRecord, length.out=numberRecords)
        midPointIndexes <- 1 + recordStartIndexes
        intervalRadiiIndexes <- 2 + recordStartIndexes
        midPoint <- records[midPointIndexes]
        intervalRadius <- records[intervalRadiiIndexes]
        allCoefficients <- records[-c(midPointIndexes, intervalRadiiIndexes)]
        chebyshevCoefficients <- matrix(allCoefficients, ncol=6*numberCoefficients, nrow=numberRecords,
                                        byrow = TRUE)
        colnames(chebyshevCoefficients) <- paste(rep(c("positionXCoeff", "positionYCoeff", "positionZCoeff",
                                                       "velocityXCoeff", "velocityYCoeff", "velocityZCoeff"), 
                                                     each=numberCoefficients), 1:numberCoefficients, sep="")
        initialEpoch <- seq(from=firstInitialEpoch, by=intervalLength, length.out=numberRecords)
        # formattedArray <- cbind(initialEpoch, midPoint, intervalRadius, chebyshevCoefficients)
        formattedArray <- list(
            polynomialDegree=numberCoefficients - 1,
            chebyshevCoefficients=cbind(initialEpoch, midPoint, intervalRadius, chebyshevCoefficients)
        )
    } else if(SPKTypeCode == 5) {
        numberRecords <- SPKArray[length(SPKArray)]
        centralBodyGM <- SPKArray[length(SPKArray)-1]
        stateVectors <- SPKArray[1:(numberRecords*6)]
        epoch <- SPKArray[(numberRecords*6 + 1):(numberRecords*7 + 1)]
        stateVectorsMatrix <- matrix(stateVectors, ncol=6, nrow=numberRecords, byrow=TRUE)
        colnames(stateVectorsMatrix) <- c("positionX", "positionY", "positionZ",
                                          "velocityX", "velocityY", "velocityZ")
        formattedArray <- list(
            centralBodyGM=centralBodyGM,
            stateVectors=cbind(epoch, stateVectorsMatrix)
        )
    } else if(SPKTypeCode == 8) {
        numberRecords <- SPKArray[length(SPKArray)]
        polynomialDegree <- SPKArray[length(SPKArray)-1]
        stepSize <- SPKArray[length(SPKArray)-2]
        epoch1 <- SPKArray[length(SPKArray)-3]
        stateVectors <- SPKArray[1:(numberRecords*6)]
        epoch <- seq(from=epoch1, by=stepSize, length.out=numberRecords)
        stateVectorsMatrix <- matrix(stateVectors, ncol=6, nrow=numberRecords, byrow=TRUE)
        colnames(stateVectorsMatrix) <- c("positionX", "positionY", "positionZ",
                                          "velocityX", "velocityY", "velocityZ")
        formattedArray <- list(
            polynomialDegree=polynomialDegree,
            stateVectors=cbind(epoch, stateVectorsMatrix)
        )
    } else if(SPKTypeCode == 9) {
        numberRecords <- SPKArray[length(SPKArray)]
        polynomialDegree <- SPKArray[length(SPKArray)-1]
        stateVectors <- SPKArray[1:(numberRecords*6)]
        epoch <- SPKArray[(numberRecords*6 + 1):(numberRecords*7 + 1)]
        stateVectorsMatrix <- matrix(stateVectors, ncol=6, nrow=numberRecords, byrow=TRUE)
        colnames(stateVectorsMatrix) <- c("positionX", "positionY", "positionZ",
                                          "velocityX", "velocityY", "velocityZ")
        formattedArray <- list(
            polynomialDegree=polynomialDegree,
            stateVectors=cbind(epoch, stateVectorsMatrix)
        )
    } else if(SPKTypeCode == 10) {
        numberMetaDataItems <- SPKArray[length(SPKArray)]
        if(numberMetaDataItems < 15) {
            stop("SPK kernel with missing metadata values")
        }
        if(numberMetaDataItems == 15) {
            numberMetaDataItems <- 16
            # This is done because according to sgmeta.c comments, some old SPK 
            # kernels did not count the number of metadata items as one of the
            # items themselves
            # From reading SPKR01, it seems files with NMETA=15 actually have 16 
            # items counting NMETA itself, so these actually have 16 metadata items
            # I guess in these cases the missing metadata item is PKTOFF (latest added
            # in sgparam.inc).
        }
        metaData <- SPKArray[(length(SPKArray) - numberMetaDataItems + 1):length(SPKArray)]
        # The following should be 8 for type 10 SPKs
        numberConstants <- metaData[2]
        # I guess this should always evaluate to 1 since generic segments start
        # with constants always, unless no constants are present. But they should
        # always be for type 10 SPKs
        firstConstantsDouble <- metaData[1] + 1 
        numberDataPackets <- metaData[12]
        firstDataPacketsDouble <- metaData[11]
        dataPacketSize <- metaData[15]
        # From my tests, the type 10 SPK written out by SPICE contain the epoch
        # of each TLE before the 14 elements that each one has as described in 
        # the required reading, so for now add 1 
        #TODO VERIFY THIS
        dataPacketSize <- dataPacketSize + 1
        # These are generically called reference values in SGMETA, but in type 10
        # SPK documentation these data blocks seem to be the epochs
        numberRefValues <- metaData[7]
        firstRefValuesDouble <- metaData[6]
        constants <- as.list(SPKArray[firstConstantsDouble:(firstConstantsDouble + numberConstants - 1)])
        names(constants) <- c("J2", "J3", "J4", "sqrtGM", "highAltBound", "lowAltBound", "earthRadius", "distUnitsPerRadius")
        stateVectors <- SPKArray[(firstDataPacketsDouble + 1):(firstDataPacketsDouble + numberDataPackets * dataPacketSize)]
        stateVectorsMatrix <- matrix(stateVectors, ncol=dataPacketSize, nrow=numberDataPackets, byrow=TRUE)
        if(ncol(stateVectorsMatrix) < 15) {
            # This is required because apparently some type 10 kernels dont contain
            # obliquity and longitude nutation angles
            stateVectorsMatrix <- cbind(stateVectorsMatrix, rep(NULL, times=numberDataPackets*(15-ncol(stateVectorsMatrix))))
        }
        # remove column 11 to avoid having epoch twice
        stateVectorsMatrix <- stateVectorsMatrix[, -11, drop=FALSE]
        # What comes in the SPK type 10 seems to be:
        # epoch
        # NDT20: 1 half of 1st derivative of mean motion in radians/minute^2
        # NDD60: 1 sixth of 2nd derivative of mean motion in radians/minute^3
        # Bstar
        # Inclination in radians
        # right ascension of the ascending node in radians
        # eccentricity
        # argument of perigee in radians
        # mean anomaly in radians
        # mean motion in radians/min
        # epoch (again)
        # NU.OBLIQUITY (radians I guess) nutation angle delta psi 
        # NU.LONGITUDE (radians I guess) nutation angle delta epsilon
        #  dOBLIQUITY/dt (radians/min I guess) rate of change of delta psi
        #  dLONGITUDE/dt (radians/min I guess) rate of change of delta epsilon
        # the latter 4 quantities are calculated from Wahr series. 
        # I guess they are useful for converting TEME to GCRF
        colnames(stateVectorsMatrix) <- c("epoch", "meanMotionDerivative", "meanMotionSecondDerivative",
                                          "Bstar", "inclination", "ascension", "eccentricity",
                                          "perigeeArgument", "meanAnomaly", "meanMotion", "deltaPsi",
                                          "deltaEpsilon", "deltaPsiDerivative", "deltaEpsilonDerivative")
        stateVectorsMatrix[, 2] <- stateVectorsMatrix[, 2] * 2
        stateVectorsMatrix[, 3] <- stateVectorsMatrix[, 3] * 6
        formattedArray <- list(
            constants=constants,
            TLEs=stateVectorsMatrix
        )
    } else if(SPKTypeCode == 12) {
        numberRecords <- SPKArray[length(SPKArray)]
        windowSize <- SPKArray[length(SPKArray) - 1] + 1
        polynomialDegree <- 2*windowSize - 1
        stepSize <- SPKArray[length(SPKArray) - 2]
        epoch1 <- SPKArray[length(SPKArray)-3]
        stateVectors <- SPKArray[1:(numberRecords*6)]
        epoch <- seq(from=epoch1, by=stepSize, length.out=numberRecords)
        stateVectorsMatrix <- matrix(stateVectors, ncol=6, nrow=numberRecords, byrow=TRUE)
        colnames(stateVectorsMatrix) <- c("positionX", "positionY", "positionZ",
                                          "velocityX", "velocityY", "velocityZ")
        formattedArray <- list(
            polynomialDegree=polynomialDegree,
            windowSize=windowSize,
            stateVectors=cbind(epoch, stateVectorsMatrix)
        )
    } else if(SPKTypeCode == 13) {
        numberRecords <- SPKArray[length(SPKArray)]
        windowSize <- SPKArray[length(SPKArray) - 1] + 1
        polynomialDegree <- 2*windowSize - 1
        stateVectors <- SPKArray[1:(numberRecords*6)]
        epoch <- SPKArray[(numberRecords*6 + 1):(numberRecords*7)]
        stateVectorsMatrix <- matrix(stateVectors, ncol=6, nrow=numberRecords, byrow=TRUE)
        colnames(stateVectorsMatrix) <- c("positionX", "positionY", "positionZ",
                                          "velocityX", "velocityY", "velocityZ")
        formattedArray <- list(
            polynomialDegree=polynomialDegree,
            windowSize=windowSize,
            stateVectors=cbind(epoch, stateVectorsMatrix)
        )
    } else if(SPKTypeCode == 14) {
        numberMetaDataItems <- SPKArray[length(SPKArray)]
        if(numberMetaDataItems < 15) {
            stop("SPK kernel with missing metadata values")
        }
        if(numberMetaDataItems == 15) {
            numberMetaDataItems <- 16
        }
        metaData <- SPKArray[(length(SPKArray) - numberMetaDataItems + 1):length(SPKArray)]
        # The following should be 1 for type 14 SPKs
        numberConstants <- metaData[2]
        # Next should also be 1
        firstConstantsDouble <- metaData[1] + 1 
        numberDataPackets <- metaData[12]
        firstDataPacketsDouble <- metaData[11]
        dataPacketSize <- metaData[15]
        dataPacketSize <- dataPacketSize + 1
        numberRefValues <- metaData[7]
        firstRefValuesDouble <- metaData[6]
        polynomialDegree <- SPKArray[firstConstantsDouble]
        numberCoefficients <- polynomialDegree + 1
        # The following should be equal to dataPacketSize
        elementsPerRecord <- numberCoefficients*6 + 2
        if(elementsPerRecord != dataPacketSize) {
            stop("Inconsistency in number of coefficients and elements per record.
                 Please report the problem and provide an example file.")
        }
        lastDataPacketsDouble <- firstDataPacketsDouble + numberDataPackets * dataPacketSize
        records <- SPKArray[(firstDataPacketsDouble + 1):lastDataPacketsDouble]
        recordStartIndexes <- seq(from=0, by=elementsPerRecord, length.out=numberDataPackets)
        midPointIndexes <- 1 + recordStartIndexes
        intervalRadiiIndexes <- 2 + recordStartIndexes
        midPoints <- records[midPointIndexes]
        intervalRadii <- records[intervalRadiiIndexes]
        allCoefficients <- records[-c(midPointIndexes, intervalRadiiIndexes)]
        chebyshevCoefficients <- matrix(allCoefficients, ncol=6*numberCoefficients, nrow=numberRecords,
                                        byrow = TRUE)
        colnames(chebyshevCoefficients) <- paste(rep(c("positionXCoeff", "positionYCoeff", "positionZCoeff", 
                                                       "velocityXCoeff", "velocityYCoeff", "velocityZCoeff"), 
                                                     each=numberCoefficients), 1:numberCoefficients, sep="")
        initialEpoch <- SPKArray[(firstRefValuesDouble + 1):(firstRefValuesDouble + numberRefValues)]
        formattedArray <- list(
            polynomialDegree=polynomialDegree,
            chebyshevCoefficients=cbind(initialEpoch, midPoint, intervalRadius, chebyshevCoefficients)
        )
    } else if(SPKTypeCode == 15) {
        if(length(SPKArray != 16)) {
            stop("Malformed type 15 segment (wrong length).")
        }
        formattedArray <- as.list(SPKArray)
        names(formattedArray) <- c("epochPeriapsis", "unitVectorTrajectoryPoleX",
                                   "unitVectorTrajectoryPoleY", "unitVectorTrajectoryPoleZ",
                                   "unitVectorPeriapsisX", "unitVectorPeriapsisY", 
                                   "unitVectorPeriapsisZ", "semiLatusRectum", 
                                   "eccentricity", "J2ProcessingFlag", "unitVectorCentralBodyPoleX",
                                   "unitVectorCentralBodyPoleY", "unitVectorCentralBodyPoleZ",
                                   "centralBodyGM", "centralBodyJ2", "centralBodyRadius")
        # Note: according to SPICE comments, all units in radians, km and seconds 
        # except J2 (dimensionless)
    } else if(SPKTypeCode == 17) {
        if(length(SPKArray != 12)) {
            stop("Malformed type 17 segment (wrong length).")
        }
        formattedArray <- as.list(SPKArray)
        names(formattedArray) <- c("epochPeriapsis", "semiMajorAxis", "equinoctialH", 
                                   "equinoctialK", "meanLongitude", "equinoctialP", 
                                   "equinoctialQ", "longitudePeriapsisDerivative",
                                   "meanLongitudeDerivative", "longitudeAscendingNodeDerivative",
                                   "equatorialPoleRightAscension", "equatorialPoleDeclination")
        # All units in km, radians and radians/second.
    } else if(SPKTypeCode == 18) {
        numberDataPackets <- SPKArray[length(SPKArray)]
        windowSize <- SPKArray[length(SPKArray) - 1]
        subTypeCode <- SPKArray[length(SPKArray) - 2]
        if(subTypeCode != 0 && subTypeCode != 1) {
            stop("Invalid SPK type 18 subtype code")
        }
        if(subTypeCode == 0){
            dataPacketSize <- 12
            polynomialDegree <- 2*windowSize - 1
            elementNames <- c("positionX", "positionY", "positionZ",
                              "firstVelocityX", "firstVelocityY", "firstVelocityZ",
                              "secondVelocityX", "secondVelocityY", "secondVelocityZ",
                              "accelerationX", "accelerationY", "accelerationZ")
            interpolationType <- "Hermite"
        } else if(subTypeCode == 1) {
            dataPacketSize <- 6
            polynomialDegree <- windowSize - 1
            elementNames <- c("positionX", "positionY", "positionZ",
                              "velocityX", "velocityY", "velocityZ")
            interpolationType <- "Lagrange"
        }
        lastDataPacketsDouble <- numberDataPackets*dataPacketSize
        allEphemerides <- SPKArray[1:lastDataPacketsDouble]
        ephemeridesMatrix <- matrix(allEphemerides, ncol=dataPacketSize*allEphemerides, 
                                    nrow=numberDataPackets, byrow = TRUE)
        colnames(ephemeridesMatrix) <- elementNames
        epoch <- SPKArray[(lastDataPacketsDouble + 1):(lastDataPacketsDouble + numberDataPackets)]
        formattedArray <- list(
            subTypeCode=subTypeCode,
            polynomialDegree=polynomialDegree,
            interpolationType=interpolationType,
            windowSize=windowSize,
            stateVectors=cbind(epoch, ephemeridesMatrix)
        )
    } else if(SPKTypeCode == 19) {
        numberIntervals <- SPKArray[length(SPKArray)]
        boundaryChoiceFlag <- SPKArray[length(SPKArray) - 1]
        lastMinisegmentEndIndex <- SPKArray[length(SPKArray) - 2]
        minisegmentStartIndexes <- SPKArray[(length(SPKArray) - 2 - numberIntervals):(length(SPKArray) - 3)]
        minisegmentEndIndexes <- c(minisegmentStartIndexes[2:numberIntervals] - 1, lastMinisegmentEndIndex)
        minisegmentIndexes <- Map(`:`, minisegmentStartIndexes, minisegmentEndIndexes)
        minisegmentArrays <- lapply(minisegmentIndexes, function(i) SPKArray[i])
        intervalStarts <- SPKArray[(lastMinisegmentEndIndex + 1):(lastMinisegmentEndIndex + numberIntervals)]
        lastIntervalEnd <- SPKArray[lastMinisegmentEndIndex + numberIntervals + 1]
        intervalEnds <- c(intervalStarts[2:numberIntervals], lastIntervalEnd)
        minisegments <- vector(mode="list", length=numberIntervals)
        for(i in 1:numberIntervals) {
            newMinisegment <- formatSPKType19Minisegment(minisegmentArrays[[i]])
            append(newMinisegment, list(intervalStartEpoch=intervalStarts[i],
                                        intervalEndEpoch=intervalEnds[i]))
            minisegments[[i]] <- newMinisegment
        }
        formattedArray <- list(
            boundaryChoiceFlag=boundaryChoiceFlag,
            minisegments=minisegments
        )
    } else if(SPKTypeCode == 20) {
        numberRecords <- SPKArray[length(SPKArray)]
        elementsPerRecord <- SPKArray[length(SPKArray)-1]
        numberCoefficients <- (elementsPerRecord - 3)/3
        polynomialDegree <- numberCoefficients+1
        # we convert interval length to TDB Julian seconds
        intervalLength <- SPKArray[length(SPKArray)-2] * 86400
        # 2 records for first initial epoch in Julian TDB days, integer and fractional parts
        # we add them and convert to TDB Julian seconds
        firstInitialEpoch <- (SPKArray[length(SPKArray)-3] + SPKArray[length(SPKArray)-4]) * 86400
        # tscale is time scale used for velocity, in TDB seconds. E.g., a value of 1 means
        # we are using velocity directly in TDB seconds. A value of 86400 means the time units
        # of velocity would be TDB Julian days, etc.
        # A similar concept applies to dscale, but dscale is the scale of distance
        # units, in km. Used for both position and velocity
        tScale <- SPKArray[length(SPKArray)-5]
        dScale <- SPKArray[length(SPKArray)-6]
        records <- SPKArray[1:(length(SPKArray)-7)]
        recordStartIndexes <- seq(from=0, by=elementsPerRecord, length.out=numberRecords)
        midPointPositionXIndexes <- recordStartIndexes + 1 + numberCoefficients
        midPointPositionYIndexes <- midPointPositionXIndexes + 1 + numberCoefficients
        midPointPositionZIndexes <- midPointPositionYIndexes + 1 + numberCoefficients
        midPointPositionX <- SPKArray[midPointPositionXIndexes]
        midPointPositionY <- SPKArray[midPointPositionYIndexes]
        midPointPositionZ <- SPKArray[midPointPositionZIndexes]
        allCoefficients <- SPKArray[-c(midPointPositionXIndexes, midPointPositionYIndexes, 
                                       midPointPositionZIndexes)]
        coefficientsMatrix <- matrix(allCoefficients, ncol=3*numberCoefficients, 
                                     nrow=numberRecords, byrow = TRUE)
        colnames(coefficientsMatrix) <- paste(rep(c("velocityXCoeff", "velocityYCoeff", "velocityZCoeff"), 
                                                  each=numberCoefficients), 1:numberCoefficients, sep="")
        initialEpoch <- seq(from=firstInitialEpoch, by=intervalLength, length.out=numberRecords)
        intervalRadius <- rep(intervalLength/2, numberRecords)
        midPoint <- initialEpoch + intervalRadius
        formattedArray <- list(
            polynomialDegree=polynomialDegree,
            dScale=dScale,
            tScale=tScale,
            chebyshevCoefficients=cbind(initialEpoch, midPoint, intervalRadius, 
                                        coefficientsMatrix, midPointPositionX,
                                        midPointPositionY, midPointPositionZ)
        )
    } else if(SPKTypeCode == 21) {
        numberRecords <- SPKArray[length(SPKArray)]
        # There is 2 parameters to take into account here, DLSIZE and MAXDIM
        # DLSIZE is the total number of elements in each record (would be 71
        # for a type 1 SPK), and MAXDIM is the number of coefficients for each
        # cartesian component in the difference arrays. It would be 15 for a type
        # 1 SPK, and the relationship between the two is:
        # DLSIZE = 4*MAXDIM + 11
        # According to SPK required reading, the second-to-last element in the
        # segment (the entire array read as a generic DAF array) contains DLSIZE
        # however, comments on the python module for type 21 SPKs indicate that
        # in spite of this, it actually contains MAXDIM. So I'm going ahead with 
        # MAXDIM, but this needs verification
        maxDim <- SPKArray[length(SPKArray) - 1]
        DLSize <- 4*maxDim + 11
        formattedArray <- vector(mode="list", length=numberRecords)
        for(i in 1:numberRecords) {
            startingIndex <- DLSize*(numberRecords - 1) + 1
            finalEpochs <- SPKArray[(DLSize*numberRecords + 1):((DLSize+1)*numberRecords)]
            formattedArray[[i]] <- list(
                referenceEpoch = SPKArray[startingIndex],
                finalEpoch = finalEpochs[i],
                stepsizeFunctionVector = SPKArray[(startingIndex+1):(startingIndex+maxDim)],
                referencePosition = SPKArray[startingIndex+maxDim+c(1,3,5)],
                referenceVelocity = SPKArray[startingIndex+maxDim+c(2,4,6)],
                MDA = matrix(SPKArray[(startingIndex+maxDim+7):(startingIndex+6+maxDim*4)], ncol=3),
                numberCoefficients=maxDim,
                maxIntegrationOrderP1 = SPKArray[startingIndex+7+maxDim*4],
                integrationOrderArray = SPKArray[(startingIndex+8+maxDim*4):(startingIndex+10+maxDim*4)]
            )
        }
    } else {
        stop("Unknown SPK type")
    }
    return(formattedArray)
}

formatSPKType19Minisegment <- function(minisegmentArray) {
    numberDataPackets <- minisegmentArray[length(minisegmentArray)]
    windowSize <- minisegmentArray[length(minisegmentArray) - 1]
    subTypeCode <- minisegmentArray[length(minisegmentArray) - 2]
    if(subTypeCode != 0 && subTypeCode != 1 && subTypeCode != 1) {
        stop("Invalid SPK type 19 subtype code")
    }
    if(subTypeCode == 0){
        dataPacketSize <- 12
        polynomialDegree <- 2*windowSize - 1
        if(polynomialDegree %% 4 != 3){
            stop("Invalid interpolation degree for minisegment of subtype 0")
        }
        elementNames <- c("positionX", "positionY", "positionZ",
                          "firstVelocityX", "firstVelocityY", "firstVelocityZ",
                          "secondVelocityX", "secondVelocityY", "secondVelocityZ",
                          "accelerationX", "accelerationY", "accelerationZ")
        interpolationType <- "Hermite"
    } else if(subTypeCode == 1) {
        dataPacketSize <- 6
        polynomialDegree <- windowSize - 1
        if(polynomialDegree %% 2 != 1){
            stop("Invalid interpolation degree for minisegment of subtype 1")
        }
        elementNames <- c("positionX", "positionY", "positionZ",
                          "velocityX", "velocityY", "velocityZ")
        interpolationType <- "Lagrange"
    } else if(subTypeCode == 2) {
        dataPacketSize <- 6
        polynomialDegree <- 2*windowSize - 1
        if(polynomialDegree %% 4 != 3){
            stop("Invalid interpolation degree for minisegment of subtype 2")
        }
        elementNames <- c("positionX", "positionY", "positionZ",
                          "velocityX", "velocityY", "velocityZ")
        interpolationType <- "Hermite-joint"
    }
    lastDataPacketsDouble <- numberDataPackets*dataPacketSize
    allEphemerides <- minisegmentArray[1:lastDataPacketsDouble]
    ephemeridesMatrix <- matrix(allEphemerides, ncol=dataPacketSize*allEphemerides, 
                                nrow=numberDataPackets, byrow = TRUE)
    colnames(ephemeridesMatrix) <- elementNames
    epoch <- minisegmentArray[(lastDataPacketsDouble + 1):(lastDataPacketsDouble + numberDataPackets)]
    formattedArray <- list(
        subTypeCode=subTypeCode,
        polynomialDegree=polynomialDegree,
        interpolationType=interpolationType,
        windowSize=windowSize,
        stateVectors=cbind(epoch, ephemeridesMatrix)
    )
}

readTextKernel <- function(filename) {
    lines <- readLines(filename)
    kernelType <- trimws(lines[1])
    if(!grepl("/", kernelType) | nchar(kernelType) > 8) {
        stop("Invalid kernel type line")
    }
    labelStartLine <- grep("beginlabel", lines)
    if(length(labelStartLine) > 0) {
        labelEndLine <- grep("endlabel", lines)
        if(length(labelEndLine) == 0) {
            stop("File contains a \\beginlabel line but no \\endlabel line")
        }
        parsedLabel <- parseTextKernelLabelBlock(lines[(labelStartLine+1):(labelEndLine-1)])
        firstCommentStartLine <- labelEndLine+1
    } else {
        parsedLabel <- NULL
        firstCommentStartLine <- 2
    }
    dataBlockStartLines <- grep("begindata", lines)
    commentStartLines <- grep("begintext", lines)
    dataBlockEndLines <- commentStartLines
    commentEndLines <- c(dataBlockStartLines[-1])
    if(length(dataBlockEndLines) == length(dataBlockStartLines) - 1) {
        dataBlockEndLines <- c(dataBlockEndLines, length(lines))
        lastCommentPresent <- FALSE
    } else if(commentStartLines[length(commentStartLines)]) {
        lastCommentPresent <- TRUE
    }
    if(dataBlockStartLines[1] > firstCommentStartLine) {
        firstComment <- lines[firstCommentStartLine:(dataBlockStartLines[1]-1)]
    } else {
        firstComment <- NULL
    }
    for(i in 1:length(dataBlockStartLines)) {
        dataBlockLines <- lines[(dataBlockStartLines[i]+1):(dataBlockEndLines[1]-1)]
        dataBlockLines <- dataBlockLines[dataBlockLines!=""]
        parsedDataBlock <- parseTextKernelDataBlock(dataBlockLines)
    }
    dataLines <- lines[(dataBlockStartLines[1]+1):(dataBlockEndLines[1]-1)]
    dataLines <- dataLines[dataLines!=""]
    return(parseTextKernelDataBlock(dataLines))
}

parseTextKernelLabelBlock <- function(labelLines) {
    labelComments <- labelLines[grepl("^;", labelLines)]
    variableAssignmentStartLines <- grep("=", labelLines)
    variableAssignmentLineLengths <- diff(c(variableAssignmentStartLines, length(labelLines)))
    labelVariables <- vector(mode="list", length=length(variableAssignmentStartLines))
    labelVariableNames <- vector(mode="character", length=length(variableAssignmentStartLines))
    splitAssignments <- strsplit(labelLines, "=")
    for(i in 1:length(variableAssignmentStartLines)) {
        if(variableAssignmentLineLengths[i] == 1){
            currentAssignment <- trimws(splitAssignments[[variableAssignmentStartLines[i]]])
            variableName <- currentAssignment[1]
            variableValue <- currentAssignment[2]
        } else {
            currentAssignment <- trimws(unlist(splitAssignments[variableAssignmentStartLines[i]:(variableAssignmentStartLines[i] 
                                                                                                 + variableAssignmentLineLengths[i] -1)]))
            variableName <- currentAssignment[1]
            variableValue <- paste(currentAssignment[-1], collapse=" ")
        }
        variableValue <- gsub("\"", "", variableValue)
        labelVariables[[i]] <- variableValue
        labelVariableNames[i] <- variableName
    }
    names(labelVariables) <- labelVariableNames
    return(list(
        labelComments=labelComments,
        labelVariables=labelVariables
    ))
}

parseTextKernelDataBlock <- function(dataLines, stringContinuationMarker=NULL) {
    variableAllAssignmentStartLines <- grep("=", dataLines)
    variableAllAssignmentLineLengths <- diff(c(variableAllAssignmentStartLines, length(dataLines)+1))
    dataBlockVariables <- list()
    for(i in 1:length(variableAllAssignmentStartLines)) {
        firstAssignmentLine <- dataLines[variableAllAssignmentStartLines[i]]
        if(grepl("\\+=", firstAssignmentLine)) {
            extensionAssignment <- TRUE
            assignmentSeparator <- "\\+="
        } else {
            extensionAssignment <- FALSE
            assignmentSeparator <- "="
        }
        firstAssignmentLineSplit <- strsplit(firstAssignmentLine, assignmentSeparator)[[1]]
        variableName <- trimws(firstAssignmentLineSplit[1])
        if(variableAllAssignmentLineLengths[i] == 1){
            variableString <- trimws(firstAssignmentLineSplit[2])
        } else {
            restOfAssignment <- trimws(dataLines[(1:(variableAllAssignmentLineLengths[i]-1))+variableAllAssignmentStartLines[i]])
            #variableString <- c(trimws(firstAssignmentLineSplit[2]), restOfAssignment)
            #variableString <- paste(trimws(firstAssignmentLineSplit[2]), restOfAssignment, sep=" ")
            variableString <- paste(c(trimws(firstAssignmentLineSplit[2]), restOfAssignment), collapse =" ")
        }
        if(grepl("^\\(", variableString)) {
            # Vector of values
            # first remove brackets
            variableString <- gsub("^\\(\\s*|\\s*\\)$", "", variableString)
            # vectorCharacterValues <- strsplit(variableString, ",|\\s+")[[1]]
            # vectorCharacterValues <- vectorCharacterValues[vectorCharacterValues != ""]
            if(grepl("^'", variableString)) {
                # Variable is of type string. Even though assignment is of length 1 line
                # need to check for continued string just in case
                vectorCharacterMatches <- gregexpr("(^|(,|\\s+))\\K'(([^']*)|([^']*'{2}[^']*)*)'(?=($|(,|\\s+)))", variableString, perl = TRUE)
                vectorCharacterValues <- regmatches(variableString, vectorCharacterMatches)
                vectorCharacterValues <- gsub("^'|'$", "", vectorCharacterValues)
                vectorCharacterValues <- gsub("''", "'", vectorCharacterValues)
                if(!is.null(stringContinuationMarker)) {
                    stringContinuationMarker <- paste(stringContinuationMarker, "\\s*", sep="")
                    numberValues <- length(vectorCharacterValues - (sum(
                        grepl(stringContinuationMarker, vectorCharacterValues)
                    )))
                    variableValue <- character(numberValues)
                    variableValue[1] <- vectorCharacterValues[1]
                    continuedString <- grepl(stringContinuationMarker, vectorCharacterValues[1])
                    currentElement <- 1
                    for(i in 2:length(vectorCharacterValues)) {
                        if(continuedString) {
                            newPartOfElement <- gsub(stringContinuationMarker, "", vectorCharacterValues[i])
                            variableValue[currentElement] <- paste(variableValue[currentElement], newPartOfElement, sep="")
                        } else {
                            currentElement <- currentElement + 1
                            variableValue[currentElement] <- vectorCharacterValues[i]
                        }
                        continuedString <- grepl(stringContinuationMarker, vectorCharacterValues[i])
                    }
                }
            } else if(grepl("^@", variableString)) {
                # Variable is of type special string specifying a date-time
                # currently will be stored just a string... 
                # TODO: implement str2et equivalent ("the monster")
                vectorCharacterValues <- strsplit(variableString, ",|\\s+")[[1]]
                vectorCharacterValues <- vectorCharacterValues[vectorCharacterValues != ""]
                vectorCharacterValues <- gsub("^@", "", vectorCharacterValues)
                variableValue <- vectorCharacterValues
            } else {
                # numeric values
                vectorCharacterValues <- strsplit(variableString, ",|\\s+")[[1]]
                vectorCharacterValues <- vectorCharacterValues[vectorCharacterValues != ""]
                variableValue <- as.numeric(vectorCharacterValues)
            }
        } else {
            # not vector assignment
            if(grepl("^'", variableString)) {
                if(!is.null(stringContinuationMarker)) {
                    stringContinuationMarker <- paste(stringContinuationMarker, "\\s*", sep="")
                    variableString <- gsub(stringContinuationMarker, "", variableString)
                    variableString <- paste(variableString, collapse="")
                }
                variableValue <- gsub("^'|'$", "", variableString)
                variableValue <- gsub("''", "'", vectorCharacterValues)
            } else if(grepl("^@", variableString)) {
                # Variable is of type special string specifying a date-time
                variableValue <- gsub("^@", "", variableString)
            } else {
                #numeric value
                variableValue <- as.numeric(variableString)
            }
        }
        if(extensionAssignment) {
            dataBlockVariables[[variableName]] <- c(dataBlockVariables[[variableName]], variableValue)
        } else {
            dataBlockVariables[[variableName]] <- variableValue
        }
    }
    return(dataBlockVariables)
}

# parseTextDEHeader <- function(headerLines) {
#     headerLines <- trimws(headerLines)
#     headerLines <- headerLines[headerLines != ""]
#     firstLine <- headerLines[1]
#     if(!(grepl("KSIZE", firstLine) & grepl("NCOEFF", firstLine))) {
#         stop("First line of header lacks definitions of KSIZE or NCOEFF")
#     }
#     ksizeNcoeff <- as.numeric(unlist(regmatches(firstLine, gregexpr('\\(?[0-9,.]+', firstLine))))
#     if(length(ksizeNcoeff) != 2) {
#         stop("First line of header contains wrong number of assignments (should be 2)")
#     }
#     ksize <- ksizeNcoeff[1]
#     ncoeff <- ksizeNcoeff[2]
#     groupLineIndexes <- grep("GROUP", headerLines)
#     groupLines <- headerLinesgroupLineIndexes
#     groupNumbers <- as.numeric(unlist(regmatches(groupLines, gregexpr('\\(?[0-9,.]+', groupLines))))
#     if(identical(groupNumbers[1:5], c(1010, 1030, 1040, 1041, 1050))) {
#         stop("Header file does not contain groups 1010, 1030, 1040, 1041 and 1050 in this order")
#     }
#     group1010Index <- groupLines[1]
#     group1030Index <- groupLines[2]
#     group1040Index <- groupLines[3]
#     group1041Index <- groupLines[4]
#     group1050Index <- groupLines[5]
#     if(length(groupNumbers == 6)) {
#         group1070Index <- groupLines[6]
#     } else {
#         group1070Index <- length(headerLines)
#     }
#     group1010Lines <- headerLines[(group1010Index+1):(group1030Index-1)]
#     group1030Lines <- headerLines[(group1030Index+1):(group1040Index-1)]
#     group1040Lines <- headerLines[(group1040Index+1):(group1041Index-1)]
#     group1041Lines <- headerLines[(group1041Index+1):(group1050Index-1)]
#     group1050Lines <- headerLines[(group1050Index+1):(group1070Index-1)]
#     group1030Values <- as.numeric(unlist(strsplit(group1030Lines, "\\s+")))
#     startJD <- group1030Values[1]
#     endJD <- group1030Values[2]
#     stepSize <- group1030Values[3]
#     numberConstants <- as.numeric(group1040Lines[1])
#     namesConstants <- unlist(strsplit(group1040Lines[-1], "\\s+"))
#     valuesConstants <- as.numeric(gsub("D", "E", unlist(strsplit(h2[28:78], "\\s+")), ignore.case = TRUE))
#     
#     
#     
# }
