
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
    julianDay <-
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

readTLE <- function(filename) {
    lines <- readLines(filename)
    lines <- lines[lines != ""]
    # Determine if TLE are actually in 2 or 3 line format
    # if(nchar(trimws(lines[1])) <= 12) {
    if(substr(lines[1], 1, 1) != "1" | nchar(trimws(lines[1])) <= 12) {
        numberLines <- 3
    } else {
        numberLines <- 2
    }
    numberTLEs <- length(lines)/numberLines
    if(numberTLEs == 1) {
        parsedTLEs <- parseTLElines(lines)
    } else {
        startingLines <- seq(1, by=numberLines, length.out = numberTLEs)
        parsedTLEs <- vector(mode = "list", length = numberTLEs)
        for(i in 1:length(startingLines)) {
            singleTLElines <- lines[startingLines[i]:(startingLines[i]+numberLines-1)]
            parsedTLEs[[i]] <- parseTLElines(singleTLElines)
        }
    }
    return(parsedTLEs)
}

parseGLONASSNavigationRINEXlines <- function(lines) {
    if(length(lines) != 4) {
        stop("Invalid GLONASS navigation RINEX file")
    }
    line1 <- gsub("D", "E", lines[1], ignore.case=TRUE)
    line2 <- gsub("D", "E", lines[2], ignore.case=TRUE)
    line3 <- gsub("D", "E", lines[3], ignore.case=TRUE)
    line4 <- gsub("D", "E", lines[4], ignore.case=TRUE)
    satelliteNumber <- as.numeric(substr(line1, 1, 2))
    epochYearShort <- as.numeric(substr(line1, 4, 5))
    epochMonth <- as.numeric(substr(line1, 7, 8))
    epochDay <- as.numeric(substr(line1, 10, 11))
    epochHour <- as.numeric(substr(line1, 13, 14))
    epochMinute <- as.numeric(substr(line1, 16, 17))
    epochSecond <- as.numeric(substr(line1, 18, 22))
    clockBias <- as.numeric(substr(line1, 23, 41))
    clockDrift <- as.numeric(substr(line1, 42, 60))
    clockDriftRate <- as.numeric(substr(line1, 61, 79))
    relativeFreqBias <- as.numeric(substr(line1, 42, 60))
    messageFrameTime <- as.numeric(substr(line1, 61, 79))
    positionX <- as.numeric(substr(line2, 4, 22))
    velocityX <- as.numeric(substr(line2, 23, 41))
    accelX <- as.numeric(substr(line2, 42, 60))
    satelliteHealthCode <- as.numeric(substr(line2, 61, 79))
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
    dateTimeString <- paste(epochYear, "-", epochMonth, "-", epochDay, " ", 
                            epochHour, ":", epochMinute, ":", epochSecond, 
                            sep="")
    dateTimePOSIXct <- as.POSIXct(dateTimeString)
    UTCDateTimePOSIXct <- dateTimePOSIXct # - 3*3600
    UTCDateTimeString <- as.character(UTCDateTimePOSIXct) # In UTC
    return(list(
        satelliteNumber=satelliteNumber,
        epochYearShort=epochYearShort,
        epochMonth=epochMonth,
        epochDay=epochDay,
        epochHour=epochHour,
        epochMinute=epochMinute,
        epochSecond=epochSecond,
        UTCepochDateTime=UTCDateTimeString,  # TODO consider clock bias difference GLONASS-UTC
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
        satelliteHealthCode=satelliteHealthCode,
        freqNumber=freqNumber,
        informationAge=informationAge
    ))
}

parseGLONASSNavigationRINEXheaderLines <- function(lines) {
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
        sysTimeCorrection <- as.numeric(gsub("D", "E", substr(systemTimeCorrectionLine, 22, 40)))
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

readGLONASSNavigationRINEX <- function(filename) {
    lines <- readLines(filename)
    lines <- lines[lines != ""]
    endOfHeader <- grep("END OF HEADER", lines)
    headerLines <- lines[1:endOfHeader]
    bodyLines <- lines[(endOfHeader+1):length(lines)]
    headerFields <- parseGLONASSNavigationRINEXheaderLines(headerLines)
    messageNumberLines <- 4
    numberMessages <- length(bodyLines)/messageNumberLines
    parsedMessages <- vector(mode = "list", length = numberMessages)
    startingLines <- seq(1, by=messageNumberLines, length.out = numberMessages)
    for(i in 1:length(startingLines)) {
        singleMessageLines <- bodyLines[startingLines[i]:(startingLines[i]+messageNumberLines-1)]
        parsedMessages[[i]] <- parseGLONASSNavigationRINEXlines(singleMessageLines)
    }
    return(list(
        header=headerFields, 
        messages=parsedMessages))
}

parseGPSNavigationRINEXlines <- function(lines, leapSeconds=0) {
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
    satellitePRNCode <- as.numeric(substr(line1, 1, 2))
    epochYearShort <- as.numeric(substr(line1, 4, 5))
    epochMonth <- as.numeric(substr(line1, 7, 8))
    epochDay <- as.numeric(substr(line1, 10, 11))
    epochHour <- as.numeric(substr(line1, 13, 14))
    epochMinute <- as.numeric(substr(line1, 16, 17))
    epochSecond <- as.numeric(substr(line1, 18, 22))
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
    timeOfEphemeris <- as.numeric(substr(line4, 4, 22))
    inclinationCorrectionCosine <- as.numeric(substr(line4, 23, 41))
    ascension <- as.numeric(substr(line4, 42, 60)) # radians, OMEGA angle
    inclinationCorrectionSine <- as.numeric(substr(line4, 61, 79))
    inclination <- as.numeric(substr(line5, 4, 22)) # radiands, initial inclination
    radiusCorrectionCosine <- as.numeric(substr(line5, 23, 41))
    perigeeArgument <- as.numeric(substr(line5, 42, 60)) # radians, omega angle
    OMEGADot <- as.numeric(substr(line5, 61, 79)) # radians/s, derivtive of OMEGA
    inclinationRate <- as.numeric(substr(line6, 4, 22)) 
    codesL2Channel <- as.numeric(substr(line6, 23, 41))
    GPSWeek <- as.numeric(substr(line6, 42, 60)) 
    dataFlagL2P <- as.numeric(substr(line6, 61, 79))
    satelliteAccuracy <- as.numeric(substr(line7, 4, 22)) 
    satelliteHealthCode <- as.numeric(substr(line7, 23, 41))
    totalGroupDelay <- as.numeric(substr(line7, 42, 60)) 
    IODC <- as.numeric(substr(line7, 61, 79))
    transmissionTime <- as.numeric(substr(line8, 4, 22)) # seconds of GPS week
    fitInterval <- as.numeric(substr(line8, 23, 41))
    if (is.numeric(epochYearShort)) if(epochYearShort >= 80) {
        epochYear <- 1900 + epochYearShort
    } else {
        epochYear <- 2000 + epochYearShort
    }
    dateTimeString <- paste(epochYear, "-", epochMonth, "-", epochDay, " ", 
                            epochHour, ":", epochMinute, ":", epochSecond, 
                            sep="")
    dateTimePOSIXct <- as.POSIXct(dateTimeString)
    UTCDateTimePOSIXct <- dateTimePOSIXct - leapSeconds
    UTCDateTimeString <- as.character(UTCDateTimePOSIXct) # In UTC
    return(list(
        satellitePRNCode=satellitePRNCode,
        epochYearShort=epochYearShort,
        epochMonth=epochMonth,
        epochDay=epochDay,
        epochHour=epochHour,
        epochMinute=epochMinute,
        epochSecond=epochSecond,
        UTCepochDateTime=UTCDateTimeString,  # TODO consider clock bias difference GPS-UTC
        clockBias=clockBias,
        clockDrift=clockDrift,
        clockDriftRate=clockDriftRate,
        IODE=IODE,
        radiusCorrectionSine=radiusCorrectionSine,
        deltaN=deltaN,
        meanAnomaly=meanAnomaly,
        latitudeCorrectionCosine=latitudeCorrectionCosine,
        eccentricity=eccentricity,
        latitudeCorrectionSine=latitudeCorrectionSine,
        semiMajorAxis=semiMajorAxis,
        timeOfEphemeris=timeOfEphemeris,
        inclinationCorrectionCosine=inclinationCorrectionCosine,
        ascension=ascension,
        inclinationCorrectionSine=inclinationCorrectionSine,
        inclination=inclination,
        radiusCorrectionCosine=radiusCorrectionCosine,
        perigeeArgument=perigeeArgument,
        OMEGADot=OMEGADot,
        inclinationRate=inclinationRate,
        codesL2Channel=codesL2Channel,
        GPSWeek=GPSWeek,
        dataFlagL2P=dataFlagL2P,
        satelliteAccuracy=satelliteAccuracy,
        satelliteHealthCode=satelliteHealthCode,
        totalGroupDelay=totalGroupDelay,
        IODC=IODC,
        transmissionTime=transmissionTime,
        fitInterval=fitInterval
    ))
}

parseGPSNavigationRINEXheaderLines <- function(lines) {
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
    deltaUTCLineIndex <- grep("DELTA UTC", lines)
    if(length(deltaUTCLineIndex) > 0) {
        deltaUTCLine <- lines[deltaUTCLineIndex]
        deltaUTCA0 <- as.numeric(gsub("D", "E", substr(deltaUTCLine, 1, 21)))
        deltaUTCA1 <- as.numeric(gsub("D", "E", substr(deltaUTCLine, 22, 42)))
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

readGPSNavigationRINEX <- function(filename) {
    lines <- readLines(filename)
    lines <- lines[lines != ""]
    endOfHeader <- grep("END OF HEADER", lines)
    headerLines <- lines[1:endOfHeader]
    bodyLines <- lines[(endOfHeader+1):length(lines)]
    headerFields <- parseGPSNavigationRINEXheaderLines(headerLines)
    messageNumberLines <- 8
    numberMessages <- length(bodyLines)/messageNumberLines
    parsedMessages <- vector(mode = "list", length = numberMessages)
    startingLines <- seq(1, by=messageNumberLines, length.out = numberMessages)
    for(i in 1:length(startingLines)) {
        singleMessageLines <- bodyLines[startingLines[i]:(startingLines[i]+messageNumberLines-1)]
        parsedMessages[[i]] <- parseGPSNavigationRINEXlines(singleMessageLines, headerFields$leapSeconds)
    }
    return(list(
        header=headerFields, 
        messages=parsedMessages))
}
