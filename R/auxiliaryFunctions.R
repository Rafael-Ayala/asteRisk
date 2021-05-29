deg2rad <- function(degrees) {
    return(degrees*pi/180)
}

rad2deg <- function(radians) {
    return(radians*180/pi)
}

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
    } else if(launchYearShort >= 56) {
        launchYear <- 1900 + launchYearShort
    } else {
        launchYear <- 2000 + launchYearShort
    }
    launchNumber = substr(line1, 12, 14)
    launchPiece = gsub(" ", "", substr(line1, 15, 17), fixed=TRUE)
    epochYearShort <- as.numeric(substr(line1, 19, 20))
    if (is.numeric(epochYearShort)) if(epochYearShort >= 56) {
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

UTCdateTimeToGMST <- function(dateTime) {
    daysToJ2000_0 <- as.numeric(julian(as.POSIXct(dateTime, tz="UTC"),
                                       origin=as.POSIXct("2000-01-01 12:00:00", tz="UTC")))
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
    a <- (earth_mu)^(1/3) / ((2*pi*meanMotion/86400)^(2/3))
    return(a)
}

semiMajorAxistToMeanMotion <- function(semiMajorAxis) {
    n <- sqrt(earth_mu/semiMajorAxis^3) * (86400/(2*pi))
    return(n)
}

vectorCrossProduct3D <- function(u, v) {
    w <- c(u[2]*v[3] - u[3]*v[2],
           u[3]*v[1] - u[1]*v[3],
           u[1]*v[2] - u[2]*v[1])
    return(w)
}

getLatestSpaceData <- function() { # TODO: MOVE TO THIS PACKAGE
    hasData()
    asteRiskData::getLatestSpaceData()
}

parseNavigationRINEXlinesGLONASSv21 <- function(lines) {
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
    relativeFreqBias <- as.numeric(substr(line1, 42, 60))
    messageTimeFrame <- as.numeric(substr(line1, 61, 79))
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
    if (is.numeric(epochYearShort)) if(epochYearShort >= 56) {
        epochYear <- 1900 + epochYearShort
    } else {
        epochYear <- 2000 + epochYearShort
    }
    dateTimeString <- paste(epochYear, "-", epochMonth, "-", epochDay, " ", 
                            epochHour, ":", epochMinute, ":", epochSecond, 
                            sep="")
    dateTimePOSIXct <- as.POSIXct(dateTimeString)
    UTCDateTimePOSIXct <- dateTimePOSIXct - 3*3600
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
        messageTimeFrame=messageTimeFrame,
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

parseNavigationRINEXheaderLines <- function(lines) {
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

readNavigationRINEXfile <- function(filename, messageNumberLines, rinexParser) {
    lines <- readLines(filename)
    lines <- lines[lines != ""]
    endOfHeader <- grep("END OF HEADER", lines)
    headerLines <- lines[1:endOfHeader]
    bodyLines <- lines[(endOfHeader+1):length(lines)]
    headerFields <- parseNavigationRINEXheaderLines(headerLines)
    numberMessages <- length(bodyLines)/messageNumberLines
    parsedMessages <- vector(mode = "list", length = numberMessages)
    startingLines <- seq(1, by=messageNumberLines, length.out = numberMessages)
    for(i in 1:length(startingLines)) {
        singleMessageLines <- bodyLines[startingLines[i]:(startingLines[i]+messageNumberLines-1)]
        parsedMessages[[i]] <- rinexParser(singleMessageLines)
    }
    return(list(
        header=headerFields, 
        messages=parsedMessages))
}


# readRINEXfileGLONASSv21 <- function(filename) {
#     lines <- readLines(filename)
#     lines <- lines[lines != ""]
#     endOfHeader <- grep("END OF HEADER", lines)
#     headerLines <- lines[1:endOfHeader]
#     bodyLines <- lines[(endOfHeader+1):length(lines)]
#     headerFields <- parseRINEXheaderLines(headerLines)
#     numberMessages <- length(bodyLines)/4
#     if(numberTLEs == 1) {
#         parsedTLEs <- parseTLElines(lines)
#     } else {
#         startingLines <- seq(1, by=numberLines, length.out = numberTLEs)
#         parsedTLEs <- vector(mode = "list", length = numberTLEs)
#         for(i in 1:length(startingLines)) {
#             singleTLElines <- lines[startingLines[i]:(startingLines[i]+numberLines-1)]
#             parsedTLEs[[i]] <- parseTLElines(singleTLElines)
#         }
#     }
# }