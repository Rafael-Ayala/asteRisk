
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

parseGLONASSNavigationRINEXlines <- function(lines, tauC=0) {
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
    clockBias <- -as.numeric(substr(line1, 23, 41))
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
    correctedEphemerisUTC <- as.nanotime(as.POSIXct(dateTimeString, tz="UTC")) + clockBias*1e9 + tauC*1e9 # AQUIAQUIAQUI
    return(list(
        satelliteNumber=satelliteNumber,
        epochYearShort=epochYearShort,
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
        sysTimeCorrection <- -as.numeric(gsub("D", "E", substr(systemTimeCorrectionLine, 22, 40)))
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
        parsedMessages[[i]] <- parseGLONASSNavigationRINEXlines(singleMessageLines, tauC=headerFields$sysTimeCorrection)
    }
    return(list(
        header=headerFields, 
        messages=parsedMessages))
}

parseGPSNavigationRINEXlines <- function(lines, leapSeconds=0, deltaUTCA0=0,
                                         deltaUTCA1=0, referenceTimeUTC,
                                         referenceWeekUTC) {
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
    tocDateTimeString <- paste(tocYear, "-", tocMonth, "-", tocDay, " ", 
                            tocHour, ":", tocMinute, ":", tocSecond, 
                            sep="")
    tocGPSTseconds <- as.double(difftime(tocDateTimeString, "1980-01-06 00:00:00", tz="UTC", units="secs")) 
    tocGPSWeek <- tocGPSTseconds%/%604800
    tocSecondsOfGPSWeek <- tocGPSTseconds%%604800
    differenceToeToc <- toeSecondsOfGPSWeek - tocSecondsOfGPSWeek
    if(differenceToeToc > 302400) differenceToeToc <- differenceToeToc - 604800
    if(differenceToeToc < -302400) differenceToeToc <- differenceToeToc + 604800
    F_constant <- -2*sqrt(GM_Earth_TCB)/c_light^2
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
    oblateEarthAccelerationFactor <- -1.5*J2_WGS84*(GM_Earth_TCB/correctedRadius^2)*(earthRadius_WGS84/correctedRadius)^2
    acceleration_ECEF <- c(
        -GM_Earth_TCB*(position_ECEF[1]/correctedRadius^3) + oblateEarthAccelerationFactor*
            ((1 - 5*(position_ECEF[3]/correctedRadius)^2)*(position_ECEF[1]/correctedRadius)) +
            2*velocity_ECEF[2]*omegaEarth + position_ECEF[1]*omegaEarth^2,
        -GM_Earth_TCB*(position_ECEF[2]/correctedRadius^3) + oblateEarthAccelerationFactor*
            ((1 - 5*(position_ECEF[3]/correctedRadius)^2)*(position_ECEF[1]/correctedRadius)) -
            2*velocity_ECEF[1]*omegaEarth + position_ECEF[2]*omegaEarth^2,
        -GM_Earth_TCB*(position_ECEF[3]/correctedRadius^3) + oblateEarthAccelerationFactor*
            ((3-5*(position_ECEF[3]/correctedRadius)^2)*(position_ECEF[3]/correctedRadius))
    )
    return(list(
        satellitePRNCode=satellitePRNCode,
        tocYearShort=tocYearShort,
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
        parsedMessages[[i]] <- parseGPSNavigationRINEXlines(singleMessageLines, headerFields$leapSeconds, 
                                                            headerFields$deltaUTCA0,
                                                            headerFields$deltaUTCA1, 
                                                            headerFields$referenceTimeUTC,
                                                            headerFields$referenceWeekUTC)
    }
    return(list(
        header=headerFields, 
        messages=parsedMessages))
}

parseOEMheader <- function(lines) {
    OEMVersion <- trimws(strsplit(lines[1], split="=")[[1]][2])
    creationDate <- gsub("T", " ", trimws(strsplit(lines[length(lines) - 1], split="=")[[1]][2]))
    creator <- trimws(strsplit(lines[length(lines)], split="=")[[1]][2])
    return(list(
        OEMVersion=OEMVersion,
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
    headerFields <- parseOEMheader(lines[1:endOfHeader])
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
