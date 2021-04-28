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

vectorCrossProduct3D <- function(u, v) {
    w <- c(u[2]*v[3] - u[3]*v[2],
           u[3]*v[1] - u[1]*v[3],
           u[1]*v[2] - u[2]*v[1])
    return(w)
}
