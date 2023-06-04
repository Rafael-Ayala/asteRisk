dateTimeToMJD <- function(dateTime, timeSystem="UTC") {
    if(!(timeSystem %in% c("UTC", "UT1", "TT", "TDB"))) {
        stop(strwrap("Please choose a time system from UTC, UT1, TT and TDB"))
    }
    nanotimeDateTime <- nanotime(dateTime, format="%Y-%m-%d %H:%M:%E9S", tz = "UTC")
    if(is.na(nanotimeDateTime)) {
        stop("Please provide a complete date-time string in 
             year-month-day hour:minute:second format.")
    }
    dateTimeComponents <- unclass(as.POSIXlt(nanotimeDateTime, tz="UTC"))
    year <- dateTimeComponents$year + 1900
    month <- dateTimeComponents$mon + 1
    day <- dateTimeComponents$mday
    hour <- dateTimeComponents$hour
    minute <- dateTimeComponents$min
    second <- dateTimeComponents$sec
    MJD <- iauCal2jd(year, month, day, hour, minute, second)$DATE
    if(timeSystem == "TT") {
        hasData()
        MJD <- MJDUTCtoMJDTT(MJD)
    } else if(timeSystem == "UT1") {
        hasData()
        MJD <- MJDUTCtoMJDUT1(MJD)
    } else if(timeSystem == "TDB") {
        hasData()
        MJD <- Mjday_TDB(MJDUTCtoMJDTT(MJD))
    }
    return(MJD)
}

delaunayVariables <- function(CenturiesJ2000_TDB, degreesOutput=FALSE) {
    rev <- if(degreesOutput) 360 else (2*pi)
    delaunayL <- 134.96340251 + (CenturiesJ2000_TDB * (1717915923.2178 + 
                    CenturiesJ2000_TDB * (31.8792 + CenturiesJ2000_TDB *
                    (0.051635 + CenturiesJ2000_TDB * - 0.00024470))))/3600
    delaunayLprime <- 357.52910918 + (CenturiesJ2000_TDB * (129596581.0481 + 
                     CenturiesJ2000_TDB * (-0.5532 + CenturiesJ2000_TDB *
                     (0.000136 + CenturiesJ2000_TDB * - 0.00001149))))/3600
    delaunayF <- 93.27209062 + (CenturiesJ2000_TDB * (1739527262.8478 + 
                   CenturiesJ2000_TDB * (-12.7512 + CenturiesJ2000_TDB *
                   (-0.001037 + CenturiesJ2000_TDB * 0.00000417))))/3600
    delaunayD <- 297.85019547 + (CenturiesJ2000_TDB * (1602961601.2090 + 
                   CenturiesJ2000_TDB * (-6.3706 + CenturiesJ2000_TDB *
                   (0.006593 + CenturiesJ2000_TDB * -0.00003169))))/3600
    delaunayOmega <- 125.04455501 + (CenturiesJ2000_TDB * (-6962890.5431 + 
                     CenturiesJ2000_TDB * (7.4722 + CenturiesJ2000_TDB *
                     (0.007702 + CenturiesJ2000_TDB * -0.00005939))))/3600
    delaunayValues <- c(
        delaunayL,
        delaunayLprime,
        delaunayF,
        delaunayD,
        delaunayOmega
    )
    if(!degreesOutput) delaunayValues <- deg2rad(delaunayValues)
    #delaunayValues <- rem(delaunayValues, rev)
    return(delaunayValues)
}

delaunayToDoodsonVariables <- function(delaunayVariables, GMST, degreesInput=FALSE) {
    halfRev <- if(degreesInput) 180 else pi
    doodsonS <- delaunayVariables[3] + delaunayVariables[5]
    doodsonH <- doodsonS - delaunayVariables[4]
    doodsonP <- doodsonS - delaunayVariables[1]
    doodsonNprime <- -delaunayVariables[5]
    doodsonPs <- doodsonS - delaunayVariables[4] - delaunayVariables[2]
    doodsonTau <- GMST + halfRev - doodsonS
    doodsonVars <- c(doodsonTau, doodsonS, doodsonH, doodsonP, doodsonNprime, doodsonPs)
    return(doodsonVars)
}

xprod <- function(...) {
    ## Code by M. Lundberg and G.V.Welland
    args <- list(...)
    # Check for valid arguments
    if (length(args) == 0) {
        stop("No data supplied")
    }
    len <- unique(sapply(args, FUN=length))
    if (length(len) > 1) {
        stop("All vectors must be the same length")
    }
    if (len != length(args) + 1) {
        stop("Must supply N-1 vectors of length N")
    }
    # Compute generalized cross product by taking the determinant of sub-matrices
    m <- do.call(rbind, args)
    sapply(seq(len),
           FUN=function(i) {
               det(m[,-i,drop=FALSE]) * (-1)^(i+1)
           })
}

JPLephemerides <- function(MJD, timeSystem="UTC", centralBody="SSB", derivatives="acceleration") {
    hasData()
    if(!(timeSystem %in% c("UTC", "UT1", "TT", "TDB"))) {
        stop(strwrap("Please choose a time system from UTC, UT1, TT and TDB"))
    }
    if(!(derivatives %in% c("none", "velocity", "acceleration"))) {
        stop(strwrap("Please choose a derivatives level from none, velocities or acceleration"))
    }
    if(length(centralBody) != 1 | !(centralBody %in% c("SSB", "Mercury", "Venus", "Earth",
                                                       "Moon", "Jupiter", "Saturn", "Uranus",
                                                       "Neptune", "Pluto")))
    if(timeSystem == "UTC") {
        MJD <- Mjday_TDB(MJDUTCtoMJDTT(MJD))
    } else if(timeSystem == "UT1") {
        MJD <- Mjday_TDB(MJDUTCtoMJDTT(MJDUT1toMJDUTC(MJD)))
    } else if(timeSystem == "TT") {
        MJD <- Mjday_TDB(MJD)
    }
    derivativesOrder <- switch(derivatives, velocity=1, acceleration=2, 0)
    return(JPLephemeridesDE440(MJD, centralBody=centralBody, derivativesOrder=derivativesOrder))
}

# doodsonVariables <- function(MJD_TT) {
#     Tu0 <- (floor(MJD_TT)-51544.5)/36525.0
#     gmst0 <- (6.0/24 + 41.0/(24*60) + 50.54841/(24*60*60))
#     gmst0 <- gmst0 + (8640184.812866/(24*60*60))*Tu0
#     gmst0 <- gmst0 + (0.093104/(24*60*60))*Tu0*Tu0
#     gmst0 <- gmst0 + (-6.2e-6/(24*60*60))*Tu0*Tu0*Tu0
#     r <- 1.002737909350795 + 5.9006e-11*Tu0 - 5.9e-15*Tu0*Tu0
#     gmst <- rem(2*pi*(gmst0 + r * rem(MJD_TT,1)), 2*pi)
#     t = (MJD_TT-51544.5)/365250.0
#     doodsonS = (218.31664562999 + (4812678.81195750 + (-0.14663889 + ( 0.00185140 + -0.00015355*t)*t)*t)*t) * pi/180
#     doodsonH = (280.46645016002 + ( 360007.69748806 + ( 0.03032222 + ( 0.00002000 + -0.00006532*t)*t)*t)*t) * pi/180
#     doodsonP = ( 83.35324311998 + (  40690.13635250 + (-1.03217222 + (-0.01249168 +  0.00052655*t)*t)*t)*t) * pi/180
#     doodsonNprime = (234.95544499000 + (  19341.36261972 + (-0.20756111 + (-0.00213942 +  0.00016501*t)*t)*t)*t) * pi/180
#     doodsonPs = (282.93734098001 + (     17.19457667 + ( 0.04568889 + (-0.00001776 + -0.00003323*t)*t)*t)*t) * pi/180
#     doodsonTau = gmst + pi - doodsonS
#     doodsonVars <- c(doodsonTau, doodsonS, doodsonH, doodsonP, doodsonNprime, doodsonPs)
#     return(doodsonVars)
# }


# alternativeGmst <- function(MJD_TT) {
#     Tu0 <- (floor(MJD_TT)-51544.5)/36525.0
#     gmst0 <- (6.0/24 + 41.0/(24*60) + 50.54841/(24*60*60))
#     gmst0 <- gmst0 + (8640184.812866/(24*60*60))*Tu0
#     gmst0 <- gmst0 + (0.093104/(24*60*60))*Tu0*Tu0
#     gmst0 <- gmst0 + (-6.2e-6/(24*60*60))*Tu0*Tu0*Tu0
#     r <- 1.002737909350795 + 5.9006e-11*Tu0 - 5.9e-15*Tu0*Tu0
#     gmst <- rem(2*pi*(gmst0 + r * rem(MJD_TT,1)), 2*pi)
#     return(gmst)
# }

legendreNormFactor <- function(n, m) {
    if(m==0) {
        delta <- 1
    } else {
        delta <- 0
    }
    sqrt((factorial(n-m)*(2*n+1)*(2-delta))/factorial(n+m))
}

acot <- function(x) {
    return(atan(1/x))
}
