
hpop_debug <- function(position, velocity, dateTime, times, satelliteMass, dragArea, 
                       radiationArea, dragCoefficient, radiationCoefficient, 
                       earthSphericalHarmonicsDegree = 130, solidEarthTides=TRUE,
                       oceanTides=TRUE, moonSphericalHarmonicsDegree = 150, 
                       solidMoonTides=TRUE, centralBody="Earth", autoCentralBodyChange=TRUE, ...) {
    if(!(centralBody %in% c("SSB", "Mercury", "Venus", "Earth", "Moon",
                            "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto")))
        hasData()
    propagationGlobalVars <- new.env()
    extraArgs <- list(...)
    UTCSecondsJ2000 <- as.numeric(difftime(as.POSIXct(dateTime, tz="UTC"), J2000_POSIXct, units="secs"))
    MJD_UTC <- UTCSecondsJ2000/86400 + MJD_J2000
    MJD_trunc <- trunc(MJD_UTC)
    if(MJD_trunc > asteRiskData::earthPositions[nrow(asteRiskData::earthPositions), 4]) {
        message(strwrap("Attempting propagation for a date for which Earth orientation
                         parameters and other space data are not currently available. 
                         Now running getLatestSpaceData() to try to retrieve updated data."), 
                initial="", prefix="\n")
        getLatestSpaceData()
        if(MJD_trunc > asteRiskData::earthPositions[nrow(asteRiskData::earthPositions), 4]) {
            stop(strwrap("The required data was not found even after updating. The
                            target propagation date is too long into the future, and
                            therefore cannot be performed."), 
                 initial="", prefix="\n")
        } else {
            message(strwrap("Required data succesfully retrieved. Proceeding with
                            calculation of trajectory."), 
                    initial="", prefix="\n")
        }
    }
    MJD_TDB <- Mjday_TDB(MJDUTCtoMJDTT(MJD_UTC))
    JPLephemerides <- JPLephemeridesDE440(MJD_TDB, centralBody=centralBody, derivativesOrder=2)
    realCentralBody <- determineCentralBody(position, JPLephemerides[3:11], JPLephemerides[[12]])
    if(centralBody != realCentralBody) {
        message(strwrap(paste(centralBody, " was selected as the central body,
                              but the object has been determined to be in the sphere of
                              influence of ", realCentralBody, ". Propagation will therefore
                              be performed with the latter as central body.", sep=""), initial="", prefix="\n"))
        position <- position - JPLephemerides[[paste("position", realCentralBody, sep="")]]
        velocity <- velocity - JPLephemerides[[paste("velocity", realCentralBody, sep="")]]
    }
    #propagationGlobalVars$latestCentralBody <- realCentralBody
    latestCentralBody <- realCentralBody
    timeOfLastCentralBodyChange <- 0
    initial_state <- c(position, velocity)
    progressBar <- txtProgressBar(min = 0, max = max(times))
    parameters = list(
        MJD_UTC = MJD_UTC,
        MJD_TDB = MJD_TDB,
        solarArea = radiationArea,
        satelliteMass = satelliteMass,
        satelliteArea = dragArea,
        Cr = radiationCoefficient,
        Cd = dragCoefficient,
        earthSPHDegree = earthSphericalHarmonicsDegree,
        moonSPHDegree = moonSphericalHarmonicsDegree,
        SETcorrections = solidEarthTides,
        OTcorrections = oceanTides,
        SMTcorrections = solidMoonTides,
        centralBody = realCentralBody,
        autoCentralBodyChange = autoCentralBodyChange,
        progressBar = progressBar,
        globalVarsEnv = environment())
    integration_results <- ode(y=initial_state, times=times, func=odeModel_debug,
                               parms=parameters, 
                               #method="radau", rtol=1e-13,
                               #atol=1e-16, hini=0.01, 
                               ...)
    close(progressBar)
    previousStepCentralBodies <- integration_results[, 8]
    oldCentralBody <- previousStepCentralBodies[1]
    if(autoCentralBodyChange & length(unique(previousStepCentralBodies)) > 1) {
        combinedResults <- integration_results
        totalChangePoint <- 0
        while(length(unique(previousStepCentralBodies)) > 1) {
            changePoint <- which(diff(previousStepCentralBodies) != 0)[1] + 1
            newCentralBody <- previousStepCentralBodies[changePoint]
            newTimes <- times[changePoint:length(times)] - times[changePoint]
            newMjd_UTC <- MJD_UTC + times[changePoint]/86400
            newMJD_TDB <- Mjday_TDB(MJDUTCtoMJDTT(newMjd_UTC))
            JPLephemerides_oldCentralBody <- JPLephemeridesDE440(newMJD_TDB, centralBody=names(centralBodiesNum[oldCentralBody]), 
                                                                 derivativesOrder=2)
            newPosition <- integration_results[changePoint, 2:4] -
                JPLephemerides_oldCentralBody[[paste("position", names(centralBodiesNum[newCentralBody]), sep="")]]
            newVelocity <- integration_results[changePoint, 5:7] -
                JPLephemerides_oldCentralBody[[paste("velocity", names(centralBodiesNum[newCentralBody]), sep="")]]
            newInitial_state <- c(newPosition, newVelocity)
            progressBar <- txtProgressBar(min = 0, max = max(newTimes))
            newParameters <- list(
                MJD_UTC = newMjd_UTC,
                MJD_TDB = newMJD_TDB,
                solarArea = radiationArea,
                satelliteMass = satelliteMass,
                satelliteArea = dragArea,
                Cr = radiationCoefficient,
                Cd = dragCoefficient,
                earthSPHDegree = earthSphericalHarmonicsDegree,
                moonSPHDegree = moonSphericalHarmonicsDegree,
                SETcorrections = solidEarthTides,
                OTcorrections = oceanTides,
                SMTcorrections = solidMoonTides,
                centralBody = names(which(centralBodiesNum == newCentralBody)),
                autoCentralBodyChange = autoCentralBodyChange,
                progressBar = progressBar,
                globalVarsEnv = environment()
            )
            newIntegration_results <- ode(y=newInitial_state, times=newTimes, func=odeModel,
                                          parms=newParameters, 
                                          #method="radau", rtol=1e-13,
                                          #atol=1e-16, hini=0.01, 
                                          ...)
            close(progressBar)
            totalChangePoint <- totalChangePoint + changePoint 
            combinedResults[totalChangePoint:nrow(combinedResults), ] <- newIntegration_results
            previousStepCentralBodies <- newIntegration_results[, 8]
            oldCentralBody <- previousStepCentralBodies[1]
        }
        combinedResults[, 1] <- times
        integration_results <- combinedResults
    }
    numeric_results <- integration_results[, c(1:7, 9:164)]
    central_bodies <- names(centralBodiesNum[integration_results[, 8]])
    output <- cbind(as.data.frame(numeric_results), central_bodies)
    colnames(output) <- c("time", "positionX", "positionY", "positionZ", 
                          "velocityX", "velocityY", "velocityZ", 
                          "accelerationX", "accelerationY", "accelerationZ",
                          "accelMoonSPHX", "accelMoonSPHY", "accelMoonSPHZ",
                          "accelEarthSPHX", "accelEarthSPHY", "accelEarthSPHZ",
                          "accelSolarRadX", "accelSolarRadY", "accelSolarRadZ",
                          "accelSunPointX", "accelSunPointY", "accelSunPointZ",
                          "accelMercuryPointX", "accelMercuryPointY", "accelMercuryPointZ",
                          "accelVenusPointX", "accelVenusPointY", "accelVenusPointZ",
                          "accelMarsPointX", "accelMarsPointY", "accelMarsPointZ",
                          "accelJupiterPointX", "accelJupiterPointY", "accelJupiterPointZ",
                          "accelSaturnPointX", "accelSaturnPointY", "accelSaturnPointZ",
                          "accelUranusPointX", "accelUranusPointY", "accelUranusPointZ",
                          "accelNeptunePointX", "accelNeptunePointY", "accelNeptunePointZ",
                          "accelPlutoPointX", "accelPlutoPointY", "accelPlutoPointZ",
                          "positionEarthX", "positionEarthY", "positionEarthZ",
                          "velocityEarthX", "velocityEarthY", "velocityEarthZ",
                          "accelerationEarthX", "accelerationEarthY", "accelerationEarthZ",
                          "positionMoonX", "positionMoonY", "positionMoonZ",
                          "velocityMoonX", "velocityMoonY", "velocityMoonZ",
                          "accelerationMoonX", "accelerationMoonY", "accelerationMoonZ",
                          "positionMercuryX", "positionMercuryY", "positionMercuryZ",
                          "velocityMercuryX", "velocityMercuryY", "velocityMercuryZ",
                          "accelerationMercuryX", "accelerationMercuryY", "accelerationMercuryZ",
                          "positionVenusX", "positionVenusY", "positionVenusZ",
                          "velocityVenusX", "velocityVenusY", "velocityVenusZ",
                          "accelerationVenusX", "accelerationVenusY", "accelerationVenusZ",
                          "positionMarsX", "positionMarsY", "positionMarsZ",
                          "velocityMarsX", "velocityMarsY", "velocityMarsZ",
                          "accelerationMarsX", "accelerationMarsY", "accelerationMarsZ",
                          "positionSaturnX", "positionSaturnY", "positionSaturnZ",
                          "velocitySaturnX", "velocitySaturnY", "velocitySaturnZ",
                          "accelerationSaturnX", "accelerationSaturnY", "accelerationSaturnZ",
                          "positionJupiterX", "positionJupiterY", "positionJupiterZ",
                          "velocityJupiterX", "velocityJupiterY", "velocityJupiterZ",
                          "accelerationJupiterX", "accelerationJupiterY", "accelerationJupiterZ",
                          "positionNeptuneX", "positionNeptuneY", "positionNeptuneZ",
                          "velocityNeptuneX", "velocityNeptuneY", "velocityNeptuneZ",
                          "accelerationNeptuneX", "accelerationNeptuneY", "accelerationNeptuneZ",
                          "positionUranusX", "positionUranusY", "positionUranusZ",
                          "velocityUranusX", "velocityUranusY", "velocityUranusZ",
                          "accelerationUranusX", "accelerationUranusY", "accelerationUranusZ",
                          "positionPlutoX", "positionPlutoY", "positionPlutoZ",
                          "velocityPlutoX", "velocityPlutoY", "velocityPlutoZ",
                          "accelerationPlutoX", "accelerationPlutoY", "accelerationPlutoZ",
                          "positionSunX", "positionSunY", "positionSunZ",
                          "velocitySunX", "velocitySunY", "velocitySunZ",
                          "accelerationSunX", "accelerationSunY", "accelerationSunZ",
                          "positionSunSSBX", "positionSunSSBY", "positionSunSSBZ",
                          "velocitySunSSBX", "velocitySunSSBY", "velocitySunSSBZ",
                          "accelerationSunSSBX", "accelerationSunSSBY", "accelerationSunSSBZ",
                          "lunarLibrationPhi", "lunarLibrationTheta", "lunarLibrationPsi",
                          "lunarLibrationPhiPrime", "lunarLibrationThetaPrime", "lunarLibrationPsiPrime",
                          "lunarLibrationPhiPrime2", "lunarLibrationThetaPrime2", "lunarLibrationPsiPrime2",
                          "Central body")
    return(output)
}


odeModel_debug <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        state_vector <- state
        results <- accel_debug(t, state_vector, MJD_UTC, solarArea, satelliteMass, 
                         satelliteArea, Cr, Cd, earthSPHDegree, SETcorrections,
                         OTcorrections, moonSPHDegree, SMTcorrections, 
                         centralBody, autoCentralBodyChange)
        centralBody <- results[[2]]
        if(autoCentralBodyChange & centralBody != globalVarsEnv$latestCentralBody) {
            if((t - globalVarsEnv$timeOfLastCentralBodyChange) > 3600) {
                message(strwrap(paste("A transition from the sphere of influence of ",
                                      globalVarsEnv$latestCentralBody, " to that of ",
                                      centralBody, " has been detected. Therefore, integration
                                  and results from this point will be in ICRF centered on
                                  the new central body.", sep=""), initial="", prefix="\n"))
                assign("latestCentralBody", centralBody, envir = globalVarsEnv)
                assign("timeOfLastCentralBodyChange", t, envir = globalVarsEnv)
            }
        }
        acceleration <- results[[1]]
        dx <- acceleration[1, 1]
        dy <- acceleration[1, 2]
        dz <- acceleration[1, 3]
        d2x <- acceleration[2, 1]
        d2y <- acceleration[2, 2]
        d2z <- acceleration[2, 3]
        setTxtProgressBar(progressBar, t)
        list(c(dx, dy, dz, d2x, d2y, d2z),
             centralBodiesNum[centralBody], c(d2x, d2y, d2z, unlist(results[3:53])))
    })
}


accel_debug <- function(t, Y, MJD_UTC, solarArea, satelliteMass, satelliteArea, Cr, Cd, 
                  earthSPHDegree, SETcorrections, OTcorrections, moonSPHDegree, SMTcorrections,
                  centralBody, autoCentralBodyChange) {
    MJD_UTC <- MJD_UTC + t/86400
    MJD_TT <- MJDUTCtoMJDTT(MJD_UTC)
    MJD_TDB <- MJDUTCtoMJDTDB(MJD_UTC)
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, "l")
    x_pole <- IERS_results$x_pole[[1]]
    y_pole <- IERS_results$y_pole[[1]]
    UT1_UTC <- IERS_results$UT1_UTC[[1]]
    MJD_UT1 <- MJD_UTC + UT1_UTC/86400
    LOD <- IERS_results$LOD[[1]]
    dpsi <- IERS_results$dpsi[[1]]
    deps <- IERS_results$deps[[1]]
    dx_pole <- IERS_results$dx_pole[[1]]
    dy_pole <- IERS_results$dy_pole[[1]]
    TAI_UTC <- IERS_results$TAI_UTC[[1]]
    timeDiffs_results <- timeDiffs(UT1_UTC, TAI_UTC)
    UT1_TAI <- timeDiffs_results$UT1_TAI[[1]]
    UTC_GPS <- timeDiffs_results$UTC_GPS[[1]]
    UT1_GPS <- timeDiffs_results$UT1_GPS[[1]]
    TT_UTC <- timeDiffs_results$TT_UTC[[1]]
    GPS_UTC <- timeDiffs_results$GPS_UTC[[1]]
    PMM <- iauPom00(x_pole, y_pole, iauSp00(2400000.5, MJD_TT))
    # NPB <- iauPnm06a(DJMJD0, TT)
    NPB <- iauPnm06a(2400000.5, MJD_TT)
    # theta <- iauRz(iauGst06(DJMJD0, UT1, DJMJD0, TT, NPB), diag(3))
    theta <- iauRz(iauGst06(2400000.5, MJD_UT1, 2400000.5, MJD_TT, NPB), diag(3))
    E <- PMM %*% theta %*% NPB
    # MJD_TDB <- Mjday_TDB(TT)
    JPL_ephemerides <- JPLephemeridesDE440(MJD_TDB, centralBody, derivativesOrder=2)
    if(autoCentralBodyChange) {
        # realCentralBody <- determineCentralBody(Y[1:3], JPL_ephemerides[-c(1, 2, 12, 13)], JPL_ephemerides[[12]])
        realCentralBody <- determineCentralBody(Y[1:3], JPL_ephemerides[c("positionMercury",
                                                                          "positionVenus",
                                                                          "positionEarth",
                                                                          "positionMars",
                                                                          "positionJupiter",
                                                                          "positionSaturn",
                                                                          "positionUranus",
                                                                          "positionNeptune",
                                                                          "positionPluto")], 
                                                JPL_ephemerides[["positionMoon"]])
    }
    else {
        realCentralBody <- centralBody
    }
    if(centralBody == "Earth") {
        # Acceleration due to Earth, with zonal harmonics
        a <- anelasticEarthAcceleration(MJD_UTC, JPL_ephemerides$positionSun,
                                        JPL_ephemerides$positionMoon, Y[1:3], E,
                                        UT1_UTC, TT_UTC,
                                        x_pole, y_pole, earthSPHDegree, SETcorrections,
                                        OTcorrections)
        # Acceleration due to Moon as point mass only if Moon spherical degree is 1
        if(moonSPHDegree == 1) {
            a <- a + pointMassAcceleration(Y[1:3],JPL_ephemerides$positionMoon,GM_Moon_DE440)
        } else {
            a <- a + anelasticMoonAcceleration(MJD_UTC, Y[1:3] - JPL_ephemerides$positionMoon, 
                                               JPL_ephemerides$positionSun - JPL_ephemerides$positionMoon, 
                                               -JPL_ephemerides$positionMoon,
                                               JPL_ephemerides$lunarLibrationAngles, 
                                               UT1_UTC, TT_UTC, moonSPHDegree, 
                                               SMTcorrections) - 
                GM_Moon_GRGM1200B * JPL_ephemerides$positionMoon/(sqrt(sum(JPL_ephemerides$positionMoon^2)))^3
        }
        # Acceleration due to solar radiation pressure
        a <- a + solarRadiationAcceleration(Y[1:3], JPL_ephemerides, "Earth", solarArea,
                                            satelliteMass, Cr, solarPressureConst, AU)
        # Acceleration due to atmospheric drag
        # Omega according to section III from https://hpiers.obspm.fr/iers/bul/bulb/explanatory.html
        Omega <- omegaEarth - 8.43994809e-10*LOD
        dens <- NRLMSISE00(MJD_UTC, E%*%Y[1:3], UT1_UTC, TT_UTC)["Total"]
        # a <- a + dragAcceleration(dens, Y[1:3], Y[4:6], NPB, satelliteArea, 
        #                           satelliteMass, Cd, Omega)
        a <- a + dragAcceleration(dens, Y[1:3], Y[4:6], E, satelliteArea, 
                                  satelliteMass, Cd, Omega)
        # Relativistic effects
        a <- a + relativity(Y[1:3], Y[4:6])
    } else if(centralBody == "Moon") {
        # Acceleration due to Moon with spherical harmonics
        a_moon <- elasticMoonAcceleration(MJD_UTC, Y[1:3], JPL_ephemerides$positionSun, 
                                            JPL_ephemerides$positionEarth,
                                            JPL_ephemerides$lunarLibrationAngles, 
                                            UT1_UTC, TT_UTC, moonSPHDegree, 
                                            SMTcorrections)
        # Acceleration due to Earth as point mass if spherical harmonics degree is set to 1
        if(earthSPHDegree == 1) {
            a_earth <- pointMassAcceleration(Y[1:3],JPL_ephemerides$positionEarth,GM_Earth_DE440)
        }
        else {
            a_earth <- anelasticEarthAcceleration(MJD_UTC, JPL_ephemerides$positionSun - JPL_ephemerides$positionEarth,
                                                JPL_ephemerides$positionMoon - JPL_ephemerides$positionEarth,
                                                Y[1:3] - JPL_ephemerides$positionEarth, E,
                                                UT1_UTC, TT_UTC, x_pole, y_pole, earthSPHDegree,
                                                #SETcorrections, OTcorrections) - GM_Earth_TT * JPL_ephemerides$positionEarth/(sqrt(sum(JPL_ephemerides$positionEarth^2)))^3
                                                SETcorrections, OTcorrections) - #JPL_ephemerides$accelerationEarth
                anelasticEarthAcceleration(MJD_UTC, JPL_ephemerides$positionSun - JPL_ephemerides$positionEarth,
                                           JPL_ephemerides$positionMoon - JPL_ephemerides$positionEarth,
                                           0 - JPL_ephemerides$positionEarth, E,
                                           UT1_UTC, TT_UTC, x_pole, y_pole, earthSPHDegree,
                                           #SETcorrections, OTcorrections) - GM_Earth_TT * JPL_ephemerides$positionEarth/(sqrt(sum(JPL_ephemerides$positionEarth^2)))^3
                                           SETcorrections, OTcorrections)
            
        }
        a_solarrad <- solarRadiationAcceleration(Y[1:3], JPL_ephemerides, "Moon", solarArea,
                                             satelliteMass, Cr, solarPressureConst, AU)
    } else {
        # Acceleration due to Earth and Moon as point masses
        a_earth <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionEarth, GM_Earth_DE440)
        a_moon <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionMoon, GM_Moon_DE440)
    }
    # Acceleration due to Sun
    a_sun <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionSun,GM_Sun_DE440)
    # Accelerations due to planets
    a_mercury <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionMercury, GM_Mercury_DE440)
    a_venus <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionVenus, GM_Venus_DE440)
    a_mars <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionMars, GM_Mars_DE440)
    a_jupiter <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionJupiter, GM_Jupiter_DE440)
    a_saturn <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionSaturn, GM_Saturn_DE440)
    a_uranus <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionUranus, GM_Uranus_DE440)    
    a_neptune <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionNeptune, GM_Neptune_DE440)
    a_pluto <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionPluto, GM_Pluto_DE440)
    a <- a_earth + a_moon + a_solarrad + a_sun + a_mercury + a_venus + a_mars + a_jupiter +
             a_saturn + a_uranus + a_neptune + a_pluto
    dY <- matrix(c(Y[4:6], a), byrow=TRUE, ncol=3, nrow=2)
    colnames(dY) <- c("X", "Y", "Z")
    rownames(dY) <- c("Velocity", "Acceleration")
    return(list(dY, realCentralBody, a_moon, a_earth, a_solarrad, a_sun, a_mercury, a_venus, a_mars, a_jupiter,
                a_saturn, a_uranus, a_neptune, a_pluto, 
                JPL_ephemerides$positionEarth, JPL_ephemerides$velocityEarth, JPL_ephemerides$accelerationEarth, 
                JPL_ephemerides$positionMoon, JPL_ephemerides$velocityMoon, JPL_ephemerides$accelerationMoon,
                JPL_ephemerides$positionMercury, JPL_ephemerides$velocityMercury, JPL_ephemerides$accelerationMercury,
                JPL_ephemerides$positionVenus, JPL_ephemerides$velocityVenus, JPL_ephemerides$accelerationVenus,
                JPL_ephemerides$positionMars, JPL_ephemerides$velocityMars, JPL_ephemerides$accelerationMars,
                JPL_ephemerides$positionSaturn, JPL_ephemerides$velocitySaturn, JPL_ephemerides$accelerationSaturn,
                JPL_ephemerides$positionJupiter, JPL_ephemerides$velocityJupiter, JPL_ephemerides$accelerationJupiter,
                JPL_ephemerides$positionNeptune, JPL_ephemerides$velocityNeptune, JPL_ephemerides$accelerationNeptune,
                JPL_ephemerides$positionUranus, JPL_ephemerides$velocityUranus, JPL_ephemerides$accelerationUranus,
                JPL_ephemerides$positionPluto, JPL_ephemerides$velocityPluto, JPL_ephemerides$accelerationPluto,
                JPL_ephemerides$positionSun, JPL_ephemerides$velocitySun, JPL_ephemerides$accelerationSun,
                JPL_ephemerides$positionSunSSBarycentric, JPL_ephemerides$velocitySunSSBarycentric, JPL_ephemerides$accelerationSunSSBarycentric,
                JPL_ephemerides$lunarLibrationAngles, JPL_ephemerides$lunarLibrationAnglesDerivatives, JPL_ephemerides$lunarLibrationAnglesSecondDerivatives
                ))
}

