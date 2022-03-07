odeModel <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        state_vector <- state
        results <- accel(t, state_vector, MJD_UTC, solarArea, satelliteMass, 
                              satelliteArea, Cr, Cd, earthSPHDegree, SETcorrections,
                              OTcorrections, moonSPHDegree, centralBody, autoCentralBodyChange)
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
        list(c(dx, dy, dz, d2x, d2y, d2z),
             centralBodiesNum[centralBody])
    })
}

hpop <- function(position, velocity, dateTime, times, satelliteMass, dragArea, 
                 radiationArea, dragCoefficient, radiationCoefficient, 
                 earthSphericalHarmonicsDegree = 130, solidEarthTides=TRUE,
                 oceanTides=TRUE, moonSphericalHarmonicsDegree = 150, 
                 centralBody="Earth", autoCentralBodyChange=TRUE, ...) {
    if(!(centralBody %in% c("SSB", "Mercury", "Venus", "Earth", "Moon",
                           "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto")))
    hasData()
    propagationGlobalVars <- new.env()
    extraArgs <- list(...)
    date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%OS", tz = "UTC")
    year <- date$year + 1900
    month <- date$mon + 1
    day <- date$mday
    hour <- date$hour
    minute <- date$min
    second <- date$sec
    Mjd_UTC <- iauCal2jd(year, month, day, hour, minute, second)$DATE
    MJD_trunc <- trunc(Mjd_UTC)
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
    MJD_TDB <- Mjday_TDB(MJDUTCtoMJDTT(Mjd_UTC))
    JPLephemerides <- JPLephemeridesDE440(MJD_TDB, centralBody=centralBody, derivatives=TRUE)
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
    parameters = list(
        MJD_UTC = Mjd_UTC,
        solarArea = radiationArea,
        satelliteMass = satelliteMass,
        satelliteArea = dragArea,
        Cr = radiationCoefficient,
        Cd = dragCoefficient,
        earthSPHDegree = earthSphericalHarmonicsDegree,
        moonSPHDegree = moonSphericalHarmonicsDegree,
        SETcorrections = solidEarthTides,
        OTcorrections = oceanTides,
        centralBody = realCentralBody,
        autoCentralBodyChange = autoCentralBodyChange,
        globalVarsEnv = environment())
    integration_results <- ode(y=initial_state, times=times, func=odeModel,
                               parms=parameters, method="radau", rtol=1e-13,
                               atol=1e-16, hini=0.01, ...)
    previousStepCentralBodies <- integration_results[, 8]
    oldCentralBody <- previousStepCentralBodies[1]
    if(autoCentralBodyChange & length(unique(previousStepCentralBodies)) > 1) {
        combinedResults <- integration_results
        totalChangePoint <- 0
        while(length(unique(previousStepCentralBodies)) > 1) {
            changePoint <- which(diff(previousStepCentralBodies) != 0)[1] + 1
            newCentralBody <- previousStepCentralBodies[changePoint]
            newTimes <- times[changePoint:length(times)] - times[changePoint]
            newMjd_UTC <- Mjd_UTC + times[changePoint]/86400
            newMJD_TDB <- Mjday_TDB(MJDUTCtoMJDTT(newMjd_UTC))
            JPLephemerides_oldCentralBody <- JPLephemeridesDE440(newMJD_TDB, centralBody=names(centralBodiesNum[oldCentralBody]), derivatives=TRUE)
            newPosition <- integration_results[changePoint, 2:4] -
                JPLephemerides_oldCentralBody[[paste("position", names(centralBodiesNum[newCentralBody]), sep="")]]
            newVelocity <- integration_results[changePoint, 5:7] -
                JPLephemerides_oldCentralBody[[paste("velocity", names(centralBodiesNum[newCentralBody]), sep="")]]
            newInitial_state <- c(newPosition, newVelocity)
            newParameters <- list(
                MJD_UTC = newMjd_UTC,
                solarArea = radiationArea,
                satelliteMass = satelliteMass,
                satelliteArea = dragArea,
                Cr = radiationCoefficient,
                Cd = dragCoefficient,
                earthSPHDegree = earthSphericalHarmonicsDegree,
                moonSPHDegree = moonSphericalHarmonicsDegree,
                SETcorrections = solidEarthTides,
                OTcorrections = oceanTides,
                centralBody = names(which(centralBodiesNum == newCentralBody)),
                autoCentralBodyChange = autoCentralBodyChange,
                globalVarsEnv = environment()
            )
            newIntegration_results <- ode(y=newInitial_state, times=newTimes, func=odeModel,
                                          parms=newParameters, method="radau", rtol=1e-13,
                                          atol=1e-16, hini=0.01, ...)
            totalChangePoint <- totalChangePoint + changePoint 
            combinedResults[totalChangePoint:nrow(combinedResults), ] <- newIntegration_results
            previousStepCentralBodies <- newIntegration_results[, 8]
            oldCentralBody <- previousStepCentralBodies[1]
        }
        combinedResults[, 1] <- times
        integration_results <- combinedResults
    }
    numeric_results <- integration_results[, 1:7]
    central_bodies <- names(centralBodiesNum[integration_results[, 8]])
    output <- cbind(as.data.frame(numeric_results), central_bodies)
    colnames(output) <- c("time", "X", "Y", "Z", "dX", "dY", "dZ", "Central body")
    return(output)
}
