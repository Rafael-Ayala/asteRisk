odeModel <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        state_vector <- state
        acceleration <- accel(t, state_vector, MJD_UTC, solarArea, satelliteMass, satelliteArea, Cr, Cd)
        dx <- acceleration[1, 1]
        dy <- acceleration[1, 2]
        dz <- acceleration[1, 3]
        d2x <- acceleration[2, 1]
        d2y <- acceleration[2, 2]
        d2z <- acceleration[2, 3]
        list(c(dx, dy, dz, d2x, d2y, d2z))
    })
}

hpop <- function(position, velocity, dateTime, times, satelliteMass, dragArea, 
                 radiationArea, dragCoefficient, radiationCoefficient, ...) {
    hasData()
    extraArgs <- list(...)
    date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
    year <- date$year + 1900
    month <- date$mon + 1
    day <- date$mday
    hour <- date$hour
    minute <- date$min
    second <- date$sec
    initial_state <- c(position, velocity)
    Mjd_UTC <- iauCal2jd(year, month, day, hour, minute, second)$DATE
    parameters = list(
        MJD_UTC = Mjd_UTC,
        solarArea = radiationArea,
        satelliteMass = satelliteMass,
        satelliteArea = dragArea,
        Cr = radiationCoefficient,
        Cd = dragCoefficient)
    integration_results <- ode(y=initial_state, times=times, func=odeModel,
                               parms=parameters, method="radau", rtol=1e-13,
                               atol=1e-16, hini=0.01, ...)
    colnames(integration_results) <- c("time", "X", "Y", "Z", "dX", "dY", "dZ")
    return(integration_results)
}
