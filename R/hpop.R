odeModel <- function(t, state, parameters) {
    if(t %% 10000 == 0){
        print(t)
    }
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

hpop <- function(position, velocity, radiationArea, dragArea, satMass,
                 radiationCoefficient, dragCoefficient, year, month, day,
                 hour, minute, second, times, ...) {
    initial_state <- c(position, velocity)
    Mjd_UTC <- MJday(year, month, day, hour, minute, second)
    parameters = list(
        MJD_UTC = Mjd_UTC,
        solarArea = radiationArea,
        satelliteMass = satMass,
        satelliteArea = dragArea,
        Cr = radiationCoefficient,
        Cd = dragCoefficient)
    integration_results <- ode(y=initial_state, times=times, func=asteRisk:::odeModel, 
                      parms=parameters, method="radau", maxsteps=100, rtol=1e-13, 
                      atol=1e-16, hini=0.01)
    return(integration_results)
}
# 
# hpop <- function(position, velocity, radiationArea, dragArea, satMass,
#                  radiationCoefficient, dragCoefficient, year, month, day,
#                  hour, minute, second, times, ...) {
#     Mjd_UTC <- MJday(year, month, day, hour, minute, second)
#     year <- 2002
#     mon <- 4
#     day <- 24
#     hour <- 21
#     min <- 55
#     sec <- 28
#     Y0 <- c( # state vector
#         +7144843.808,
#         +0217687.110,
#         -0506463.296,
#         +0562.650611,
#         -1616.516697,
#         +7358.157263)
#     area_solar <- 110.5
#     area_drag <- 62.5
#     mass <- 8000
#     Cr <- 1.0
#     Cd <- 2.2
#     Mjd_UTC <- MJday(year, mon, day, hour, min, sec)
#     Y0 <- ECEFtoECI(Mjd_UTC, Y0)
#     Mjd0 <- Mjd_UTC
#     step <- 60 # integration step in s
#     n_step <- 1588 # number of integration steps to perform
#     # propagation phase starts here
#     Eph <- Ephemeris(Y0, n_step, step)
# 
# }