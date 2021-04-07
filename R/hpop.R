hpop <- function(position, velocity, radiationArea, dragArea, satMass,
                 radiationCoefficient, dragCoefficient) {
    # load epoch state
    year <- 2002
    mon <- 4
    day <- 24
    hour <- 21
    min <- 55
    sec <- 28
    Y0 <- c( # state vector
        +7144843.808,
        +0217687.110,
        -0506463.296,
        +0562.650611,
        -1616.516697,
        +7358.157263)
    area_solar <- 110.5
    area_drag <- 62.5
    mass <- 8000
    Cr <- 1.0
    Cd <- 2.2
    Mjd_UTC <- MJday(year, mon, day, hour, min, sec)
    Y0 <- ECEFtoECI(Mjd_UTC, Y0)
    Mjd0 <- Mjd_UTC
    step <- 60 # integration step in s
    n_step <- 1588 # number of integration steps to perform
    # propagation phase starts here
    Eph <- Ephemeris(Y0, n_step, step)

}