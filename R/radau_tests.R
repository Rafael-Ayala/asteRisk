# # # # library(deSolve)
# # # # library(profvis)
# # # # library(asteRisk)
# # # # library(plotly)
# # # # t = 0
# # # # # Y = matrix(c(-7136143.23821097, -415951.978574370, -505103.379544365, -575.578554821520, 1079.95441989431, 7358.26105845415), nrow=6, ncol = 6)
# # # # Y = c(-7136143.23821097, -415951.978574370, -505103.379544365, -575.578554821520, 1079.95441989431, 7358.26105845415)
# # # # MJD_UTC = 52388.9135185187
# # # # solarArea = 110.5
# # # # satelliteMass = 8000
# # # # satelliteArea = 62.500000000000000
# # # # Cr = 1
# # # # Cd = 2.200000000000000
# # # #
# # # # parameters = list(
# # # #     MJD_UTC = 52388.9135185187,
# # # #     solarArea = 110.5,
# # # #     satelliteMass = 8000,
# # # #     satelliteArea = 62.500000000000000,
# # # #     Cr = 1,
# # # #     Cd = 2.200000000000000
# # # # )
# # # #
# # # # asteRisk:::accel(10000, Y, MJD_UTC, solarArea, satelliteMass, satelliteArea, Cr, Cd)
# # # #
# # # # times <- c(0,1)
# # # #
# # # # test_ode <- profvis(ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters, method="ode45"))
# # # #
# # # # times <- seq(0, by=60, length.out = 1588+1)
# # # # # test_ode2 <- profvis(ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters, method="radau", maxsteps=1, rtol=1e-13, atol=1e-16))
# # # # test_ode2 <- profvis(ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters, method="radau", maxsteps=1, rtol=1e-13, atol=1e-16, hini=0.01))
# # # #
# # # # test_radau <- ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters, method="radau", maxsteps=100, rtol=1e-13, atol=1e-16, hini=0.01)
# # # #
# # # # #test_ode2 <- profvis(ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters))
# # # # asteRisk:::hpop(c(-7136143.23821097, -415951.978574370, -505103.379544365), c(-575.578554821520, 1079.95441989431, 7358.26105845415), 110.5, 62.500000000000000, 8000, 1, 2.2, 2002, 4, 24, 21, 55, 28, times=seq(0, by=60, length.out=10))
# # # #
# # # # asteRisk:::hpop(initialPosition, initialVelocity, molniyaCrossSection,
# # # #                 molniyaCrossSection, molniyaMass, molniyaCr,
# # # #                 molniyaCd, 2006, 6, 25, 0, 33, 42.835, targetTimes)
# # # #
# # # # x <- test_radau[,2]
# # # # y <- test_radau[,3]
# # # # z <- test_radau[,4]
# # # # c <- test_radau[,1]
# # # # data <- data.frame(x,y,z,c)
# # # #
# # # # fig <- plot_ly(data, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines+markers',
# # # #                line = list(width = 6, color = ~c, colorscale = 'Viridis'),
# # # #                marker = list(size = 3.5, color = ~c, colorscale = 'Greens', cmin = -20, cmax = 50))
# # # #
# # # # fig
# # # #
# # #
# # #
# # # library(bvpSolve)
# # #
# # # state_a_r <- c(-14471729.582, -4677558.558, 9369.461)
# # # state_a_v <- c(-3251.691, -3276.008, 4009.228)
# # #
# # # state_b_r <- c(-14471729.582 * 1.5, -4677558.558 * 1.5, 9369.461)
# # #
# # # times <- seq(0, 7200, 60)
# # #
# # # initial_cond <- c(x=state_a_r[1], y=state_a_r[2], z=state_a_r[3], dx=NA, dy=NA, dz=NA)
# # # end_cond <- c(x=state_b_r[1], y=state_b_r[2], z=state_b_r[3], dx=NA, dy=NA, dz=NA)
# # # dateTime <- "2006-06-25 0:33:42.83"
# # # molniyaMass <- 1600
# # # molniyaCrossSection <- 15
# # # molniyaCd <- 2.2
# # # molniyaCr <- 1.2
# # #
# # # date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
# # # year <- date$year + 1900
# # # month <- date$mon + 1
# # # day <- date$mday
# # # hour <- date$hour
# # # minute <- date$min
# # # second <- date$sec
# # # Mjd_UTC <- asteRisk:::MJday(year, month, day, hour, minute, second)
# # #
# # # parameters = list(
# # #     MJD_UTC = Mjd_UTC,
# # #     solarArea = molniyaCrossSection,
# # #     satelliteMass = molniyaMass,
# # #     satelliteArea = molniyaCrossSection,
# # #     Cr = molniyaCr,
# # #     Cd = molniyaCd)
# # #
# # # test_bvp_sol <- bvptwp(yini=initial_cond, x=times, func=asteRisk:::odeModel,
# # #                        yend=end_cond, parms=parameters)
# # #
# # #
# #
# #
# # ## TESTS RINEX
# #
# #
# 
# 
# test=asteRisk:::readGLONASSNavigationRINEXfile("inst/testRINEX.06g")
# testGPS=asteRisk:::readGPSNavigationRINEXfile("inst/brdc0190.18n")
# 
# 
# test_197 <- test$messages[[197]]
# test_197_positions <- c(test_197$positionX, test_197$positionY, test_197$positionZ) * 1000
# test_197_vel <- c(test_197$velocityX, test_197$velocityY, test_197$velocityZ) * 1000
# test_koe <- ECItoKOE(test_197_positions, test_197_vel)
# 
# n <- asteRisk:::semiMajorAxisToMeanMotion(test_koe$semiMajorAxis)
# 
# 
# targetTimes <- seq(0, 720, by=10)
# 
# results_position_matrix <- matrix(nrow=length(targetTimes), ncol=3)
# results_velocity_matrix <- matrix(nrow=length(targetTimes), ncol=3)
# glonass_mass <- 1415
# glonass_area <- 30.85 + 4 # http://navigation-office.esa.int/attachments_12649497_1_Flohrer_ILRS_WS_2015.pdf
# ref_rho <- 0.1570
# glonass_bstar <- (ref_rho*2.2*glonass_area)/(glonass_mass)
# epochDateTime <- "2006-01-07 23:45:00"
# glonass_cd <- 2.2
# glonass_cr <- 1.2
# 
# for(i in 1:length(targetTimes)) {
#     new_result <- sgdp4(n0=n*((2*pi)/(1440)),
#                         e0=test_koe$eccentricity,
#                         i0=test_koe$inclination,
#                         M0=test_koe$meanAnomaly,
#                         omega0=test_koe$argumentPerigee,
#                         OMEGA0=test_koe$longitudeAscendingNode,
#                         Bstar=glonass_bstar,
#                         initialDateTime=epochDateTime, targetTime = targetTimes[i])
#     results_position_matrix[i,] <- new_result[[1]]
#     results_velocity_matrix[i,] <- new_result[[2]]
# }
# last_molniya_propagation <- new_result
# results_position_matrix = cbind(results_position_matrix, targetTimes)
# colnames(results_position_matrix) <- c("x", "y", "z", "time")
# 
# # Let´s verify that the SDP4 algorithm was automatically chosen
# 
# last_molniya_propagation$algorithm
# 
# 
# # Let´s use the HPOP to calculate the position each 2 minutes during a period
# # of 3 hours
# 
# targetTimes <- seq(0, 43200, by=120)
# 
# hpop_results <- hpop(test_197_positions, test_197_vel, epochDateTime,
#                      targetTimes, glonass_mass, glonass_area,
#                      glonass_area, glonass_cr, glonass_cd)
# 
# # Now we can calculate and plot the corresponding geodetic coordinates
# 
# geodetic_matrix_hpop <- matrix(nrow=nrow(hpop_results), ncol=3)
# 
# for(i in 1:nrow(geodetic_matrix_hpop)) {
#     new_dateTime <- as.character(as.POSIXct(epochDateTime, tz="UTC") + targetTimes[i])
#     new_geodetic <- GCRFtoLATLON(hpop_results[i, 2:4], new_dateTime)
#     geodetic_matrix_hpop[i,] <- new_geodetic
# }
# 
# colnames(geodetic_matrix_hpop) <- c("latitude", "longitude", "altitude")
# 
# library(ggmap)
# 
# ggmap(get_map(c(left=-180, right=180, bottom=-80, top=80))) +
#     geom_segment(data=as.data.frame(geodetic_matrix_hpop),
#                  aes(x=longitude, y=latitude,
#                      xend=c(tail(longitude, n=-1), NA),
#                      yend=c(tail(latitude, n=-1), NA)),
#                  na.rm=TRUE) +
#     geom_point(data=as.data.frame(geodetic_matrix_hpop), aes(x=longitude, y=latitude),
#                color="blue", size=0.3, alpha=0.8)
# 
# 
# 
# 
# targetTimes <- seq(0, 86400, by=120)
# 
# epochDateTime = as.character(as.POSIXct(testMes$ephemerisUTCTime, tz="UTC))
# GCRF_posvel=ECEFtoGCRF(testMes$position_ECEF, velocity_ECEF = testMes$velocity_ECEF, dateTime = epochDateTime)
# 
# hpop_results <- hpop(GCRF_posvel$position, GCRF_posvel$velocity, epochDateTime,
#                      targetTimes, 1080, 22+4.7,
#                      22+4.7, 2.2, 1.2)

# # Now we can calculate and plot the corresponding geodetic coordinates
# 
# geodetic_matrix_hpop <- matrix(nrow=nrow(hpop_results), ncol=3)
# 
# for(i in 1:nrow(geodetic_matrix_hpop)) {
#     new_dateTime <- as.character(as.POSIXct(epochDateTime, tz="UTC") + targetTimes[i])
#     new_geodetic <- GCRFtoLATLON(hpop_results[i, 2:4], new_dateTime)
#     geodetic_matrix_hpop[i,] <- new_geodetic
# }
# 
# colnames(geodetic_matrix_hpop) <- c("latitude", "longitude", "altitude")
# 
# library(ggmap)
# 
# ggmap(get_map(c(left=-180, right=180, bottom=-80, top=80))) +
#     geom_segment(data=as.data.frame(geodetic_matrix_hpop),
#                  aes(x=longitude, y=latitude,
#                      xend=c(tail(longitude, n=-1), NA),
#                      yend=c(tail(latitude, n=-1), NA)),
#                  na.rm=TRUE) +
#     geom_point(data=as.data.frame(geodetic_matrix_hpop), aes(x=longitude, y=latitude),
#                color="blue", size=0.3, alpha=0.8)
# 
# 
# GCRF_posvel=KOEtoECI(testMes$semiMajorAxis, testMes$eccentricity, testMes$inclination, testMes$meanAnomaly, testMes$perigeeArgument, testMes$ascension)
# 
# hpop_results <- hpop(GCRF_posvel$position, GCRF_posvel$velocity, epochDateTime,
#                      targetTimes, 1080, 22+4.7,
#                      22+4.7, 2.2, 1.2)
# # 
# test_TLEs <- readTLE(paste0(path.package("asteRisk"), "/testTLE.txt"))
# test_TLE1=test_TLEs[[1]]
# 
# sgdp4(n0=test_TLE1$meanMotion*((2*pi)/(1440)),
#       e0=test_TLE1$eccentricity,
#       i0=test_TLE1$inclination*pi/180,
#       M0=test_TLE1$meanAnomaly*pi/180,
#       omega0=test_TLE1$perigeeArgument*pi/180,
#       OMEGA0=test_TLE1$ascension*pi/180,
#       Bstar=test_TLE1$Bstar,
#       targetTime = 0)
# 
# results_position_matrix <- matrix(nrow=length(target_times_1), ncol=3)
# results_velocity_matrix <- matrix(nrow=length(target_times_1), ncol=3)
# 
# for(i in 1:length(target_times_1)) {
#     new_result <- sgdp4(n0=test_TLE1$meanMotion*((2*pi)/(1440)),
#                         e0=test_TLE1$eccentricity,
#                         i0=test_TLE1$inclination*pi/180,
#                         M0=test_TLE1$meanAnomaly*pi/180,
#                         omega0=test_TLE1$perigeeArgument*pi/180,
#                         OMEGA0=test_TLE1$ascension*pi/180,
#                         Bstar=test_TLE1$Bstar,
#                         initialDateTime=test_TLE1$dateTime, targetTime = target_times_1[i])
#     results_position_matrix[i,] <- new_result[[1]]
#     results_velocity_matrix[i,] <- new_result[[2]]
# }
# 
# test_TLE11=test_TLEs[[11]]
# target_times_11 = seq(0, 2880, by=120)
# 
# sgdp4(n0=test_TLE11$meanMotion*((2*pi)/(1440)),
#       e0=test_TLE11$eccentricity,
#       i0=test_TLE11$inclination*pi/180,
#       M0=test_TLE11$meanAnomaly*pi/180,
#       omega0=test_TLE11$perigeeArgument*pi/180,
#       OMEGA0=test_TLE11$ascension*pi/180,
#       Bstar=test_TLE11$Bstar,
#       targetTime = 2880, initialDateTime = test_TLE11$dateTime)
# 
# results_position_matrix <- matrix(nrow=length(target_times_11), ncol=3)
# results_velocity_matrix <- matrix(nrow=length(target_times_11), ncol=3)
# 
# for(i in 1:length(target_times_11)) {
#     new_result <- sgdp4(n0=test_TLE11$meanMotion*((2*pi)/(1440)),
#                         e0=test_TLE11$eccentricity,
#                         i0=test_TLE11$inclination*pi/180,
#                         M0=test_TLE11$meanAnomaly*pi/180,
#                         omega0=test_TLE11$perigeeArgument*pi/180,
#                         OMEGA0=test_TLE11$ascension*pi/180,
#                         Bstar=test_TLE11$Bstar,
#                         initialDateTime=test_TLE11$dateTime, targetTime = target_times_11[i])
#     results_position_matrix[i,] <- new_result[[1]]
#     results_velocity_matrix[i,] <- new_result[[2]]
# }
# 


# # testline1 = '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991'
# # testline2 = '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'
# # 
# # testTLEpy = parseTLElines(c(testline1, testline2))
# targetTimePy = "2019-12-9 20:42:9.07"
# 
# sgp4(n0=testTLEpy$meanMotion*((2*pi)/(1440)),
#      e0=testTLEpy$eccentricity,
#      i0=testTLEpy$inclination*pi/180,
#      M0=testTLEpy$meanAnomaly*pi/180,
#      omega0=testTLEpy$perigeeArgument*pi/180,
#      OMEGA0=testTLEpy$ascension*pi/180,
#      Bstar = testTLEpy$Bstar,
#      initialDateTime = testTLEpy$dateTime,
#      targetTime = targetTimePy)


# Checks for planet lab cubesats propagation
# Planetlabs uses 3U cubesats (10cmx10cmx30cm, 3.99 kg)
# Planetlabs epochs are given in s since J2000 TT epoch, which is  11:58:55.816 UTC 1st of January 2000
# Using hardware ID 0505 ephemeris at 2021-04-10, target 2021-04-25
# # 
# #
# getLatestSpaceData()
# 
# planet_lab_mass <- 3.99
# planet_lab_area <- 0.5*(10+10+30)*0.0001
# planet_lab_cd <- 2.2
# planet_lab_cr <- 1.2
# planet_lab_initial_position <- c(-6.791120798515737e+06,  1.298567476543681e+06, -4.919967206741248e+05)
# planet_lab_initial_velocity <- c(-1.083847203010712e+03, -3.073534729296989e+03,  6.844313262593339e+03)
# planet_lab_initial_dateTime <- as.POSIXct(671365399.184000, origin="2000-01-01 11:58:55.816", tz="UTC")
# planet_lab_target_dateTime <- as.POSIXct(672661404.184000, origin="2000-01-01 11:58:55.816", tz="UTC")
# target_time <- 672661404.184000 - 671365399.184000
# 
# propagation_planet_lab <- hpop(position = planet_lab_initial_position,
#                                velocity = planet_lab_initial_velocity,
#                                dateTime = planet_lab_initial_dateTime,
#                                satelliteMass = planet_lab_mass,
#                                dragArea = planet_lab_area,
#                                radiationArea = planet_lab_area,
#                                dragCoefficient = planet_lab_cd,
#                                radiationCoefficient = planet_lab_cr,
#                                times = seq(0, target_time, by=10800),
#                                maxsteps = 100000)
# 
# # saveRDS(propagation_planet_lab, "propagation_planet_lab.rds")
# # 
# planet_lab_initial_dateTime2 <- as.POSIXct(671365399.184000 + 1296000, origin="2000-01-01 11:58:55.816", tz="UTC")
# planet_lab_initial_position2 <- propagation_planet_lab[121, 2:4]
# planet_lab_initial_velocity2 <- propagation_planet_lab[121, 5:7]
# 
# propagation_planet_lab_2 <- hpop(position = planet_lab_initial_position2,
#                                  velocity = planet_lab_initial_velocity2,
#                                  dateTime = planet_lab_initial_dateTime2,
#                                  satelliteMass = planet_lab_mass,
#                                  dragArea = planet_lab_area,
#                                  radiationArea = planet_lab_area,
#                                  dragCoefficient = planet_lab_cd,
#                                  radiationCoefficient = planet_lab_cr,
#                                  times = seq(0, 5, by=1),
#                                  maxsteps = 100000)
# 
# # Compare with TLE
# 
# planet_lab_TLE1 <- readTLE(paste0(path.package("asteRisk"), "/planet_mc_20210410.tle"))
# 
# planet_lab_0505_TLE1 <- planet_lab_TLE1[[1]]
# 
# sgp4_epoch_position <- sgp4(n0=planet_lab_0505_TLE1$meanMotion*((2*pi)/(1440)),
#      e0=planet_lab_0505_TLE1$eccentricity,
#      i0=planet_lab_0505_TLE1$inclination*pi/180,
#      M0=planet_lab_0505_TLE1$meanAnomaly*pi/180,
#      omega0=planet_lab_0505_TLE1$perigeeArgument*pi/180,
#      OMEGA0=planet_lab_0505_TLE1$ascension*pi/180,
#      Bstar = planet_lab_0505_TLE1$Bstar,
#      initialDateTime = planet_lab_0505_TLE1$dateTime,
#      targetTime = 0)
# 
# sgp4_GCRF_epoch_position <- TEMEtoGCRF(sgp4_epoch_position$position*1000, sgp4_epoch_position$velocity*1000, planet_lab_0505_TLE1$dateTime)
# 
# sgp4_final_position <- sgp4(n0=planet_lab_0505_TLE1$meanMotion*((2*pi)/(1440)),
#                             e0=planet_lab_0505_TLE1$eccentricity,
#                             i0=planet_lab_0505_TLE1$inclination*pi/180,
#                             M0=planet_lab_0505_TLE1$meanAnomaly*pi/180,
#                             omega0=planet_lab_0505_TLE1$perigeeArgument*pi/180,
#                             OMEGA0=planet_lab_0505_TLE1$ascension*pi/180,
#                             Bstar = planet_lab_0505_TLE1$Bstar,
#                             initialDateTime = planet_lab_0505_TLE1$dateTime,
#                             targetTime = (672661404.184000 - 671365399.184000)/60)
# 
# sgp4_GCRF_final_position <- TEMEtoGCRF(sgp4_final_position$position*1000, sgp4_final_position$velocity*1000, planet_lab_0505_TLE1$dateTime)
# 
# 
# 
# # Testing proapgation of GPS RINEX for 1 day interval
# 
# RINEXfile_2021010 <- readGPSNavigationRINEX(paste0(path.package("asteRisk"), "/brdc0100.21n"))
# RINEXfile_2021011 <- readGPSNavigationRINEX(paste0(path.package("asteRisk"), "/brdc0110.21n"))
# 
# RINEX_2021010_msg1 <- RINEXfile_2021010$messages[[1]]
# RINEX_2021011_msg1 <- RINEXfile_2021011$messages[[1]]
# 
# GPS_initial_dateTime <- "2021-01-09 23:59:41.999217172"
# GPS_initial_GCRF <- ECEFtoGCRF(RINEX_2021010_msg1$position_ECEF, RINEX_2021010_msg1$velocity_ECEF, GPS_initial_dateTime)
# GPS_target_dateTime <- "2021-01-10 23:59:41.999217712"
# GPS_target_GCRF <-  ECEFtoGCRF(RINEX_2021011_msg1$position_ECEF, RINEX_2021011_msg1$velocity_ECEF, GPS_target_dateTime)
# 
# GPS_initial_position <- GPS_initial_GCRF$position
# GPS_initial_velocity <- GPS_initial_GCRF$velocity
# gps_mass <- 1633 # satellite of block IIF
# gps_area <- 0.5*(5.4 + 7.01 + 5.72 + 22.25)*0.85 # from https://mediatum.ub.tum.de/doc/1188612/719708.pdf, the 0.85 is shadowing factor
# gps_cd <- 2.2
# gps_cr <- 1.2
# 
# propagation_GPS <- hpop(position = GPS_initial_position,
#                         velocity = GPS_initial_velocity,
#                         dateTime = GPS_initial_dateTime,
#                         satelliteMass = gps_mass,
#                         dragArea = gps_area,
#                         radiationArea = gps_area,
#                         dragCoefficient = gps_cd,
#                         radiationCoefficient = gps_cr,
#                         times = seq(0, 86400, by=3600))
# 
# test_GCRF_propagation_position <- propagation_GPS[25, 2:4]
# test_GCRF_propagation_velocity <- propagation_GPS[25, 5:7]
# 
# sqrt(sum((GPS_target_GCRF$position - test_GCRF_propagation_position)^2))



##### CHECKS FOR CONVERSION BETWEEN FRAMES
# TEME to GCRF
# 
# r_teme <- c(5094.18016210, 6127.64465950, 6380.34453270)
# dateTime <- "2004-04-06 7:51:28.386009"
# TEMEtoGCRF(r_teme, c(0,0,0), dateTime)
