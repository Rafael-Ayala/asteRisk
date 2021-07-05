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
