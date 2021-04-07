elasticEarthAcceleration <- function(Mjd_UTC, r_sun, r_moon, r, E, UT1_UTC,
                                     TT_UTC, x_pole, y_pole) {
    C <- Cnm
    S <- Snm
    r_moon <- E %*% r_moon
    moon_polar <- CartesianToPolar(r_moon)
    r_sun <- E %*% r_sun
    sun_polar <- CartesianToPolar(r_sun)
    Mjd_TT <- Mjd_UTC + TT_UTC/86400
    Mjd_UT1 <- Mjd_UTC + UT1_UTC/86400
    Time <- (Mjd_TT-MJD_J2000)/36525
    # Start corrections for solid Earth tides
    legendre_moonTheta <- legendre(4, 4, moon_polar[2])
    legendre_sunTheta <- legendre(4, 4, sun_polar[2])
    # step 1 of corrections for solid Earth tides
    dCnm20 <- (0.29525/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)
                            *legendre_moonTheta$normLegendreValues[3,1] + 
                               (GM_Sun/gm)*((r_ref/sun_polar[3])^3)
                            *legendre_sunTheta$normLegendreValues[3,1] )
    dCnm21 <- (0.29470/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)
                            *legendre_moonTheta$normLegendreValues[3,2]*cos(moon_polar[1]) + 
                                (GM_Sun/gm)*((r_ref/sun_polar[3])^3)
                            *legendre_sunTheta$normLegendreValues[3,2]*cos(sun_polar[1]) )
    dSnm21 <- (0.29470/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)
                            *legendre_moonTheta$normLegendreValues[3,2]*(sin(moon_polar[1]))+ 
                                (GM_Sun/gm)*((r_ref/sun_polar[3])^3)
                            *legendre_sunTheta$normLegendreValues[3,2]*(sin(sun_polar[1])) )
    dCnm22 <- (0.29801/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)
                            *legendre_moonTheta$normLegendreValues[3,3]*cos(2*moon_polar[1])+ 
                                (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*
                                legendre_sunTheta$normLegendreValues[3,3]*cos(2*sun_polar[1]) )
    dSnm22 <- (0.29801/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*
                                legendre_moonTheta$normLegendreValues[3,3]*(sin(2*moon_polar[1]))+
                                (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*
                                legendre_sunTheta$normLegendreValues[3,3]*(sin(2*sun_polar[1])) )
    dCnm30 <- (0.093/7)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*
                              legendre_moonTheta$normLegendreValues[4,1]+ (GM_Sun/gm)*
                              ((r_ref/sun_polar[3])^4)*legendre_sunTheta$normLegendreValues[4,1] )
    dCnm31 <- (0.093/7)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*
                              legendre_moonTheta$normLegendreValues[4,2]*cos(moon_polar[1])+ 
                              (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*
                              legendre_sunTheta$normLegendreValues[4,2]*cos(sun_polar[1]) )
    dSnm31 <- (0.093/7)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*
                              legendre_moonTheta$normLegendreValues[4,2]*sin(moon_polar[1])+ 
                              (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*
                              legendre_sunTheta$normLegendreValues[4,2]*sin(sun_polar[1]) )
    dCnm32 <- (0.093/7)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*
                              legendre_moonTheta$normLegendreValues[4,3]*cos(2*moon_polar[1])+ 
                              (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*
                              legendre_sunTheta$normLegendreValues[4,3]*cos(2*sun_polar[1]) )
    dSnm32 <- (0.093/7)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*
                              legendre_moonTheta$normLegendreValues[4,3]*sin(2*moon_polar[1])+ 
                              (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*
                              legendre_sunTheta$normLegendreValues[4,3]*sin(2*sun_polar[1]) )
    dCnm33 <- (0.094/7)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*
                              legendre_moonTheta$normLegendreValues[4,4]*cos(3*moon_polar[1])+ 
                              (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*
                              legendre_sunTheta$normLegendreValues[4,4]*cos(3*sun_polar[1]) )
    dSnm33 <- (0.094/7)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*
                              legendre_moonTheta$normLegendreValues[4,4]*sin(3*moon_polar[1])+ 
                              (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*
                              legendre_sunTheta$normLegendreValues[4,4]*sin(3*sun_polar[1]) )
    dCnm40 <- (-0.00087/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*
                                 legendre_moonTheta$normLegendreValues[5,1]+ (GM_Sun/gm)*
                                 ((r_ref/sun_polar[3])^3)*legendre_sunTheta$normLegendreValues[5,1])
    dCnm41 <- (-0.00079/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*
                                 legendre_moonTheta$normLegendreValues[5,2]*cos(moon_polar[1])+
                                 (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*
                                 legendre_sunTheta$normLegendreValues[5,2]*cos(sun_polar[1]) )
    dSnm41 <- (-0.00079/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*
                                 legendre_moonTheta$normLegendreValues[5,2]*sin(moon_polar[1])+
                                 (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*
                                 legendre_sunTheta$normLegendreValues[5,2]*sin(sun_polar[1]) )
    dCnm42 <- (-0.00057/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*
                                 legendre_moonTheta$normLegendreValues[5,3]*cos(2*moon_polar[1])+
                                 (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*
                                 legendre_sunTheta$normLegendreValues[5,3]*cos(2*sun_polar[1]) )
    dSnm42 <- (-0.00057/5)*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*
                                 legendre_moonTheta$normLegendreValues[5,3]*sin(2*moon_polar[1])+
                                 (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*
                                 legendre_sunTheta$normLegendreValues[5,3]*sin(2*sun_polar[1]) )
    # step 2 of corrections for solid Earth tides
    invjday_results <- invjday(Mjd_UTC+2400000.5)
    iauCal2jd_results <- iauCal2jd(invjday_results$year,
                                   invjday_results$month,
                                   invjday_results$day)
    TIME <- (60*(60*invjday_results$hour+invjday_results$min)+invjday_results$sec)/86400
    UTC <- iauCal2jd_results$DATE+TIME
    TT <- UTC+TT_UTC/86400
    TUT <- TIME+UT1_UTC/86400
    UT1 <- iauCal2jd_results$DATE+TUT
    theta_g <- iauGmst06(iauCal2jd_results$DJMJD0, UT1, iauCal2jd_results$DJMJD0, TT)
    dC21 <- 0
    dS21 <- 0
    for (i in 1:48) {
        dC21 <- dC21 + 1e-12*solidEarthTides_dC21dS21[i, 6]*sin(theta_g+pi)
        dS21 <- dS21 + 1e-12*solidEarthTides_dC21dS21[i, 6]*cos(theta_g+pi)
    }
    dCnm21 <- dCnm21 + dC21
    dSnm21 <- dSnm21 + dS21
    dC22 <- 0
    dS22 <- 0
    for (i in 1:2) {
        dC22 <- dC22 + 1e-12*solidEarthTides_dC22dS22[i, 6]*sin(theta_g+pi)
        dS22 <- dS22 + 1e-12*solidEarthTides_dC22dS22[i, 6]*cos(theta_g+pi)
    }
    dCnm22 <- dCnm22 + dC22
    dSnm22 <- dSnm22 + dS22
    # Effect of permanent tide (elastic Earth)
    dC20 <- 4.4228e-8*(-0.31460)*0.29525
    dCnm20 <- dCnm20 - dC20
    # Effect of solid Earth pole tide (elastic Earth)
    dC21 <- -1.290e-9*x_pole
    dS21 <- 1.290e-9*y_pole
    dCnm21 <- dCnm21 + dC21
    dSnm21 <- dSnm21 + dS21
    C[3,1] <- C[3,1] + dCnm20
    C[3,2] <- C[3,2] + dCnm21
    C[3,3] <- C[3,3] + dCnm22
    S[3,2] <- S[3,2] + dSnm21
    S[3,3] <- S[3,3] + dSnm22
    C[4,1] <- C[4,1] + dCnm30
    C[4,2] <- C[4,2] + dCnm31
    C[4,3] <- C[4,3] + dCnm32
    C[4,4] <- C[4,4] + dCnm33
    S[4,2] <- S[4,2] + dSnm31
    S[4,3] <- S[4,3] + dSnm32
    S[4,4] <- S[4,4] + dSnm33
    C[5,1] <- C[5,1] + dCnm40
    C[5,2] <- C[5,2] + dCnm41
    C[5,3] <- C[5,3] + dCnm42
    S[5,2] <- S[5,2] + dSnm41
    S[5,3] <- S[5,3] + dSnm42    
    # End corrections for solid Earth tides
    # Start block for corrections due to ocean tides
    legendre_moonTheta <- legendre(6, 6, moon_polar[2])
    legendre_sunTheta <- legendre(6, 6, sun_polar[2])
    
    dCnm20 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*legendre_moonTheta$normLegendreValues[3,1] + (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*legendre_sunTheta$normLegendreValues[3,1] )
    dCnm21 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*legendre_moonTheta$normLegendreValues[3,2]*cos(moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*legendre_sunTheta$normLegendreValues[3,2]*cos(sun_polar[1]) )
    dSnm21 <- -0.3075/5*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*legendre_moonTheta$normLegendreValues[3,2]*sin(moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*legendre_sunTheta$normLegendreValues[3,2]*sin(sun_polar[1]) )
    dCnm22 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*legendre_moonTheta$normLegendreValues[3,3]*cos(2*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*legendre_sunTheta$normLegendreValues[3,3]*cos(2*sun_polar[1]) )
    dSnm22 <- -0.3075/5*( (GM_Moon/gm)*((r_ref/moon_polar[3])^3)*legendre_moonTheta$normLegendreValues[3,3]*sin(2*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^3)*legendre_sunTheta$normLegendreValues[3,3]*sin(2*sun_polar[1]) )
    dCnm30 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*legendre_moonTheta$normLegendreValues[4,1] + (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*legendre_sunTheta$normLegendreValues[4,1] )
    dCnm31 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*legendre_moonTheta$normLegendreValues[4,2]*cos(moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*legendre_sunTheta$normLegendreValues[4,2]*cos(sun_polar[1]) )
    dSnm31 <- -0.195/7*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*legendre_moonTheta$normLegendreValues[4,2]*sin(moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*legendre_sunTheta$normLegendreValues[4,2]*sin(sun_polar[1]) )
    dCnm32 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*legendre_moonTheta$normLegendreValues[4,3]*cos(2*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*legendre_sunTheta$normLegendreValues[4,3]*cos(2*sun_polar[1]) )
    dSnm32 <- -0.195/7*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*legendre_moonTheta$normLegendreValues[4,3]*sin(2*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*legendre_sunTheta$normLegendreValues[4,3]*sin(2*sun_polar[1]) )
    dCnm33 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*legendre_moonTheta$normLegendreValues[4,4]*cos(3*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*legendre_sunTheta$normLegendreValues[4,4]*cos(3*sun_polar[1]) )
    dSnm33 <- -0.195/7*( (GM_Moon/gm)*((r_ref/moon_polar[3])^4)*legendre_moonTheta$normLegendreValues[4,4]*sin(3*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^4)*legendre_sunTheta$normLegendreValues[4,4]*sin(3*sun_polar[1]) )
    dCnm40 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,1] + (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,1] )
    dCnm41 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,2]*cos(moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,2]*cos(sun_polar[1]) )
    dSnm41 <- -0.132/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,2]*sin(moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,2]*sin(sun_polar[1]) )
    dCnm42 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,3]*cos(2*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,3]*cos(2*sun_polar[1]) )
    dSnm42 <- -0.132/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,3]*sin(2*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,3]*sin(2*sun_polar[1]) )
    dCnm43 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,4]*cos(3*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,4]*cos(3*sun_polar[1]) )
    dSnm43 <- -0.132/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,4]*sin(3*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,4]*sin(3*sun_polar[1]) )
    dCnm44 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,5]*cos(4*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,5]*cos(4*sun_polar[1]) )
    dSnm44 <- -0.132/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^5)*legendre_moonTheta$normLegendreValues[5,5]*sin(4*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^5)*legendre_sunTheta$normLegendreValues[5,5]*sin(4*sun_polar[1]) )
    dCnm50 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,1]+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,1] )
    dCnm51 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,2]*cos(moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,2]*cos(sun_polar[1]) )
    dSnm51 <- -0.1032/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,2]*sin(moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,2]*sin(sun_polar[1]) )
    dCnm52 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,3]*cos(2*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,3]*cos(2*sun_polar[1]) )
    dSnm52 <- -0.1032/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,3]*sin(2*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,3]*sin(2*sun_polar[1]) )
    dCnm53 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,4]*cos(3*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,4]*cos(3*sun_polar[1]) )
    dSnm53 <- -0.1032/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,4]*sin(3*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,4]*sin(3*sun_polar[1]) )
    dCnm54 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,5]*cos(4*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,5]*cos(4*sun_polar[1]) )
    dSnm54 <- -0.1032/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,5]*sin(4*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,5]*sin(4*sun_polar[1]) )
    dCnm55 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,6]*cos(5*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,6]*cos(5*sun_polar[1]) )
    dSnm55 <- -0.1032/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^6)*legendre_moonTheta$normLegendreValues[6,6]*sin(5*moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^6)*legendre_sunTheta$normLegendreValues[6,6]*sin(5*sun_polar[1]) )
    dCnm60 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,1]+ (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,1] )
    dCnm61 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,2]*cos(moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,2]*cos(sun_polar[1]) )
    dSnm61 <- -0.0892/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,2]*sin(moon_polar[1])+ (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,2]*sin(sun_polar[1]) )
    dCnm62 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,3]*cos(2*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,3]*cos(2*sun_polar[1]) )
    dSnm62 <- -0.0892/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,3]*sin(2*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,3]*sin(2*sun_polar[1]) )
    dCnm63 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,4]*cos(3*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,4]*cos(3*sun_polar[1]) )
    dSnm63 <- -0.0892/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,4]*sin(3*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,4]*sin(3*sun_polar[1]) )
    dCnm64 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,5]*cos(4*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,5]*cos(4*sun_polar[1]) )
    dSnm64 <- -0.0892/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,5]*sin(4*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,5]*sin(4*sun_polar[1]) )
    dCnm65 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,6]*cos(5*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,6]*cos(5*sun_polar[1]) )
    dSnm65 <- -0.0892/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,6]*sin(5*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,6]*sin(5*sun_polar[1]) )
    dCnm66 <- 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,7]*cos(6*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,7]*cos(6*sun_polar[1]) )
    dSnm66 <- -0.0892/9*( (GM_Moon/gm)*((r_ref/moon_polar[3])^7)*legendre_moonTheta$normLegendreValues[7,7]*sin(6*moon_polar[1]) + (GM_Sun/gm)*((r_ref/sun_polar[3])^7)*legendre_sunTheta$normLegendreValues[7,7]*sin(6*sun_polar[1]) )
    C[3,1] <- C[3,1] + dCnm20
    C[3,2] <- C[3,2] + dCnm21
    C[3,3] <- C[3,3] + dCnm22
    S[3,2] <- S[3,2] + dSnm21
    S[3,3] <- S[3,3] + dSnm22
    C[4,1] <- C[4,1] + dCnm30
    C[4,2] <- C[4,2] + dCnm31
    C[4,3] <- C[4,3] + dCnm32
    C[4,4] <- C[4,4] + dCnm33
    S[4,2] <- S[4,2] + dSnm31
    S[4,3] <- S[4,3] + dSnm32
    S[4,4] <- S[4,4] + dSnm33
    C[5,1] <- C[5,1] + dCnm40
    C[5,2] <- C[5,2] + dCnm41
    C[5,3] <- C[5,3] + dCnm42
    C[5,4] <- C[5,4] + dCnm43
    C[5,5] <- C[5,5] + dCnm44
    S[5,2] <- S[5,2] + dSnm41
    S[5,3] <- S[5,3] + dSnm42
    S[5,4] <- S[5,4] + dSnm43
    S[5,5] <- S[5,5] + dSnm44
    C[6,1] <- C[6,1] + dCnm50
    C[6,2] <- C[6,2] + dCnm51
    C[6,3] <- C[6,3] + dCnm52
    C[6,4] <- C[6,4] + dCnm53
    C[6,5] <- C[6,5] + dCnm54
    C[6,6] <- C[6,6] + dCnm55
    S[6,2] <- S[6,2] + dSnm51
    S[6,3] <- S[6,3] + dSnm52
    S[6,4] <- S[6,4] + dSnm53
    S[6,5] <- S[6,5] + dSnm54
    S[6,6] <- S[6,6] + dSnm55
    C[7,1] <- C[7,1] + dCnm60
    C[7,2] <- C[7,2] + dCnm61
    C[7,3] <- C[7,3] + dCnm62
    C[7,4] <- C[7,4] + dCnm63
    C[7,5] <- C[7,5] + dCnm64
    C[7,6] <- C[7,6] + dCnm65
    C[7,7] <- C[7,7] + dCnm66
    S[7,2] <- S[7,2] + dSnm61
    S[7,3] <- S[7,3] + dSnm62
    S[7,4] <- S[7,4] + dSnm63
    S[7,5] <- S[7,5] + dSnm64
    S[7,6] <- S[7,6] + dSnm65
    S[7,7] <- S[7,7] + dSnm66    
    # End of ocean tide correction
    r_bf <- E %*% r
    d <- sqrt(sum(r_bf^2))
    latgc <- asin(r_bf[3]/d)
    lon <- atan2(r_bf[2], r_bf[1])
    legendre_latgc <- legendre(n, m, latgc)
    dUdr <- 0
    dUdlatgc <- 0
    dUdlon <- 0
    q1 <- 0 
    q2 <- 0
    q3 <- 0
    for(n in 0:n) {
        b1 <- (-gm/d^2)*(r_ref/d)^n*(n+1)
        b2 <-  (gm/d)*(r_ref/d)^n
        b3 <-  (gm/d)*(r_ref/d)^n
        for(m in 0:m) {
            q1 <- q1 + legendre_latgc$normLegendreValues[n+1,m+1]*(C[n+1,m+1]*cos(m*lon)+S[n+1,m+1]*sin(m*lon))
            q2 <- q2 + legendre_latgc$normLegendreDerivativeValues[n+1,m+1]*
                (C[n+1,m+1]*cos(m*lon)+S[n+1,m+1]*sin(m*lon))
            q3 <- q3 + m*legendre_latgc$normLegendreValues[n+1,m+1]*(S[n+1,m+1]*cos(m*lon)-C[n+1,m+1]*sin(m*lon));
        }
        dUdr <- dUdr + q1*b1
        dUdlatgc <- dUdlatgc + q2*b2
        dUdlon <- dUdlon + q3*b3
        q3 <- 0
        q2 <- 0
        q1 <- 0
    }
    r2xy <- r_bf[1]^2+r_bf[2]^2
    ax <- (1/d*dUdr-r_bf[3]/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf[1]-(1/r2xy*dUdlon)*r_bf[2]
    ay <- (1/d*dUdr-r_bf[3]/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf[2]+(1/r2xy*dUdlon)*r_bf[1]
    az <- 1/d*dUdr*r_bf[3]+sqrt(r2xy)/d^2*dUdlatgc
    a_bf <- c(ax, ay, az)
    a <- t(E) %*% a_bf
    return(a)
}

# asteRisk:::elasticEarthAcceleration(
#     52388.9135185187,
#     c(          123994411196.822,
#                 78220316016.5757,
#                 33912542828.469),
#     c(         -358572728.821806,
#                -32247804.1300803,
#                20051273.5230284),
#     r=c(         -7136143.23821097,
#                -415951.97857437,
#                -505103.379544365),
#     matrix(c(        -0.999614581251245       ,-0.0277606328990702      ,0.000190291726216941,
#                      0.0277606336661,        -0.999614599341428,      1.39017700274997e-06,
#                      0.000190179795466898,      6.67226010371812e-06,         0.999999981893563), byrow=TRUE, ncol=3),
#     -0.204492811852054,
#     64.184,
#     2.15823965046995e-07,
#     2.6892390799488e-06
# )

pointMassAcceleration <- function(r, s, GM) {
    # difference position vector
    d <- r - s
    a <- -GM * (d/(sqrt(sum(d^2)))^3 + s/(sqrt(sum(s^2)))^3)
    return(a)
}

solarRadiationAcceleration <- function(r, r_earth, r_moon, r_sun, r_sunSSB, 
                                       area, mass, Cr, P0, AU, shm) {
    pccor <- r_moon
    ccor <- r_earth - r_sunSSB
    pscor <- r_moon - r_sun
    sbcor <- r
    bcor <- r - r_sun
    sbpcor <- r - r_moon
    if(shm == "cylindrical") {
        warning("Using cylindrical shadow model for solar radiation pressure")
        nu <- cylindricalShadow(r, r_sun)
    } else {
        nu <- geometricShadow(pccor,ccor,pscor,sbcor,bcor,sbpcor)$lambda
    }
    a <- nu*Cr*(area/mass)*P0*(AU*AU)*bcor/(norm(bcor)^3)
    return(a)
}

dragAcceleration <- function(dens, r, v, T, area, mass, CD, Omega) {
    omega = c(0, 0, Omega)
    r_tod <- T * r
    v_tod <- T * v
    v_rel <- v_tod - c(omega[2]*r_tod[3] - omega[3]*r_tod[2],
                       omega[3]*r_tod[1] - omega[1]*r_tod[3],
                       omega[1]*r_tod[2] - omega[2]*r_tod[1])
    v_abs <- sqrt(sum(v_rel^2))
    a_tod <- -0.5*CD*(area/mass)*dens*v_abs*v_rel
    a <- t(T) %*% a_tod
    return(a)
}

relativity <- function(r, v) {
    r_Sat <- sqrt(sum(r^2))
    v_Sat <- sqrt(sum(v^2))
    a <- GM_Earth/(c_light^2*r_Sat^3)*((4*GM_Earth/r_Sat-v_Sat^2)*r+4*(r %*% v)*v)
    return(a)
}

accel <- function(t, Y, MJD_UTC, satelliteMass, satelliteArea, Cr, Cd) {
    MJD_UTC <- MJD_UTC + t/86400
    IERS_results <- IERS(earthPositions, MJD_UTC, "l")
    timeDiffs_results <- timeDiffs(IERS_results$UT1_UTC, IERS_results$TAI_UTC)
    JD <- MJD_UTC + 2400000.5
    invjday_results <- invjday(MJD_UTC+2400000.5)
    iauCal2jd_results <- iauCal2jd(invjday$year, invjday$month, invjday$day)
    
    TIME <- (60*(60*invjday$hour + invjday$minute) + invjday$sec)/86400
    UTC <- iauCal2jd_results$DATE + TIME
    TT <- UTC + timeDiffs_results$TT_UTC/86400
    TUT <- TIME + IERS_results$UT1_UTC/86400
    UT1 <- iauCal2jd_results$DATE + TUT
    
    PMM <- iauPom00(IERS_results$x_pole, IERS_results$y_pole, iauSp00(iauCal2jd$DJMJD0, TT))
    NPB <- iauPnm06a(iauCal2jd$DJMJD0, TT)

    theta <- iauRz(iauGst06(iauCal2jd$DJMJD0, UT1, iauCal2jd$DJMJD0, TT, NPB), diag(3))
    
    E <- PMM * theta * NPB
    
    MJD_TDB <- Mjday_TDB(TT)
    JPL_ephemerides <- JPL_Eph_DE436(MJD_TDB)
    
    # Acceleration due to Earth, with harmonic effects
    a <- elasticEarthAcceleration(MJD_UTC, JPL_ephemerides$positionSunGeocentric,
                                  JPL_ephemerides$positionMoon, Y[1:3], E,
                                  IERS_results$UT1_UTC, timeDiffs_results$TT_UTC,
                                  IERS_results$x_pole, IERS_results$y_pole)
    # Acceleration due to Sun
    a <- a + pointMassAcceleration(Y[1:3],JPL_ephemerides$positionSunGeocentric,GM_Sun)
    # Acceleration due to Moon
    a <- a + pointMassAcceleration(Y[1:3],JPL_ephemerides$positionMoon,GM_Moon)
    # Accelerations due to planets
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionMercury, GM_Mercury)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionVenus, GM_Venus)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionMars, GM_Mars)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionJupiter, GM_Jupiter)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionSaturn, GM_Saturn)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionUranus, GM_Uranus)    
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionNeptune, GM_Neptune)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionPluto, GM_Pluto)
    # Acceleration due to solar radiation pressure
    a <- a + solarRadiationAcceleration(Y[1:3], JPL_ephemerides$positionEarth, 
                                        JPL_ephemerides$positionMoon, 
                                        JPL_ephemerides$positionSunGeocentric,
                                        JPL_ephemerides$positionSunBarycentric, 
                                        satelliteArea,
                                        satelliteMass,
                                        Cr, solarPressureConst, AU, "geometrical")
    # Acceleration due to atmospheric drag
    Omega <- omegaEarth - 8.43994809e-10*IERS_results$LOD
    dens <- NRLMSISE00(MJD_UTC, E*Y[1:3], timeDiffs_results$UT1_UTC,
                      timeDiffs_results$TT_UTC)
    a <- a + dragAcceleration(dens, Y[1:3], Y[4:6], NPB, satelliteArea, 
                              satelliteMass, Cd, Omega)
    # Relativistic effects
    a <- a + relativity(Y[1:3], Y[4:6])
    ##### ATENCION: ESTA ES LA CLAVE, VERIFICAR QUE EL VALOR DE a PARA LA PRIMERA
    ##### LLAMADA A ESTA FUNCION SEA EL MISMO QUE EN MATLAB
    dY <- matrix(c(Y[4:6], a), byrow=TRUE, ncol=3, nrow=2)
    return(dY)
}

