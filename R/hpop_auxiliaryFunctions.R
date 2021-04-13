MJday <- function(year, month, day, hour=0, min=0, sec=0) {
    y <- year
    m <- month
    b <- 0
    c <- 0
    if (m <= 2) {
        y <- y -1
        m <- m + 12
    }
    if (y < 0) {
        c <- -0.75
    }
    if (year < 1582) {
        # if statement only for subsequent conditionals. TODO CHANGE
    } else if (year > 1582) {
        a <- trunc(y/100)
        b <- 2 - a + floor(a / 4)
    } else if (month < 10) {
        # TODO CHANGE.
    } else if (month > 10) { ## TODO VERIFY introduction of gregorian in 1582
        a <- trunc(y/100)
        b <- 2 - a + floor(a / 4)
    } else if (day <= 4) {
        # TODO CHANGE
    } else if (day > 14) {
        a <- trunc(y/100)
        b <- 2 - a + floor(a / 4)
    } else {
        stop("Please enter a valid calendar date")
    }
    jd <- trunc(365.25*y + c) + trunc(30.6001 * (m+1))
    jd <- jd + day + b + 1720994.5
    jd <- jd + (hour+min/60+sec/3600)/24
    Mjd <- jd - 2400000.5
    return(Mjd)
}


IERS <- function(eop,Mjd_UTC,interp="n") {
    if(interp == "l") {
        mjd <- floor(Mjd_UTC)
        i <- which(mjd == eop[, 4])[1]
        preeop <- eop[i, ]
        nexteop <- eop[i+1, ]
        mfme <- 1440*(Mjd_UTC - floor(Mjd_UTC))
        fixf <- mfme/1440
        # setting IERS Earth rotation parameters 
        # (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole <- preeop[5]+(nexteop[5]-preeop[5])*fixf
        y_pole <- preeop[6]+(nexteop[6]-preeop[6])*fixf
        UT1_UTC <- preeop[7]+(nexteop[7]-preeop[7])*fixf
        LOD <- preeop[8]+(nexteop[8]-preeop[8])*fixf
        dpsi <- preeop[9]+(nexteop[9]-preeop[9])*fixf
        deps <- preeop[10]+(nexteop[10]-preeop[10])*fixf
        dx_pole <- preeop[11]+(nexteop[11]-preeop[11])*fixf
        dy_pole = preeop[12]+(nexteop[12]-preeop[12])*fixf
        TAI_UTC = preeop[13]
        
        x_pole <- x_pole/const_Arcs  # Pole coordinate (rad)
        y_pole <- y_pole/const_Arcs  # Pole coordinate (rad)
        dpsi <- dpsi/const_Arcs
        deps <- deps/const_Arcs
        dx_pole <- dx_pole/const_Arcs  # Pole velocity? (rad)
        dy_pole <- dy_pole/const_Arcs  # Pole velocity? (rad)
    } else if(interp == "n") {
        mjd = (floor(Mjd_UTC))
        i <- which(mjd == eop[, 4])[1]
        eop <- eop[i, ]
        # setting IERS Earth rotation parameters 
        # (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole <- eop[5]/const_Arcs # Pole coordinate (rad)
        y_pole <- eop[6]/const_Arcs # Pole coordinate (rad)
        UT1_UTC <- eop[7] # UT1-UTC time difference (s)
        LOD <- eop[8] # Length of day (s)
        dpsi <- eop[9]/const_Arcs
        deps <- eop[10]/const_Arcs
        dx_pole <- eop[11]/const_Arcs # Pole velocity? (rad)
        dy_pole <- eop[12]/const_Arcs # Pole velocity? (rad)
        TAI_UTC <- eop[13] # TAI-UTC time difference (s)
    }
    return(list(
        x_pole=x_pole,
        y_pole=y_pole,
        UT1_UTC=UT1_UTC,
        LOD=LOD,
        dpsi=dpsi,
        deps=deps,
        dx_pole=dx_pole,
        dy_pole=dy_pole,
        TAI_UTC=TAI_UTC
    ))
}

timeDiffs <- function(UT1_UTC,TAI_UTC) {
    TT_TAI <- 32.184 # TT-TAI time difference (s)
    GPS_TAI <- -19.0 # GPS-TAI time difference (s)
    TT_GPS <- TT_TAI-GPS_TAI # TT-GPS time difference (s)
    TAI_GPS <- -GPS_TAI # TAI-GPS time difference (s)
    UT1_TAI <- UT1_UTC-TAI_UTC # UT1-TAI time difference (s)
    UTC_TAI <- -TAI_UTC # UTC-TAI time difference (s)
    UTC_GPS <- UTC_TAI-GPS_TAI # UTC_GPS time difference (s)
    UT1_GPS <- UT1_TAI-GPS_TAI # UT1-GPS time difference (s)
    TT_UTC <- TT_TAI-UTC_TAI # TT-UTC time difference (s)
    GPS_UTC <- GPS_TAI-UTC_TAI # GPS-UTC time difference (s)
    return(list(
        TT_TAI=TT_TAI,
        GPS_TAI=GPS_TAI,
        TT_GPS=TT_GPS,
        TAI_GPS=TAI_GPS,
        UT1_TAI=UT1_TAI,
        UTC_TAI=UTC_TAI,
        UTC_GPS=UTC_GPS,
        UT1_GPS=UT1_GPS,
        TT_UTC=TT_UTC,
        GPS_UTC=GPS_UTC
    ))
}

invjday <- function(jd) {
    z <- trunc(jd + 0.5)
    fday <- jd + 0.5 - z
    if (fday < 0) {
        fday <- fday + 1
        z <- z - 1
    }
    if (z < 2299161) {
        a <- z
    } else {
        alpha <- floor((z-1867216.25) / 36524.25)
        a <- z + 1 + alpha - floor(alpha/4)
    }
    b <- a + 1524
    c <- trunc((b - 122.1) / 365.25)
    d <- trunc(365.25 * c)
    e <- trunc((b-d) / 30.6001)
    day <- b - d - trunc(30.6001 * e) + fday
    if (e < 14) {
        month <- e - 1
    } else {
        month <- e - 13
    }
    if (month > 2) {
        year <- c - 4716
    } else {
        year <- c - 4715
    }
    hour <- abs(day-floor(day))*24
    min <- abs(hour-floor(hour))*60
    sec <- abs(min-floor(min))*60
    day <- floor(day)
    hour <- floor(hour)
    min <- floor(min)
    return(list(
        year=year, 
        month=month, 
        day=day, 
        hour=hour, 
        min=min, 
        sec=sec
    ))
}

ECEFtoECI <- function(MJD_UTC, Y0) {
    IERS_results <- IERS(earthPositions, MJD_UTC)
    timeDiffs_results <- timeDiffs(IERS_results$UT1_UTC, IERS_results$TAI_UTC)
    invjday_results <- invjday(MJD_UTC+2400000.5)
    iauCal2jd_results <- iauCal2jd(invjday_results$year, invjday_results$month, invjday_results$day)
    TIME <- (60*(60*invjday_results$hour + invjday_results$minute) + invjday_results$sec)/86400
    UTC <- iauCal2jd_results$DATE + TIME
    TT <- UTC + timeDiffs_results$TT_UTC/86400
    TUT <- TIME + IERS_results$UT1_UTC/86400
    UT1 <- iauCal2jd_results$DATE + TUT
    NPB <- iauPnm06a(iauCal2jd_results$DJMJD0, TT)
    theta <- iauRz(iauGst06(iauCal2jd_results$DJMJD0, UT1, iauCal2jd_results$DJMJD0, TT, NPB), diag(3))
    PMM <- iauPom00(IERS_results$x_pole, IERS_results$y_pole, iauSp00(iauCal2jd_results$DJMJD0, TT))
    S <- matrix(c(0, 1, 0,
                  1, 0, 0,
                  0, 0, 0),
                byrow = TRUE)
    omega <- 7292115.8553e-11+4.3e-15*((MJD_UTC-MJD_J2000)/36525)
    dTheta <- omega*S*theta
    U <- PMM*theta*NPB
    dU <- PMM*dTheta*NPB
    r <- t(U)*t(Y0(1:3))
    v <- t(U)*t(Y0(4:6)) + t(dU)*t(Y0(1:3))
    return (list(
        position=r,
        velocity=v
    ))
}

Mjday_TDB <- function(Mjd_TT) {
    # Given Modified julian date (TT), compute Modified julian date (TDB)
    T_TT <- (Mjd_TT - 51544.5)/36525
    Mjd_TDB <- Mjd_TT + ( 0.001658*sin(628.3076*T_TT +6.2401)+
                              0.000022*sin(575.3385*T_TT+4.2970)+
                              0.000014*sin(1256.6152*T_TT + 6.1969)+
                              0.000005*sin(606.9777*T_TT+4.0212)+  
                              0.000005*sin(52.9691*T_TT+0.4444) +   
                              0.000002*sin(21.3299*T_TT+5.5431)+   
                              0.000010*sin(628.3076*T_TT+4.2490) )/86400
    return(Mjd_TDB)
}

cheb3D <- function(t, N, Ta, Tb, Cx, Cy, Cz) {
    if ((t<Ta) | (t>Tb)) {
        stop('Time out of range for Chebyshev approximation')
    }
    tau <- (2*t-Ta-Tb)/(Tb-Ta)
    f1 <- rep(0, 3)
    f2 <- f1
    for(i in N:2) {
        old_f1 <- f1
        f1 <- 2*tau*f1-f2+c(Cx[[i]],Cy[[i]],Cz[[i]])
        f2 <- old_f1
    }
    cheb_approximation <- tau*f1-f2+c(Cx[[1]],Cy[[1]],Cz[[1]])
    return(cheb_approximation)
}

JPL_Eph_DE436 <- function(Mjd_TDB) {
    # calculate equatorial position of sun, moon, and nine major planets 
    # using JPL Ephemerides
    JD <- Mjd_TDB + 2400000.5
    i <- which(JD > DE436coeffs[, 1] & JD <= DE436coeffs[, 2])[1]
    current_DE436coeffs <- DE436coeffs[i, ]
    t1 <- current_DE436coeffs[[1]] - 2400000.5
    dt <- Mjd_TDB - t1
    indexes <- c(231, 244, 257, 270)
    Cx_earth <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_earth <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_earth <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    indexes <- indexes+39
    Cx <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    Cx_earth <- c(Cx_earth, Cx)
    Cy_earth <- c(Cy_earth, Cy)
    Cz_earth <- c(Cz_earth, Cz)
    if(dt >= 0 & dt <= 16) {
        j <- 0
        Mjd0 <- t1
    } else if (dt > 16 & dt <= 32) {
        j <- 1
        Mjd0 <- t1 + 16*j
    }
    r_earth <- 1000*cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, 
                           Cx_earth[(13*j+1):(13*j+13)], 
                           Cy_earth[(13*j+1):(13*j+13)], 
                           Cz_earth[(13*j+1):(13*j+13)])
    indexes <- c(441, 454, 467, 480)
    Cx_moon <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_moon <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_moon <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    for (i in 1:7) {
        indexes <- indexes + 39
        Cx <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
        Cy <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
        Cz <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
        Cx_moon <- c(Cx_moon,Cx)
        Cy_moon <- c(Cy_moon,Cy)
        Cz_moon <- c(Cz_moon,Cz) 
    }
    if (dt>=0 & dt<=4) {
        j <- 0
        Mjd0 <- t1
    } else if(dt>4 & dt<=8) {
        j <- 1
        Mjd0 <- t1+4*j
    } else if(dt>8 & dt<=12) {
        j <- 2
        Mjd0 <- t1+4*j
    } else if(dt>12 & dt<=16) {
        j <- 3
        Mjd0 <- t1+4*j
    } else if(dt>16 & dt<=20) {
        j <- 4
        Mjd0 <- t1+4*j
    } else if(dt>20 & dt<=24) {
        j <- 5
        Mjd0 <- t1+4*j
    } else if(dt>24 & dt<=28) {
        j <- 6
        Mjd0 <- t1+4*j
    } else if(dt>28 & dt<=32) {
        j <- 7
        Mjd0 <- t1+4*j
    }
    r_moon <- 1000*cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, 
                          Cx_moon[(13*j+1):(13*j+13)],
                          Cy_moon[(13*j+1):(13*j+13)], 
                          Cz_moon[(13*j+1):(13*j+13)])
    indexes <- c(753, 764, 775, 786)
    Cx_sun <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_sun <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_sun <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    indexes <- indexes + 33
    Cx <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    Cx_sun <- c(Cx_sun, Cx)
    Cy_sun <- c(Cy_sun, Cy)
    Cz_sun <- c(Cz_sun, Cz)
    if(dt >= 0 & dt <= 16) {
        j <- 0
        Mjd0 <- t1
    } else if (dt > 16 & dt <= 32) {
        j <- 1
        Mjd0 <- t1 + 16*j
    }
    r_sun <- 1000*cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, 
                         Cx_sun[(11*j+1):(11*j+11)],
                         Cy_sun[(11*j+1):(11*j+11)], 
                         Cz_sun[(11*j+1):(11*j+11)])
    indexes <- c(3, 17, 31, 45)
    Cx_mercury <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_mercury <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_mercury <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    for (i in 1:3) {
        indexes <- indexes + 42
        Cx <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
        Cy <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
        Cz <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
        Cx_mercury <- c(Cx_mercury,Cx)
        Cy_mercury <- c(Cy_mercury,Cy)
        Cz_mercury <- c(Cz_mercury,Cz) 
    }
    if (dt>=0 & dt<=8) {
        j <- 0
        Mjd0 <- t1
    } else if(dt>8 & dt<=16) {
        j <- 1
        Mjd0 <- t1+8*j
    } else if(dt>16 & dt<=24) {
        j <- 2
        Mjd0 <- t1+8*j
    } else if(dt>24 & dt<=32) {
        j <- 3
        Mjd0 <- t1+8*j
    }
    r_mercury <- 1e3*cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, 
                            Cx_mercury[(14*j+1):(14*j+14)],
                            Cy_mercury[(14*j+1):(14*j+14)], 
                            Cz_mercury[(14*j+1):(14*j+14)])
    indexes <- c(171, 181, 191, 201)
    Cx_venus <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_venus <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_venus <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    indexes <- indexes+30
    Cx <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    Cx_venus <- c(Cx_venus, Cx)
    Cy_venus <- c(Cy_venus, Cy)
    Cz_venus <- c(Cz_venus, Cz)
    if(dt >= 0 & dt <= 16) {
        j <- 0
        Mjd0 <- t1
    } else if (dt > 16 & dt <= 32) {
        j <- 1
        Mjd0 <- t1 + 16*j
    }
    r_venus <- 1000*cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, 
                           Cx_venus[(10*j+1):(10*j+10)],
                           Cy_venus[(10*j+1):(10*j+10)], 
                           Cz_venus[(10*j+1):(10*j+10)])
    indexes <- c(309, 320, 331, 342)
    Cx_mars <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_mars <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_mars <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    j <- 0
    Mjd0 <- t1
    r_mars = 1000*cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, 
                         Cx_mars[(11*j+1):(11*j+11)],
                         Cy_mars[(11*j+1):(11*j+11)], 
                         Cz_mars[(11*j+1):(11*j+11)])
    indexes <- c(342, 350, 358, 366)
    Cx_jupiter <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_jupiter <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_jupiter <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    j <- 0
    Mjd0 <- t1
    r_jupiter <- 1000*cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, 
                             Cx_jupiter[(8*j+1):(8*j+8)],
                             Cy_jupiter[(8*j+1):(8*j+8)], 
                             Cz_jupiter[(8*j+1):(8*j+8)])
    indexes <- c(366, 373, 380, 387)
    Cx_saturn <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_saturn <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_saturn <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    j <- 0
    Mjd0 <- t1
    r_saturn <- 1000*cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, 
                            Cx_saturn[(7*j+1):(7*j+7)],
                            Cy_saturn[(7*j+1):(7*j+7)], 
                            Cz_saturn[(7*j+1):(7*j+7)])
    indexes <- c(387, 393, 399, 405)
    Cx_uranus <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_uranus <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_uranus <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    j <- 0
    Mjd0 <- t1
    r_uranus <- 1000*cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, 
                            Cx_uranus[(6*j+1):(6*j+6)],
                            Cy_uranus[(6*j+1):(6*j+6)],
                            Cz_uranus[(6*j+1):(6*j+6)])
    indexes <- c(405, 411, 417, 423)
    Cx_neptune <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_neptune <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_neptune <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    j <- 0
    Mjd0 <- t1
    r_neptune <- 1000*cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, 
                             Cx_neptune[(6*j+1):(6*j+6)],
                             Cy_neptune[(6*j+1):(6*j+6)], 
                             Cz_neptune[(6*j+1):(6*j+6)])
    indexes <- c(423, 429, 435, 441)
    Cx_pluto <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_pluto <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_pluto <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    j <- 0
    Mjd0 <- t1
    r_pluto <- 1000*cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, 
                           Cx_pluto[(6*j+1):(6*j+6)],
                           Cy_pluto[(6*j+1):(6*j+6)], 
                           Cz_pluto[(6*j+1):(6*j+6)])
    indexes <- c(819, 829, 839)
    Cx_nutations <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_nutations <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    for (i in 1:3) {
        indexes <- indexes + 20
        Cx <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
        Cy <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
        Cx_nutations <- c(Cx_nutations,Cx)
        Cy_nutations <- c(Cy_nutations,Cy)
    }
    if (dt>=0 & dt<=8) {
        j <- 0
        Mjd0 <- t1
    } else if(dt>8 & dt<=16) {
        j <- 1
        Mjd0 <- t1+8*j
    } else if(dt>16 & dt<=24) {
        j <- 2
        Mjd0 <- t1+8*j
    } else if(dt>24 & dt<=32) {
        j <- 3
        Mjd0 <- t1+8*j
    }
    nutations <- cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, 
                        Cx_nutations[(10*j+1):(10*j+10)],
                        Cy_nutations[(10*j+1):(10*j+10)],
                        rep(0,10))
    indexes <- c(899, 909, 919, 929)
    Cx_libration <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
    Cy_libration <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
    Cz_libration <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
    for (i in 1:3) {
        indexes <- indexes + 30
        Cx <- current_DE436coeffs[indexes[1]:(indexes[2]-1)]
        Cy <- current_DE436coeffs[indexes[2]:(indexes[3]-1)]
        Cz <- current_DE436coeffs[indexes[3]:(indexes[4]-1)]
        Cx_libration <- c(Cx_libration,Cx)
        Cy_libration <- c(Cy_libration,Cy)
        Cz_libration <- c(Cz_libration,Cz) 
    }
    if (dt>=0 & dt<=8) {
        j <- 0
        Mjd0 <- t1
    } else if(dt>8 & dt<=16) {
        j <- 1
        Mjd0 <- t1+8*j
    } else if(dt>16 & dt<=24) {
        j <- 2
        Mjd0 <- t1+8*j
    } else if(dt>24 & dt<=32) {
        j <- 3
        Mjd0 <- t1+8*j
    }
    librations <- cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, 
                         Cx_libration[(10*j+1):(10*j+10)], 
                         Cy_libration[(10*j+1):(10*j+10)], 
                         Cz_libration[(10*j+1):(10*j+10)])
    r_earth <- r_earth-EMRAT1*r_moon
    r_mercury <- -r_earth+r_mercury
    r_venus <- -r_earth+r_venus
    r_mars <- -r_earth+r_mars
    r_jupiter <- -r_earth+r_jupiter
    r_saturn <- -r_earth+r_saturn
    r_uranus <- -r_earth+r_uranus
    r_neptune <- -r_earth+r_neptune
    r_pluto <- -r_earth+r_pluto
    r_sunSSB <- r_sun
    r_sun <- -r_earth+r_sun
    return(list(
        positionMercury = r_mercury,
        positionVenus = r_venus,
        positionEarth = r_earth,
        positionMars = r_mars,
        positionJupiter = r_jupiter,
        positionSaturn = r_saturn,
        positionUranus = r_uranus,
        positionNeptune = r_neptune,
        positionPluto = r_pluto,
        positionMoon = r_moon,
        positionSunGeocentric = r_sun,
        positionSunBarycentric = r_sunSSB
    ))
}

CartesianToPolar <- function(cartesianVector) {
    rho2 <- cartesianVector[1]^2 + cartesianVector[2]^2
    # rho - length of projection on XY plane
    rho <- sqrt(rho2)
    # Norm - r
    r <- sqrt(rho2 + cartesianVector[3]^2)
    # Azimuth - phi
    if ((cartesianVector[1]==0) & (cartesianVector[2]==0)) {
        phi <- 0
    } else {
        phi <- atan2(cartesianVector[2], cartesianVector[1])
    }
    if (phi < 0) {
        phi <- phi + 2*pi
    }
    # Elevation - theta
    if ((cartesianVector[3]==0) & (rho==0)) {
        theta <- 0
    } else {
        theta <- atan2(cartesianVector[3], rho)
    }
    return(c(phi, theta, r))
}

legendre <- function(n, m, angle) {
    pnm <- matrix(0, nrow=n+1, ncol=m+1)
    dpnm <- matrix(0, nrow=n+1, ncol=m+1)
    pnm[1,1] <- 1
    pnm[2,2] <- sqrt(3)*cos(angle)
    dpnm[2,2] <- -sqrt(3)*sin(angle)
    # diagonal coefficients
    for (i in 2:n) {
        pnm[i+1, i+1] <- sqrt((2*i+1)/(2*i))*cos(angle)*pnm[i,i]
        dpnm[i+1, i+1] <- sqrt((2*i+1)/(2*i))*((cos(angle)*dpnm[i,i]) - (sin(angle)*pnm[i,i]))
    }
    # horizontal first step coefficients
    for (i in 1:n) {
        pnm[i+1, i] <- sqrt(2*i+1)*sin(angle)*pnm[i,i]
        dpnm[i+1, i] <- sqrt(2*i+1)*((cos(angle)*pnm[i,i])+(sin(angle)*dpnm[i,i]))
    }
    # horizontal second step coefficients
    j <- 0
    k <- 2
    while(j <= m & k <= n) {
        for(i in k:n) {
            pnm[i+1,j+1] <- sqrt((2*i+1)/((i-j)*(i+j)))*
                ((sqrt(2*i-1)*sin(angle)*pnm[i,j+1]) - 
                     (sqrt(((i+j-1)*(i-j-1))/(2*i-3))*pnm[i-1,j+1]))
            dpnm[i+1,j+1] <- sqrt((2*i+1)/((i-j)*(i+j)))*
                ((sqrt(2*i-1)*sin(angle)*dpnm[i,j+1]) + 
                     (sqrt(2*i-1)*cos(angle)*pnm[i,j+1]) - 
                     (sqrt(((i+j-1)*(i-j-1))/(2*i-3))*dpnm[i-1,j+1]))
        }
        j <- j + 1
        k <- k + 1
    }
    return(list(
        normLegendreValues=pnm, 
        normLegendreDerivativeValues=dpnm))
}

cylindricalShadow <- function(r, positionSun) {
    sunDirection <- positionSun/sqrt(sum(positionSun^2))
    satelliteProjection <- r%*%sunDirection
    if(satelliteProjection > 0 | sqrt(sum((r - satelliteProjection*sunDirection)^2))) {
        nu <- 1
    } else {
        nu <- 0
    }
    return(nu)
}

calculateLambda <- function(rs, rp, sep) {
    # rs : apparent radius of sun as viewed from satellite (radians)
    # rp : apparent radius of eclipsing body as viewed from satellite (radians)
    # sep: apparent separation of the center of the Sun and eclipsing body (radians)
    if (rs + rp <= sep) {
        lambda <- 1
    } else if ((rp - rs) >= sep) {
        lambda <- 0
    } else if (sep > (rs-rp)){
        if (rs > rp) { # assign r1 to be smaller disc, r2 larger disc
            r1 <- rp
            r2 <- rs
        } else {
            r1 <- rs
            r2 <- rp
        }
        phi <- acos((r1*r1+sep*sep-r2*r2)/(2*r1*sep))
        if (phi < 0) {
            phi <- phi + pi
        }
        if (r2/r1 > 5) {
            hgt <- sqrt(r1^2-(sep-r2)^2)
            area2 <- hgt*(sep-r2)
            area3 <- 0
        } else {
            hgt <- r1*sin(phi)
            theta <- asin(hgt/r2)
            area2 <- sep*hgt
            area3 <- theta*r2^2
        }
        area1 <- (pi-phi)*r1^2
        ari <- area1 + area2 - area3 # non-overlapped section of smaller disc
        area1 <- pi*rs^2
        if (rs>rp) {
            area2 <- pi*rp^2
            lambda <- (area1+ari-area2)/area1
        } else {
            lambda <- ari/area1
        }
    } else {
        lambda <- (rs^2-rp^2)/rs^2 
    }
    return(lambda)
}

geometricShadow <- function(pccor,ccor,pscor,sbcor,bcor,sbpcor) {
    lambda <- 1 # lambda ranges 0 to 1, 1 = no shadow, 0 = full shadow
    eclipse <- "none" # eclipse type
    ubcor <- rep(0, 3)
    rb <- sqrt(sum(bcor^2))
    rc <- sqrt(sum(ccor^2))
    if (rb > rc) {
        ubcor <- bcor/rb # direction of satellite relative to sun
        sepp <- c(sbcor[2]*ubcor[3] - sbcor[3]*ubcor[2],
                  sbcor[3]*ubcor[1] - sbcor[1]*ubcor[3],
                  sbcor[1]*ubcor[2] - sbcor[2]*ubcor[1]) # perpendicular
        rsbx <- sbcor%*%ubcor # projection of sbcor along ubcor
        rs <- R_sun/rb # apparent radius of Sun from satellite
        rp <- r_ref/rsbx # apparent radius of Earth from satellite
        sep <- sqrt(sum(sepp^2))/rsbx # apparent separation between Sun and Earth
        lambda <- calculateLambda(rs,rp,sep)
    }
    if (lambda < 1) {
        eclipse <-"E" # if lambda <1 here, there is at least partial shadowing by Earth = earth eclipse
    }  else { # checking Moon eclipse
        pscor <- pccor + ccor
        sbpcor <- sbcor - pccor
        rps <- sqrt(pscor[1]^2 + pscor[2]^2 + pscor[3]^2)
        if (rb > rps) {
            ubcor <- bcor/rb
            rsbx <- sbpcor %*% ubcor
            rs <- R_sun/rb
            rp <- R_moon/rsbx
            sepp <- c(sbpcor[2]*ubcor[3] - sbpcor[3]*ubcor[2],
                      sbpcor[3]*ubcor[1] - sbpcor[1]*ubcor[3],
                      sbpcor[1]*ubcor[2] - sbpcor[2]*ubcor[1])
            sep <- sqrt(sum(sepp^2))/rsbx # apparent separation between Sun and Moon
            lambda <- calculateLambda(rs, rp, sep)
            if (lambda < 1) {
                eclipse <- "M" # if lambda <1, moon eclipse
            }
        }
    }
    return(list(
        lambda = lambda,
        eclipseType = eclipse
    ))
}

fractionalDays <- function(year, month, day, hour, min, sec) {
    monthDays <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    if((year%%4 == 0) & (year%%100 == 0) & (year%%400 != 0)) {
        monthDays[2] <- 29
    }
    days <- sum(monthDays[1:(month - 1)]) + day + hour/24 + min/1440 + sec/86400
    return(days)
}

ECEFtoGeodetic <- function(r) { # Verify why it gives different results than other function
    R_equ <- 6378.1363e3
    f <- 1/298.257223563
    epsRequ <- 2.2204e-16*R_equ
    e2 <- f*(2 - f)
    X <- r[1]
    Y <- r[2]
    Z <- r[3]
    rho2 <- X^2 + Y^2
    dZ <- e2*Z
    dZ_new <- dZ + 10*epsRequ
    while(TRUE) {
        ZdZ <- Z + dZ
        Nh <- sqrt(rho2 + ZdZ^2);
        SinPhi <- ZdZ/Nh
        N <- R_equ / sqrt(1-e2*SinPhi^2)
        dZ_new <- N*e2*SinPhi
        if (abs(dZ-dZ_new) < epsRequ) {
            break
        }
        dZ <- dZ_new
    }
    lon <- atan2(Y, X)
    lat <- atan2(ZdZ, sqrt(rho2))
    h <- Nh - N
    return(list(
        lon = lon,
        lat = lat,
        h = h
    ))
}

ephemeris <- function(Y0, n_step, step) {
    
}