IERS <- function(eop,Mjd_UTC,interp="n") {
    if(interp == "l") {
        mjd <- floor(Mjd_UTC)
        i <- which(mjd == eop[, 4])[1]
        if(is.na(i)) {
            warning(strwrap("No Earth Orientation Parameters found for the
                            specified date. Try to obtain latest space data
                            by running getLatestSpaceData()"))
        }
        preeop <- as.numeric(eop[i, ])
        nexteop <- as.numeric(eop[i+1, ])
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
        dx_pole <- dx_pole/const_Arcs  # Pole velocity (rad)
        dy_pole <- dy_pole/const_Arcs  # Pole velocity (rad)
    } else if(interp == "n") {
        mjd = (floor(Mjd_UTC))
        i <- which(mjd == eop[, 4])[1]
        eop <- as.numeric(eop[i, ])
        # IERS Earth rotation parameters 
        # (UT1-UTC in s, TAI-UTC in s, x in rads, y in rads)
        x_pole <- eop[5]/const_Arcs # Pole coordinate (rad)
        y_pole <- eop[6]/const_Arcs # Pole coordinate (rad)
        UT1_UTC <- eop[7] # UT1-UTC time difference (s)
        LOD <- eop[8] # Length of day (s)
        dpsi <- eop[9]/const_Arcs
        deps <- eop[10]/const_Arcs
        dx_pole <- eop[11]/const_Arcs # Pole velocity (rad)
        dy_pole <- eop[12]/const_Arcs # Pole velocity (rad)
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
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "l")
    timeDiffs_results <- timeDiffs(IERS_results$UT1_UTC, IERS_results$TAI_UTC)
    invjday_results <- invjday(MJD_UTC+2400000.5)
    iauCal2jd_results <- iauCal2jd(invjday_results$year, invjday_results$month, invjday_results$day)
    TIME <- (60*(60*invjday_results$hour + invjday_results$min) + invjday_results$sec)/86400
    UTC <- iauCal2jd_results$DATE + TIME
    TT <- UTC + timeDiffs_results$TT_UTC/86400
    TUT <- TIME + IERS_results$UT1_UTC/86400
    UT1 <- iauCal2jd_results$DATE + TUT
    NPB <- iauPnm06a(iauCal2jd_results$DJMJD0, TT, IERS_results$dpsi, IERS_results$deps)
    theta <- iauRz(iauGst06(iauCal2jd_results$DJMJD0, UT1, iauCal2jd_results$DJMJD0, TT, NPB), diag(3))
    PMM <- iauPom00(IERS_results$x_pole, IERS_results$y_pole, iauSp00(iauCal2jd_results$DJMJD0, TT))
    S <- matrix(c(0, 1, 0,
                  -1, 0, 0,
                  0, 0, 0),
                nrow = 3, byrow = TRUE)
    # omega <- 7292115.8553e-11+4.3e-15*((MJD_UTC-MJD_J2000)/36525)
    omega <- omegaEarth - 8.43994809e-10*IERS_results$LOD
    dTheta <- omega*S%*%theta
    U <- PMM%*%theta%*%NPB
    dU <- PMM%*%dTheta%*%NPB
    r <- t(U)%*%Y0[1:3]
    v <- t(U)%*%Y0[4:6] + t(dU)%*%Y0[1:3]
    return (list(
        position=r,
        velocity=v
    ))
}

ECItoECEF <- function(MJD_UTC, Y0) {
    # This and the inverse follow the protocol established by IERS 2010 
    # technical note to convert between GCRF and ITRF
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "l")
    timeDiffs_results <- timeDiffs(IERS_results$UT1_UTC, IERS_results$TAI_UTC)
    invjday_results <- invjday(MJD_UTC+2400000.5)
    iauCal2jd_results <- iauCal2jd(invjday_results$year, invjday_results$month, invjday_results$day)
    TIME <- (60*(60*invjday_results$hour + invjday_results$min) + invjday_results$sec)/86400
    UTC <- iauCal2jd_results$DATE + TIME
    TT <- UTC + timeDiffs_results$TT_UTC/86400
    TUT <- TIME + IERS_results$UT1_UTC/86400
    UT1 <- iauCal2jd_results$DATE + TUT
    NPB <- iauPnm06a(iauCal2jd_results$DJMJD0, TT, IERS_results$dpsi, IERS_results$deps)
    theta <- iauRz(iauGst06(iauCal2jd_results$DJMJD0, UT1, iauCal2jd_results$DJMJD0, TT, NPB), diag(3))
    PMM <- iauPom00(IERS_results$x_pole, IERS_results$y_pole, iauSp00(iauCal2jd_results$DJMJD0, TT))
    S <- matrix(c(0, 1, 0,
                  -1, 0, 0,
                  0, 0, 0),
                nrow = 3, byrow = TRUE)
    # omega <- 7292115.8553e-11+4.3e-15*((MJD_UTC-MJD_J2000)/36525)
    omega <- omegaEarth - 8.43994809e-10*IERS_results$LOD
    dTheta <- omega*S%*%theta
    U <- PMM%*%theta%*%NPB
    dU <- PMM%*%dTheta%*%NPB
    r <- U%*%Y0[1:3]
    v <- U%*%Y0[4:6] + dU%*%Y0[1:3]
    return (list(
        position=r,
        velocity=v
    ))
}

Mjday_TDB <- function(Mjd_TT) {
    # Given Modified julian date (TT compatible), compute Modified julian date (TDB compatible)
    T_TT <- (Mjd_TT - 51544.5)/36525 # Julian centuries in TT
    Mjd_TDB <- Mjd_TT + ( 0.001658 * sin(628.3076 * T_TT + 6.2401) +
                              0.000022 * sin(575.3385 * T_TT + 4.2970) +
                              0.000014 * sin(1256.6152 * T_TT + 6.1969) +
                              0.000005 * sin(606.9777 * T_TT + 4.0212)+
                              0.000005 * sin(52.9691 * T_TT + 0.4444) +
                              0.000002 * sin(21.3299 * T_TT + 5.5431)+
                              0.000010 * sin(628.3076 * T_TT + 4.2490) )/86400
    # Mjd_TDB <- Mjd_TT + iauDtdb(Mjd_TT)
    return(Mjd_TDB)
}

clenshaw_old <- function(t, N, Ta, Tb, Coeffs, derivatives=FALSE) {
    # See page 75 from chapter 3 of Vallado's book. 
    # Uses Clenshaw algorithm for sum of Chebyshev polynomial
    if ((t<Ta) | (t>Tb)) {
        stop('Time out of range for Chebyshev approximation')
    }
    tau <- (2*t-Ta-Tb)/(Tb-Ta)
    f1 <- rep(0, ncol(Coeffs))
    f2 <- f1
    if(derivatives) {
        df1 <- rep(0, ncol(Coeffs))
        df2 <- df1
    }
    for(i in N:2) {
        old_f1 <- f1
        f1 <- 2*tau*f1-f2+Coeffs[i, ]
        f2 <- old_f1
    }
    chebyshevSum <- tau*f1-f2+ Coeffs[1, ]
    return(chebyshevSum)
}

clenshaw <- function(t, N, Ta, Tb, Coeffs, derivatives=FALSE) {
    # See page 75 from chapter 3 of Vallado's book. 
    # Uses Clenshaw algorithm for sum of Chebyshev polynomial
    if ((t<Ta) | (t>Tb)) {
        stop('Time out of range for Chebyshev approximation')
    }
    tau <- (2*t-Ta-Tb)/(Tb-Ta)
    f0 <- rep(0, ncol(Coeffs))
    f1 <- f0
    if(derivatives) {
        df0 <- rep(0, ncol(Coeffs))
        df1 <- df0
    }
    for(i in N:2) {
        f2 <- f1
        f1 <- f0
        f0 <- 2*tau*f1-f2+Coeffs[i, ]
        if(derivatives) {
            df2 <- df1
            df1 <- df0
            df0 <- 2*f1 + 2*tau*df1 - df2
        }
    }
    chebyshevSum <- tau*f0 - f1 + Coeffs[1, ]
    if(derivatives) {
        chebyshevSumDerivative <- (tau*df0 - df1 + f0)/((Tb-Ta)/2 * 86400)
        return(list(
            chebyshevSum=chebyshevSum,
            chebyshevSumDerivative=chebyshevSumDerivative
        ))
    } else {
        return(chebyshevSum)
    }
}

clenshawDerivative <- function(t, N, Ta, Tb, Coeffs) {
    # See page 75 from chapter 3 of Vallado's book. 
    # Uses Clenshaw algorithm for sum of Chebyshev polynomial
    if ((t<Ta) | (t>Tb)) {
        stop('Time out of range for Chebyshev approximation')
    }
    tau <- (2*t-Ta-Tb)/(Tb-Ta)
    f1 <- rep(0, ncol(Coeffs))
    f2 <- f1
    for(i in N:2) {
        old_f1 <- f1
        f1 <- 2*tau*f1-f2+Coeffs[i, ]
        f2 <- old_f1
    }
    chebyshevSum <- tau*f1-f2+ Coeffs[1, ]
    return(chebyshevSum)
}

JPLephemeridesDE440 <- function(MJD_TDB, centralBody = "Earth", derivativesOrder=2) {
    JD <- MJD_TDB + 2400000.5
    ephStartJD <- asteRiskData::DE440coeffs[1,1]
    targetRow <- (JD - ephStartJD)%/%32 + 1
    coeffs <- asteRiskData::DE440coeffs[targetRow, ]
    rowStartJD <- coeffs[1]
    rowStartMJD <- rowStartJD - 2400000.5
    dt <- JD - rowStartJD
    ## Sun ##
    subInt <- dt %/% 16 + 1
    subIntStartMJD <- rowStartMJD + 16*(subInt - 1)
    idxs <- seq(753 + 33*(subInt-1), by=11, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+10)]
    Cy <- coeffs[idxs[2]:(idxs[2]+10)]
    Cz <- coeffs[idxs[3]:(idxs[3]+10)]
    clenshawSun <- cbind(
        clenshawAllDerivatives(MJD_TDB, 11, subIntStartMJD, subIntStartMJD+16, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 11, subIntStartMJD, subIntStartMJD+16, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 11, subIntStartMJD, subIntStartMJD+16, Cz, derivativesOrder)
    )
    # clenshawSun <- clenshaw(MJD_TDB, 11, subIntStartMJD, subIntStartMJD+16, cbind(Cx, Cy, Cz), derivatives)
    ## Mercury ##
    subInt <- dt %/% 8 + 1
    subIntStartMJD <- rowStartMJD + 8*(subInt - 1)
    idxs <- seq(3 + 42*(subInt-1), by=14, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+13)]
    Cy <- coeffs[idxs[2]:(idxs[2]+13)]
    Cz <- coeffs[idxs[3]:(idxs[3]+13)]
    clenshawMercury <- cbind(
        clenshawAllDerivatives(MJD_TDB, 14, subIntStartMJD, subIntStartMJD+8, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 14, subIntStartMJD, subIntStartMJD+8, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 14, subIntStartMJD, subIntStartMJD+8, Cz, derivativesOrder)
    )
    # clenshawMercury <- clenshaw(MJD_TDB, 14, subIntStartMJD, subIntStartMJD+8, cbind(Cx, Cy, Cz), derivatives)
    ## Venus ##
    subInt <- dt %/% 16 + 1
    subIntStartMJD <- rowStartMJD + 16*(subInt - 1)
    idxs <- seq(171 + 30*(subInt-1), by=10, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+9)]
    Cy <- coeffs[idxs[2]:(idxs[2]+9)]
    Cz <- coeffs[idxs[3]:(idxs[3]+9)]
    clenshawVenus <- cbind(
        clenshawAllDerivatives(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+16, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+16, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+16, Cz, derivativesOrder)
    )
    # clenshawVenus <- clenshaw(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+16, cbind(Cx, Cy, Cz), derivatives)
    ## Earth ## (actually Earth-Moon barycenter)
    subInt <- dt %/% 16 + 1 # 16 derived from 32/numberSubInts
    subIntStartMJD <- rowStartMJD + 16*(subInt - 1) 
    idxs <- seq(231 + 39*(subInt-1), by=13, length.out=3) 
    # 231 is the start index of EMB coeffs, 13 is the number of coeffs,
    # 39 coming from 13*3 (3 parameters X Y Z for each coeff)
    Cx <- coeffs[idxs[1]:(idxs[1]+12)]
    Cy <- coeffs[idxs[2]:(idxs[2]+12)]
    Cz <- coeffs[idxs[3]:(idxs[3]+12)]
    clenshawEMB <- cbind(
        clenshawAllDerivatives(MJD_TDB, 13, subIntStartMJD, subIntStartMJD+16, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 13, subIntStartMJD, subIntStartMJD+16, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 13, subIntStartMJD, subIntStartMJD+16, Cz, derivativesOrder)
    )
    # clenshawEMB <- clenshaw(MJD_TDB, 13, subIntStartMJD, subIntStartMJD+16, cbind(Cx, Cy, Cz), derivatives)
    ## Moon ## (already geocentric)
    subInt <- dt %/% 4 + 1
    subIntStartMJD <- rowStartMJD + 4*(subInt - 1)
    idxs <- seq(441 + 39*(subInt-1), by=13, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+12)]
    Cy <- coeffs[idxs[2]:(idxs[2]+12)]
    Cz <- coeffs[idxs[3]:(idxs[3]+12)]
    clenshawMoon <- cbind(
        clenshawAllDerivatives(MJD_TDB, 13, subIntStartMJD, subIntStartMJD+4, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 13, subIntStartMJD, subIntStartMJD+4, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 13, subIntStartMJD, subIntStartMJD+4, Cz, derivativesOrder)
    )
    # clenshawMoon <- clenshaw(MJD_TDB, 13, subIntStartMJD, subIntStartMJD+4, cbind(Cx, Cy, Cz), derivatives)
    ## Mars ## No subintervals for Mars, Jupiter, Saturn, Uranus, Neptune and Pluto
    idxs <- seq(309, by=11, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+10)]
    Cy <- coeffs[idxs[2]:(idxs[2]+10)]
    Cz <- coeffs[idxs[3]:(idxs[3]+10)]
    clenshawMars <- cbind(
        clenshawAllDerivatives(MJD_TDB, 11, rowStartMJD, rowStartMJD+32, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 11, rowStartMJD, rowStartMJD+32, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 11, rowStartMJD, rowStartMJD+32, Cz, derivativesOrder)
    )
    # clenshawMars <- clenshaw(MJD_TDB, 11, rowStartMJD, rowStartMJD+32, cbind(Cx, Cy, Cz), derivatives)
    ## Jupiter ##
    idxs <- seq(342, by=8, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+7)]
    Cy <- coeffs[idxs[2]:(idxs[2]+7)]
    Cz <- coeffs[idxs[3]:(idxs[3]+7)]
    clenshawJupiter <- cbind(
        clenshawAllDerivatives(MJD_TDB, 8, rowStartMJD, rowStartMJD+32, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 8, rowStartMJD, rowStartMJD+32, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 8, rowStartMJD, rowStartMJD+32, Cz, derivativesOrder)
    )
    # clenshawJupiter <- clenshaw(MJD_TDB, 8, rowStartMJD, rowStartMJD+32, cbind(Cx, Cy, Cz), derivatives)
    ## Saturn ##
    idxs <- seq(366, by=7, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+6)]
    Cy <- coeffs[idxs[2]:(idxs[2]+6)]
    Cz <- coeffs[idxs[3]:(idxs[3]+6)]
    clenshawSaturn <- cbind(
        clenshawAllDerivatives(MJD_TDB, 7, rowStartMJD, rowStartMJD+32, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 7, rowStartMJD, rowStartMJD+32, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 7, rowStartMJD, rowStartMJD+32, Cz, derivativesOrder)
    )
    # clenshawSaturn <- clenshaw(MJD_TDB, 7, rowStartMJD, rowStartMJD+32, cbind(Cx, Cy, Cz), derivatives)
    ## Uranus ##
    idxs <- seq(387, by=6, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+5)]
    Cy <- coeffs[idxs[2]:(idxs[2]+5)]
    Cz <- coeffs[idxs[3]:(idxs[3]+5)]
    clenshawUranus <- cbind(
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cz, derivativesOrder)
    )
    # clenshawUranus <- clenshaw(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, cbind(Cx, Cy, Cz), derivatives)
    ## Neptune ##
    idxs <- seq(405, by=6, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+5)]
    Cy <- coeffs[idxs[2]:(idxs[2]+5)]
    Cz <- coeffs[idxs[3]:(idxs[3]+5)]
    clenshawNeptune <- cbind(
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cz, derivativesOrder)
    )
    # clenshawNeptune <- clenshaw(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, cbind(Cx, Cy, Cz), derivatives)
    ## Pluto ##
    idxs <- seq(423, by=6, length.out=3)
    Cx <- coeffs[idxs[1]:(idxs[1]+5)]
    Cy <- coeffs[idxs[2]:(idxs[2]+5)]
    Cz <- coeffs[idxs[3]:(idxs[3]+5)]
    clenshawPluto <- cbind(
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cx, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cy, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, Cz, derivativesOrder)
    )
    # clenshawPluto <- clenshaw(MJD_TDB, 6, rowStartMJD, rowStartMJD+32, cbind(Cx, Cy, Cz), derivatives)
    ## The following 2 sections are commented out because they are not currently used. Actually lunar librations are now used
    # ## Earth Nutations ##
    # subInt <- dt %/% 8 + 1
    # subIntStartMJD <- rowStartMJD + 8*(subInt - 1)
    # idxs <- seq(819 + 20*(subInt-1), by=10, length.out=2)
    # Cdpsi <- coeffs[idxs[1]:(idxs[1]+9)]
    # Cdeps <- coeffs[idxs[2]:(idxs[2]+9)]
    # earthNutation <- 1000 * clenshaw(MJD_TDB, 14, subIntStartMJD, subIntStartMJD+8, cbind(Cdpsi, Cdeps))
    ## Lunar Mantle Librations ##
    subInt <- dt %/% 8 + 1
    subIntStartMJD <- rowStartMJD + 8*(subInt - 1)
    idxs <- seq(899 + 30*(subInt-1), by=10, length.out=3)
    Cphi <- coeffs[idxs[1]:(idxs[1]+9)]
    Ctheta <- coeffs[idxs[2]:(idxs[2]+9)]
    Cpsi <- coeffs[idxs[3]:(idxs[3]+9)]
    # units directly in radians. Order of angles is phi, theta and psi
    # Derivatives in radians/s and radians/s2
    clenshawLunarLibration <- cbind(
        clenshawAllDerivatives(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+8, Cphi, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+8, Ctheta, derivativesOrder),
        clenshawAllDerivatives(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+8, Cpsi, derivativesOrder)
    )
    # clenshawLunarLibration <- clenshaw(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+8, cbind(Cphi, Ctheta, Cpsi), derivatives)
    # lunarLibration <- clenshaw(MJD_TDB, 10, subIntStartMJD, subIntStartMJD+8, cbind(Cphi, Ctheta, Cpsi), derivatives)
    positionSun <- clenshawSun[1, ] * 1000
    positionMercury <- clenshawMercury[1, ] * 1000
    positionVenus <- clenshawVenus[1, ] * 1000
    positionEMB <- clenshawEMB[1, ] * 1000
    positionMoon <- clenshawMoon[1, ] * 1000
    positionMars <- clenshawMars[1, ] * 1000
    positionJupiter <- clenshawJupiter[1, ] * 1000
    positionSaturn <- clenshawSaturn[1, ] * 1000
    positionUranus <- clenshawUranus[1, ] * 1000
    positionNeptune <- clenshawNeptune[1, ] * 1000
    positionPluto <- clenshawPluto[1, ] * 1000
    lunarLibration <- clenshawLunarLibration[1, ]
    if(derivativesOrder >= 1) {
        velocitySun <- clenshawSun[2, ] * 1000
        velocityMercury <- clenshawMercury[2, ] * 1000
        velocityVenus <- clenshawVenus[2, ] * 1000
        velocityEMB <- clenshawEMB[2, ] * 1000
        velocityMoon <- clenshawMoon[2, ] * 1000
        velocityMars <- clenshawMars[2, ] * 1000
        velocityJupiter <- clenshawJupiter[2, ] * 1000
        velocitySaturn <- clenshawSaturn[2, ] * 1000
        velocityUranus <- clenshawUranus[2, ] * 1000
        velocityNeptune <- clenshawNeptune[2, ] * 1000
        velocityPluto <- clenshawPluto[2, ] * 1000
        derivativesLunarLibration <- clenshawLunarLibration[2, ]
        if(derivativesOrder >= 2) {
            accelerationSun <- clenshawSun[3, ] * 1000
            accelerationMercury <- clenshawMercury[3, ] * 1000
            accelerationVenus <- clenshawVenus[3, ] * 1000
            accelerationEMB <- clenshawEMB[3, ] * 1000
            accelerationMoon <- clenshawMoon[3, ] * 1000
            accelerationMars <- clenshawMars[3, ] * 1000
            accelerationJupiter <- clenshawJupiter[3, ] * 1000
            accelerationSaturn <- clenshawSaturn[3, ] * 1000
            accelerationUranus <- clenshawUranus[3, ] * 1000
            accelerationNeptune <- clenshawNeptune[3, ] * 1000
            accelerationPluto <- clenshawPluto[3, ] * 1000
            secondDerivativesLunarLibration <- clenshawLunarLibration[3, ]
        }
    }
    positionEarth <- positionEMB - positionMoon*EMRAT1_DE440
    centralBodySSBPosition <- switch(EXPR = centralBody,
                                     SSB = 0,
                                     Mercury = positionMercury,
                                     Venus = positionVenus,
                                     Earth = positionEarth,
                                     Moon = positionMoon + positionEarth,
                                     Mars = positionMars,
                                     Jupiter = positionJupiter,
                                     Saturn = positionSaturn,
                                     Uranus = positionUranus,
                                     Neptune = positionNeptune,
                                     Pluto = positionPluto,
                                     positionEarth)
    positionEarth <- positionEarth - centralBodySSBPosition
    positionMercury <- positionMercury - centralBodySSBPosition
    positionVenus <- positionVenus - centralBodySSBPosition
    positionMars <- positionMars - centralBodySSBPosition
    positionJupiter <- positionJupiter - centralBodySSBPosition
    positionSaturn <- positionSaturn - centralBodySSBPosition
    positionUranus <- positionUranus - centralBodySSBPosition
    positionNeptune <- positionNeptune - centralBodySSBPosition
    positionPluto <- positionPluto - centralBodySSBPosition
    positionSunBodycentric <- positionSun - centralBodySSBPosition
    positionMoon <- positionMoon + positionEarth
    if(derivativesOrder >= 1) {
        velocityEarth <- velocityEMB - velocityMoon*EMRAT1_DE440
        centralBodySSBVelocity <- switch(EXPR = centralBody,
                                         SSB = 0,
                                         Mercury = velocityMercury,
                                         Venus = velocityVenus,
                                         Earth = velocityEarth,
                                         Moon = velocityMoon + velocityEarth,
                                         Mars = velocityMars,
                                         Jupiter = velocityJupiter,
                                         Saturn = velocitySaturn,
                                         Uranus = velocityUranus,
                                         Neptune = velocityNeptune,
                                         Pluto = velocityPluto,
                                         velocityEarth)
        velocityEarth <- velocityEarth - centralBodySSBVelocity
        velocityMercury <- velocityMercury - centralBodySSBVelocity
        velocityVenus <- velocityVenus - centralBodySSBVelocity
        velocityMars <- velocityMars - centralBodySSBVelocity
        velocityJupiter <- velocityJupiter - centralBodySSBVelocity
        velocitySaturn <- velocitySaturn - centralBodySSBVelocity
        velocityUranus <- velocityUranus - centralBodySSBVelocity
        velocityNeptune <- velocityNeptune - centralBodySSBVelocity
        velocityPluto <- velocityPluto - centralBodySSBVelocity
        velocitySunBodycentric <- velocitySun - centralBodySSBVelocity
        velocityMoon <- velocityMoon + velocityEarth
        if(derivativesOrder >= 2) {
            accelerationEarth <- accelerationEMB - accelerationMoon*EMRAT1_DE440
            centralBodySSBacceleration <- switch(EXPR = centralBody,
                                             SSB = 0,
                                             Mercury = accelerationMercury,
                                             Venus = accelerationVenus,
                                             Earth = accelerationEarth,
                                             Moon = accelerationMoon + accelerationEarth,
                                             Mars = accelerationMars,
                                             Jupiter = accelerationJupiter,
                                             Saturn = accelerationSaturn,
                                             Uranus = accelerationUranus,
                                             Neptune = accelerationNeptune,
                                             Pluto = accelerationPluto,
                                             accelerationEarth)
            accelerationEarth <- accelerationEarth - centralBodySSBacceleration
            accelerationMercury <- accelerationMercury - centralBodySSBacceleration
            accelerationVenus <- accelerationVenus - centralBodySSBacceleration
            accelerationMars <- accelerationMars - centralBodySSBacceleration
            accelerationJupiter <- accelerationJupiter - centralBodySSBacceleration
            accelerationSaturn <- accelerationSaturn - centralBodySSBacceleration
            accelerationUranus <- accelerationUranus - centralBodySSBacceleration
            accelerationNeptune <- accelerationNeptune - centralBodySSBacceleration
            accelerationPluto <- accelerationPluto - centralBodySSBacceleration
            accelerationSunBodycentric <- accelerationSun - centralBodySSBacceleration
            accelerationMoon <- accelerationMoon + accelerationEarth
        }
    }
    output <- list(
        positionSun = positionSunBodycentric,
        positionSunSSBarycentric = positionSun,
        positionMercury = positionMercury,
        positionVenus = positionVenus,
        positionEarth = positionEarth,
        positionMars = positionMars,
        positionJupiter = positionJupiter,
        positionSaturn = positionSaturn,
        positionUranus = positionUranus,
        positionNeptune = positionNeptune,
        positionPluto = positionPluto,
        positionMoon = positionMoon,
        lunarLibrationAngles = lunarLibration
    )
    if(derivativesOrder >= 1) {
        firstDerivativesOutput <- list(
            velocitySun = velocitySunBodycentric,
            velocitySunSSBarycentric = velocitySun,
            velocityMercury = velocityMercury,
            velocityVenus = velocityVenus,
            velocityEarth = velocityEarth,
            velocityMars = velocityMars,
            velocityJupiter = velocityJupiter,
            velocitySaturn = velocitySaturn,
            velocityUranus = velocityUranus,
            velocityNeptune = velocityNeptune,
            velocityPluto = velocityPluto,
            velocityMoon = velocityMoon,
            lunarLibrationAnglesDerivatives = derivativesLunarLibration
        )
        output <- c(output, firstDerivativesOutput)
        if(derivativesOrder >= 2) {
            secondDerivativesOutput <- list(
                accelerationSun = accelerationSunBodycentric,
                accelerationSunSSBarycentric = accelerationSun,
                accelerationMercury = accelerationMercury,
                accelerationVenus = accelerationVenus,
                accelerationEarth = accelerationEarth,
                accelerationMars = accelerationMars,
                accelerationJupiter = accelerationJupiter,
                accelerationSaturn = accelerationSaturn,
                accelerationUranus = accelerationUranus,
                accelerationNeptune = accelerationNeptune,
                accelerationPluto = accelerationPluto,
                accelerationMoon = accelerationMoon,
                lunarLibrationAnglesSecondDerivatives = secondDerivativesLunarLibration
            )
            output <- c(output, secondDerivativesOutput)
        }
    }
    return(output)
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

# 
# .legendre <- function(n, m, phi) {
#     # Lear normalized algorithm
#     pnm <- matrix(0, nrow=n+1, ncol=m+1)
#     dpnm <- matrix(0, nrow=n+1, ncol=m+1)
#     pnm[1,1] <- 1
#     dpnm[1,1] <- 0
#     pnm[2,1] <- sqrt(3) * sin(phi)
#     dpnm[2,1] <- sqrt(3) * cos(phi)
#     pnm[2,2] <- sqrt(3) * cos(phi)
#     # normalization factors are named N1, N2, etc
#     # these refer to the lambda normalization parameters described in equations
#     # 3.1 to 3.5 of NASA document
#     # Zonal terms, m=0
#     for(i in 2:n) {
#         N1 <- sqrt((2*n+1) / (2*n-1)) 
#         N2 <- sqrt((2*i + 1) / (2*i-3))
#         pnm[i+1, 1] <- (N1 * (2*i-1) * sin(phi) * pnm[i, 1] - 
#                             N2 * (n-1) * pnm[i-1, 1])/i
#         dpnm[i+1, 1] <- N1 * (sin(phi) * dpnm[i, 1] + i * pnm[i, 1])
#     }
#     # Sectorial terms, m=n
#     for(i in 2:n) {
#         N3 <- sqrt((2*i+1) / (2*i))
#         # original parameter 3.3 is divided by 2n-1, but after multiplied by
#         # same value, so we can remove it
#         pnm[i+1, i+1] <- N3 * cos(phi) * pnm[i, i]
#     } TODO
# }

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
    return(as.numeric(lambda))
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
        rs <- sunRadius/rb # apparent radius of Sun from satellite
        rp <- earthRadius_EGM96/rsbx # apparent radius of Earth from satellite
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
            rs <- sunRadius/rb
            rp <- moonRadius/rsbx
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

determineCentralBody <- function(position, planetEphemerides, moonEphemeris) {
    distances <- sapply(planetEphemerides, FUN=function(x) {sqrt(sum((x - position)^2))})
    centralBody <- substring(names(which(distances < SOIradii)), 9)
    if(length(centralBody) == 0) {
        centralBody <- "SSB"
    }
    if(centralBody=="Earth") {
        if(sqrt(sum((moonEphemeris - position)^2)) < SOIMoon) {
            centralBody <- "Moon"
        }
    }
    return(centralBody)
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

# 
# .legendre <- function(n, m, phi) {
#     # Lear normalized algorithm
#     pnm <- matrix(0, nrow=n+1, ncol=m+1)
#     dpnm <- matrix(0, nrow=n+1, ncol=m+1)
#     pnm[1,1] <- 1
#     dpnm[1,1] <- 0
#     pnm[2,1] <- sqrt(3) * sin(phi)
#     dpnm[2,1] <- sqrt(3) * cos(phi)
#     pnm[2,2] <- sqrt(3) * cos(phi)
#     # normalization factors are named N1, N2, etc
#     # these refer to the lambda normalization parameters described in equations
#     # 3.1 to 3.5 of NASA document
#     # Zonal terms, m=0
#     for(i in 2:n) {
#         N1 <- sqrt((2*n+1) / (2*n-1)) 
#         N2 <- sqrt((2*i + 1) / (2*i-3))
#         pnm[i+1, 1] <- (N1 * (2*i-1) * sin(phi) * pnm[i, 1] - 
#                             N2 * (n-1) * pnm[i-1, 1])/i
#         dpnm[i+1, 1] <- N1 * (sin(phi) * dpnm[i, 1] + i * pnm[i, 1])
#     }
#     # Sectorial terms, m=n
#     for(i in 2:n) {
#         N3 <- sqrt((2*i+1) / (2*i))
#         # original parameter 3.3 is divided by 2n-1, but after multiplied by
#         # same value, so we can remove it
#         pnm[i+1, i+1] <- N3 * cos(phi) * pnm[i, i]
#     } TODO
# }

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
    return(as.numeric(lambda))
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
        rs <- sunRadius/rb # apparent radius of Sun from satellite
        rp <- earthRadius_EGM96/rsbx # apparent radius of Earth from satellite
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
            rs <- sunRadius/rb
            rp <- moonRadius/rsbx
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