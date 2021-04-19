iauCal2jd <- function(year, month, day) {
    djm0 <- 2400000.5
    b <- 0
    c <- 0
    if (month <= 2) {
        year <- year -1
        month <- month + 12
    }
    if (year < 0) {
        c <- -0.75
    }
    if (year < 1582) {
        # if statement only for subsequent conditionals. TODO CHANGE
    } else if (year > 1582) {
        a <- trunc(year/100)
        b <- 2 - a + floor(a / 4)
    } else if (month < 10) {
        # TODO CHANGE.
    } else if (month > 10) { ## TODO VERIFY introduction of gregorian in 1582
        a <- trunc(year/100)
        b <- 2 - a + floor(a / 4)
    } else if (day <= 4) {
        # TODO CHANGE
    } else if (day > 14) {
        a <- trunc(year/100)
        b <- 2 - a + floor(a / 4)
    } else {
        stop("Please enter a valid calendar date")
    }
    jd <- trunc(365.25 * year + c) + trunc(30.6001 * (month + 1))
    jd <- jd + day + b + 1720994.5
    djm <- jd-djm0
    return(list(
        DJMJD0=djm0,
        DATE=djm
    ))
}

iauObl06 <- function(date1, date2) {
    # Interval between fundamental date J2000.0 and given date (JC).
    t <- ((date1 - DJ00) + date2) / DJC
    # Mean obliquity
    eps0 <- (84381.406+(-46.836769+(-0.0001831+(0.00200340+(-0.000000576-0.0000000434*t)*t)*t)*t)*t)*DAS2R
    return(eps0)
}

iauPfw06 <- function(date1, date2) {
    # Interval between fundamental date J2000.0 and given date (JC).
    t <- ((date1 - DJ00) + date2) / DJC
    # P03 bias+precession angles
    gamb <- (-0.052928+(10.556378+(0.4932044+(-0.00031238+(-0.000002788+0.0000000260*t)*t)*t)*t)*t)*DAS2R
    phib <- (84381.412819+(-46.811016+(0.0511268+(0.00053289+(-0.000000440-0.0000000176*t)*t)*t)*t)*t)*DAS2R
    psib <- (-0.041775+(5038.481484+(1.5584175+(-0.00018522+(-0.000026452-0.0000000148*t)*t)*t)*t)*t)*DAS2R
    epsa <- iauObl06(date1, date2)
    return(list(
        gamb=gamb, 
        phib=phib, 
        psib=psib, 
        epsa=epsa
    ))
}

iauFal03 <- function(t) {
    #  t is TDB, Julian centuries since J2000.0 (Note 1)
    # Mean anomaly of the Moon (IERS Conventions 2003).
    a <- ((485868.249036+t*(1717915923.2178+t*(31.8792+t*(0.051635+t*(-0.00024470)))))%%TURNAS)*DAS2R
    return(a)
}

iauFaf03 <- function(t) {
    # Mean longitude of the Moon minus that of the ascending node
    # (IERS Conventions 2003)
    a <- ((335779.526232+t*(1739527262.8478+t*(-12.7512+t*(-0.001037+t*(0.00000417)))))%%TURNAS)*DAS2R
    return(a)
}

iauFaom03 <- function(t) {
    # Mean longitude of the Moon's ascending node (IERS Conventions 2003).
    a <- ((450160.398036+t*(-6962890.5431+t*(7.4722+t*(0.007702+t*(-0.00005939)))))%%TURNAS)*DAS2R
    return(a)
}

iauFalp03 <- function(t) {
    a <- ((1287104.79305+t*(129596581.0481+t*(-0.5532+t*(0.000136+t*(-0.00001149)))))%%TURNAS)*DAS2R
    return(a)
}

iauFad03 <- function(t) {
    a <- ((1072260.70369+t*(1602961601.2090+t*(-6.3706+t*(0.006593+t*(-0.00003169)))))%%TURNAS)*DAS2R
    return(a)
}

iauFave03 <- function(t) {
    a <- (3.176146697 + 1021.3285546211 * t) %% (2*pi)
    return(a)
}

iauFae03 <- function(t) {
    a <- (1.753470314 + 628.3075849991 * t) %% (2*pi)
    return(a)
}

iauFapa03 <- function(t) {
    a <- (0.024381750 + 0.00000538691 * t) * t
    return(a)
}

iauNut00a <- function(date1, date2){
    # Units of 0.1 microarcsecond to radians
    U2R <- DAS2R / 1e7
    # -------------------------
        # Luni-Solar nutation model
    # -------------------------
        
        # The units for the sine and cosine coefficients are
    # 0.1 microarcsecond and the same per Julian century
    
    # nl,nlp,nf,nd,nom:   coefficients of l,l',F,D,Om
    #  sp,spt,cp:          longitude sin, t*sin, cos coefficients
    #  ce,cet,se:          obliquity cos, t*cos, sin coefficients
    #  moved the following 4 lines to constants definitions
    # xls <- read.csv("R/hpop_files/iauNut00a_xls.csv", header=FALSE)
    # NLS <- nrow(xls)
    # xpl <- read.csv("R/hpop_files/iauNut00a_xpl.csv", header=FALSE)
    # NPL <- nrow(xpl)
    # Interval between fundamental date J2000.0 and given date (JC). 
    t <- ((date1 - DJ00) + date2) / DJC
    # Mean anomaly of the Moon (IERS 2003)
    el <- iauFal03(t)
    # Mean anomaly of the Sun (MHB2000) AQUIAQUIAQUI
    elp <- iauFalp03(t) 
    # Mean longitude of the Moon minus that of the ascending node (IERS 2003).
    f <- iauFaf03(t)
    # Mean elongation of the Moon from the Sun (MHB2000). 
    d <- iauFad03(t)
    # Mean longitude of the ascending node of the Moon (IERS 2003).
    om <- iauFaom03(t)
    #  Summation of luni-solar nutation series (in reverse order).
    #  replaced by vectorized code
    arg <- (xls[NLS:1,1]*el+xls[NLS:1,2]*elp+xls[NLS:1,3]*f+xls[NLS:1,4]*d+xls[NLS:1,5]*om)%%(2*pi)
    sarg <- sin(arg)
    carg <- cos(arg)
    dp <- sum((xls[NLS:1, 6] + xls[NLS:1, 7] * t) * sarg + xls[NLS:1, 8] * carg)
    de <- sum((xls[NLS:1, 9] + xls[NLS:1, 10] * t) * carg + xls[NLS:1, 11] * sarg)
    # Old non vectorized code
    # for (i in NLS:1) {
    #     # argument and its sine and cosine
    #     arg <- (xls[i,1]*el+xls[i,2]*elp+xls[i,3]*f+xls[i,4]*d+xls[i,5]*om)%%(2*pi)
    #     sarg <- sin(arg)
    #     carg <- cos(arg)
    #     # term
    #     dp <- dp + (xls[i,6] + xls[i,7] * t) * sarg + xls[i,8] * carg
    #     de <- de + (xls[i,9] + xls[i,10] * t) * carg + xls[i,11] * sarg
    # }
    # Convert from decimes microarcseconds to radians
    dpsils <- dp * U2R
    depsls <- de * U2R
    # PLANETARY NUTATION 
    #  n.b.  The MHB2000 code computes the luni-solar and planetary nutation 
    #  in different functions, using slightly different Delaunay 
    #  arguments in the two cases.  This behaviour is faithfully 
    #  reproduced here.  Use of the IERS 2003 expressions for both 
    #  cases leads to negligible changes, well below 
    #  0.1 microarcsecond. 
    # Mean anomaly of the Moon (MHB2000). 
    al <- (2.35555598 + 8328.6914269554 * t) %% (2*pi)
    # Mean longitude of the Moon minus that of the ascending node (MHB2000). 
    af <- (1.627905234 + 8433.466158131 * t) %% (2*pi)
    # Mean elongation of the Moon from the Sun (MHB2000)
    ad <- (5.198466741 + 7771.3771468121 * t) %% (2*pi)
    # Mean longitude of the ascending node of the Moon (MHB2000)
    aom <- (2.18243920 - 33.757045 * t) %% (2*pi)
    # General accumulated precession in longitude (IERS 2003)
    apa <- iauFapa03(t)
    # Mean longitude of Mercury (IERS Conventions 2003)
    alme <- (4.402608842 + 2608.7903141574 * t) %% (2*pi)
    # Mean longitude of Venus (IERS Conventions 2003)
    alve <- iauFave03(t)
    # Mean longitude of Earth (IERS Conventions 2003)
    alea <- iauFae03(t)
    # Mean longitude of Mars (IERS Conventions 2003)
    alma <- (6.203480913 + 334.0612426700 * t) %% (2*pi)
    # Mean longitude of Jupiter (IERS Conventions 2003)
    alju <- (0.599546497 + 52.9690962641 * t) %% (2*pi)
    # Mean longitude of Saturn (IERS Conventions 2003)
    alsa <- (0.874016757 + 21.3299104960 * t) %% (2*pi)
    # Mean longitude of Uranus (IERS Conventions 2003)
    alur <- (5.481293872 + 7.4781598567 * t) %% (2*pi)
    # Neptune longitude (MHB2000). 
    alne <- (5.321159000 + 3.8127774000 * t) %% (2*pi)
    ## Calculate again nutation values (temporary variables)
    arg <- (xpl[NPL:1, 1] * al + xpl[NPL:1, 2] * af + xpl[NPL:1, 3] * ad + xpl[NPL:1, 4] * aom +
                xpl[NPL:1, 5] * alme + xpl[NPL:1, 6] * alve + xpl[NPL:1, 7] * alea +
                xpl[NPL:1, 8] * alma + xpl[NPL:1, 9] * alju + xpl[NPL:1, 10] * alsa + 
                xpl[NPL:1, 11] * alur + xpl[NPL:1, 12] * alne + xpl[NPL:1, 13] * apa) %% (2*pi)
    sarg <- sin(arg)
    carg <- cos(arg)
    dp <- sum(xpl[NPL:1, 14] * sarg + xpl[NPL:1, 15] * carg)
    de <- sum(xpl[NPL:1, 16] * sarg + xpl[NPL:1, 17] * carg)
    # old non-vectorized code
    # for (i in NPL:1) {
    #     # Argument and its sine and cosine
    #     arg <- (xpl[i,1] * al + xpl[i,2] * af + xpl[i,3] * ad + xpl[i,4] * aom +
    #             xpl[i,5] * alme + xpl[i,6] * alve +xpl[i,7] * alea +
    #             xpl[i,8] * alma + xpl[i,9] * alju +xpl[i,10] * alsa + 
    #             xpl[i,11] * alur + xpl[i,12] * alne +xpl[i,13] * apa) %% (2*pi)
    #     sarg <- sin(arg)
    #     carg <- cos(arg)
    #     # Update terms
    #     dp <- dp + xpl[i,14] * sarg + xpl[i,15] * carg
    #     de <- de + xpl[i,16] * sarg + xpl[i,17] * carg
    # }
    # Convert from 0.1 microarcsecs  to radians
    dpsipl <- dp * U2R
    depspl <- de * U2R
    # Add lunar-solar and planetary componentes of nutation
    dpsi <- dpsils + dpsipl
    deps <- depsls + depspl
    return(list(
        dpsi=dpsi,
        deps=deps
    ))
}

iauNut06a <- function(date1, date2){
    # Interval between fundamental date J2000.0 and given date (JC).
    t <- ((date1 - DJ00) + date2) / DJC
    # Factor correcting for secular variation of J2.
    fj2 = -2.7774e-6 * t
    # Obtain IAU 2000A nutation.
    iau2000aNutation <- iauNut00a(date1, date2)
    # Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5).
    dpsi <- iau2000aNutation$dpsi + iau2000aNutation$dpsi * (0.4697e-6 + fj2)
    deps <- iau2000aNutation$deps + iau2000aNutation$deps * fj2
    return(list(
        dpsi=dpsi,
        deps=deps
    ))
}

iauRz <- function(psi, r) {
    sin_psi <- sin(psi)
    cos_psi <- cos(psi)
    r <- matrix(c(cos_psi*r[1,1] + sin_psi*r[2,1], # start new row 1
                  cos_psi*r[1,2] + sin_psi*r[2,2], 
                  cos_psi*r[1,3] + sin_psi*r[2,3],
                  cos_psi*r[2,1] - sin_psi*r[1,1], # start new row 2
                  cos_psi*r[2,2] - sin_psi*r[1,2], 
                  cos_psi*r[2,3] - sin_psi*r[1,3],
                  r[3,1], r[3,2], r[3,3]), 
                byrow=TRUE, nrow=3)
    return(r)
}

iauRx <- function(phi, r) {
    sin_phi <- sin(phi)
    cos_phi <- cos(phi)
    r <- matrix(c(r[1,1], r[1,2], r[1,3],
                  cos_phi*r[2,1] + sin_phi*r[3,1], # start new row 2
                  cos_phi*r[2,2] + sin_phi*r[3,2], 
                  cos_phi*r[2,3] + sin_phi*r[3,3],
                  cos_phi*r[3,1] - sin_phi*r[2,1], # start new row 3
                  cos_phi*r[3,2] - sin_phi*r[2,2], 
                  cos_phi*r[3,3] - sin_phi*r[2,3]),
                byrow=TRUE, nrow=3)
    return(r)
}

iauRy <- function(theta, r) {
    sin_theta <- sin(theta)
    cos_theta <- cos(theta)
    r <- matrix(c(cos_theta*r[1,1] - sin_theta*r[3,1], # start new row 1
                  cos_theta*r[1,2] - sin_theta*r[3,2], 
                  cos_theta*r[1,3] - sin_theta*r[3,3],
                  r[2,1], r[2,2], r[2,3],
                  cos_theta*r[3,1] + sin_theta*r[1,1], # start new row 3
                  cos_theta*r[3,2] + sin_theta*r[1,2], 
                  cos_theta*r[3,3] + sin_theta*r[1,3]),
                byrow=TRUE, nrow=3)
    return(r)
}



iauFw2m <- function(gamb, phib, psi, eps) {
    r <- diag(3)
    r <- iauRz(gamb, r)
    r <- iauRx(phib, r)
    r <- iauRz(-psi, r)
    r <- iauRx(-eps, r)
    return(r)
}

iauPnm06a <- function(date1, date2) {
    # Fukushima-Williams angles for frame bias and precession
    iauPfw06_result <- iauPfw06(date1, date2)
    # Nutation components
    iauNut06a_result <- iauNut06a(date1, date2)
    # Equinox based nutation x precession x bias matrix
    iauFw2m_result <- iauFw2m(iauPfw06_result$gamb, iauPfw06_result$phib, 
                       iauPfw06_result$psib + iauNut06a_result$dpsi, 
                       iauPfw06_result$epsa + iauNut06a_result$deps)
    return(iauFw2m_result)
}

iauS06 <- function(date1, date2, x, y) {
    # All of the following have been moved to constants definitions
    # sp <- c(94.00e-6, 3808.65e-6, -122.68e-6, -72574.11e-6, 27.98e-6, 15.62e-6)
    # s_xyD2_coefs <- read.csv("R/hpop_files/s_xyD2_terms.csv", header=FALSE)
    # s0 <- s_xyD2_coefs[s_xyD2_coefs[, 1] == "ORDER0", -1]
    # s1 <- s_xyD2_coefs[s_xyD2_coefs[, 1] == "ORDER1", -1]
    # s2 <- s_xyD2_coefs[s_xyD2_coefs[, 1] == "ORDER2", -1]
    # s3 <- s_xyD2_coefs[s_xyD2_coefs[, 1] == "ORDER3", -1]
    # s4 <- s_xyD2_coefs[s_xyD2_coefs[, 1] == "ORDER4", -1]
    # w0 <- sp[1]
    # w1 <- sp[2]
    # w2 <- sp[3]
    # w3 <- sp[4]
    # w4 <- sp[5]
    # w5 <- sp[6]
    t <- ((date1 - DJ00) + date2) / DJC
    meanAnomalyMoon <- iauFal03(t)
    meanAnomalySun <- iauFalp03(t)
    meanLongitudeMoonMinusAN <- iauFaf03(t)
    meanElongationMoonSun <- iauFad03(t)
    meanLongitudeMoonAN <- iauFaom03(t)
    meanLongitudeVenus <- iauFave03(t)
    meanLongitudeEarth <- iauFae03(t)
    generalLongitudeAccumulatedPrecesion <- iauFapa03(t)
    fundamentalArguments <- c(meanAnomalyMoon, meanAnomalySun, meanLongitudeMoonMinusAN,
                              meanElongationMoonSun, meanLongitudeMoonAN, meanLongitudeVenus,
                              meanLongitudeEarth, generalLongitudeAccumulatedPrecesion)
    for (i in nrow(s0):1) {
        # a <- 0
        a <- sum(s0[i, 1:8] * fundamentalArguments)
        # for (j in 1:8) {
        #     a <- a + s0[i, j] * fundamentalArguments[j]
        # }
        w0 <-  w0 + s0[i, 9] * sin(a) + s0[i, 10] * cos(a)
    }
    for (i in nrow(s1):1) {
        a <- 0
        for (j in 1:8) {
            a <- a + s1[i, j] * fundamentalArguments[j]
        }
        w1 <-  w1 + s1[i, 9] * sin(a) + s1[i, 10] * cos(a)
    }
    for (i in nrow(s2):1) {
        a <- 0
        for (j in 1:8) {
            a <- a + s2[i, j] * fundamentalArguments[j]
        }
        w2 <-  w2 + s2[i, 9] * sin(a) + s2[i, 10] * cos(a)
    }
    for (i in nrow(s3):1) {
        a <- 0
        for (j in 1:8) {
            a <- a + s3[i, j] * fundamentalArguments[j]
        }
        w3 <-  w3 + s3[i, 9] * sin(a) + s3[i, 10] * cos(a)
    }
    for (i in nrow(s4):1) {
        a <- 0
        for (j in 1:8) {
            a <- a + s4[i, j] * fundamentalArguments[j]
        }
        w4 <- w4 + s4[i, 9] * sin(a) + s4[i, 10] * cos(a)
    }
    s = (w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t) * DAS2R - x * y / 2.0
    return(s)
}

iauEra00 <- function(dj1, dj2) {
    if (dj1 < dj2) {
        d1 <- dj1
        d2 <- dj2
    } else {
        d1 <- dj2
        d2 <- dj1
    }
    t <- d1 + (d2 - DJ00)
    f <- (d1 - trunc(d1)) + (d2 - trunc(d2))
    theta <- rem((f + 0.7790572732640 + 0.00273781191135448 * t) * 2*pi, (2*pi))
    if (theta < 0) {
        theta <- theta + 2*pi
    }
    return(theta)
}

iauEors <- function(rnpb, s) {
    x <- rnpb[3,1]
    ax <-  x / (1 + rnpb[3,3])
    xs <- 1 - ax * x
    ys <- -ax * rnpb[3,2]
    zs <- -x
    p <- rnpb[1,1] * xs + rnpb[1,2] * ys + rnpb[1,3] * zs
    q <- rnpb[2,1] * xs + rnpb[2,2] * ys + rnpb[2,3] * zs
    eo <- if (p != 0 | q != 0) (s - atan2(q, p)) else s
    return(eo)
}

iauGst06 <- function(uta, utb, tta, ttb, rnpb) {
    cip_x <- rnpb[3, 1]
    cip_y <- rnpb[3, 2]
    s <- iauS06(tta, ttb, cip_x, cip_y)
    era <- iauEra00(uta, utb)
    eors <- iauEors(rnpb, s)
    gst <- rem(era - eors, 2*pi)
    if (gst < 0) {
        gst <- gst + 2*pi
    }
    return(gst)
}

iauPom00 <- function(xp, yp, sp) {
    rpom <- iauRx(-yp,
                  iauRy(-xp,
                        iauRz(sp, diag(3))))
    return(rpom)
}

iauSp00 <- function(date1, date2){
    t <- ((date1 - DJ00) + date2) / DJC
    sp <- -47e-6 * t * DAS2R
    return(sp)
}

iauGmst06 <- function(uta, utb, tta, ttb) {
    # TT Julian centuries since J2000.0
    t <- ((tta - DJ00) + ttb) / DJC
    gmst <-
        iauEra00(uta, utb) + (0.014506 + (4612.156534 + (1.3915817 + (
            -0.00000044 + (-0.000029956 + (-0.0000000368) * t) * t
        ) * t) * t) * t) * DAS2R
    gmst <- rem(gmst, 2*pi)
    if (gmst < 0) {
        gmst <- gmst + 2*pi
    }
    return(gmst)
}
