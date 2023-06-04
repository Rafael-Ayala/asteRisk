anelasticEarthAcceleration <- function(Mjd_UTC, r_sun, r_moon, r, E, UT1_UTC,
                                     TT_UTC, x_pole, y_pole, earthSPHDegree, SETcorrections,
                                     OTcorrections) {
    if((earthSPHDegree+1) > nrow(asteRiskData::GGM05C_Cnm)) {
        warning(strwrap(paste("Spherical harmonics coefficients for Earth are only available
                            up to degree ", nrow(asteRiskData::GGM05C_Cnm)-1, ", but ",
                              earthSPHDegree, " was input. Defaulting to the maximum
                            available degree and order of ", nrow(asteRiskData::GGM05C_Cnm)-1, ".")))
        earthSPHDegree <- nrow(asteRiskData::GGM05C_Cnm)-1
    }
    C <- asteRiskData::GGM05C_Cnm[1:(earthSPHDegree+1), 1:(earthSPHDegree+1)]
    S <- asteRiskData::GGM05C_Snm[1:(earthSPHDegree+1), 1:(earthSPHDegree+1)]
    r_moon <- E %*% r_moon
    moon_polar <- CartesianToPolar(r_moon)
    r_sun <- E %*% r_sun
    sun_polar <- CartesianToPolar(r_sun)
    Mjd_TT <- Mjd_UTC + TT_UTC/86400
    Mjd_UT1 <- Mjd_UTC + UT1_UTC/86400
    Time <- (Mjd_TT-MJD_J2000)/36525
    Time_years <- (Mjd_TT-MJD_J2000)/365.25
    if(earthSPHDegree >= 4) {
        ## Secular variations of low order coefficients
        C[3,1] <- C[3,1] + 11.6e-12*Time_years
        C[4,1] <- C[4,1] + 4.9e-12*Time_years
        C[5,1] <- C[5,1] + 4.7e-12*Time_years
        if(SETcorrections) {
            dCnmSET <- matrix(0, nrow=5, ncol=5)
            dSnmSET <- matrix(0, nrow=5, ncol=5)
            # First, let's apply the secular long-term effects
            # Secular variation of low degree coefficients (section 6.1 of IERS note, eq. 6.4 and Table 6.2)
            dCnmSET[3, 1] <- 11.6e-12*(Time_years)
            dCnmSET[4, 1] <- 4.9e-12*(Time_years)
            dCnmSET[5, 1] <- 4.7e-12*(Time_years)
            # Start corrections for solid Earth tides
            legendre_moonTheta <- legendre(4, 4, moon_polar[2])[[1]]
            legendre_sunTheta <- legendre(4, 4, sun_polar[2])[[1]]
            # step 1 of corrections for solid Earth tides
            # matrixes for corrections due to solid Earth tides (SET)
            # not all coefficients (n 0,1) are corrected, but using 5x5 matrix because
            # up to n=4 is corrected. note not all m for each n are corrected either
            # (for n=4, only m 0, 1 and 2). See section 6.2 of IERS notes on Geopotential
            for(n in c(2,3)) {
                correctionFactorMoon <- GMratioMoonEarth*((earthRadius_EGM96/moon_polar[3])^(n+1))
                correctionFactorSun <- GMratioSunEarth*((earthRadius_EGM96/sun_polar[3])^(n+1))
                for(m in 0:n) {
                    correctionFactorLoveRe <- reKnm0An[n-1, m+1]/(2*n+1)
                    correctionFactorLoveIm <- imKnm0An[n-1, m+1]/(2*n+1)
                    dCnm <- correctionFactorLoveRe * (correctionFactorMoon*legendre_moonTheta[n+1, m+1]*cos(m*moon_polar[1]) + 
                                                          correctionFactorSun*legendre_sunTheta[n+1, m+1]*cos(m*sun_polar[1])) +
                        correctionFactorLoveIm * (correctionFactorMoon*legendre_moonTheta[n+1, m+1]*sin(m*moon_polar[1]) + 
                                                      correctionFactorSun*legendre_sunTheta[n+1, m+1]*sin(m*sun_polar[1]))
                    dSnm <- correctionFactorLoveRe * (correctionFactorMoon*legendre_moonTheta[n+1, m+1]*sin(m*moon_polar[1]) + 
                                                          correctionFactorSun*legendre_sunTheta[n+1, m+1]*sin(m*sun_polar[1])) -
                        correctionFactorLoveIm * (correctionFactorMoon*legendre_moonTheta[n+1, m+1]*cos(m*moon_polar[1]) + 
                                                      correctionFactorSun*legendre_sunTheta[n+1, m+1]*cos(m*sun_polar[1]))
                    dCnmSET[n+1, m+1] <- dCnm
                    dSnmSET[n+1, m+1] <- dSnm
                }
            }
            n <- 4
            correctionFactorMoon <- GMratioMoonEarth*((earthRadius_EGM96/moon_polar[3])^3)
            correctionFactorSun <- GMratioSunEarth*((earthRadius_EGM96/sun_polar[3])^3)
            for(m in 0:2) {
                correctionFactorLove <- knmplusAn[m+1]/(5)
                dCnm <- correctionFactorLove * (correctionFactorMoon*legendre_moonTheta[n+1, m+1]*cos(m*moon_polar[1]) + 
                                                    correctionFactorSun*legendre_sunTheta[n+1, m+1]*cos(m*sun_polar[1]))
                dSnm <- correctionFactorLove * (correctionFactorMoon*legendre_moonTheta[n+1, m+1]*sin(m*moon_polar[1]) + 
                                                    correctionFactorSun*legendre_sunTheta[n+1, m+1]*sin(m*sun_polar[1]))
                dCnmSET[n+1, m+1] <- dCnm
                dSnmSET[n+1, m+1] <- dSnm
            }
            # step 2 of corrections for solid Earth tides
            ## Calculation of theta_g
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
            #theta_g <- alternativeGmst(Mjd_UT1)
            theta_gPi <- theta_g + pi
            delaunayVars <- delaunayVariables(Time)
            theta_f0 <- vector(mode="numeric", length=nrow(asteRiskData::solidEarthTides_dC20))
            for(i in 1:nrow(asteRiskData::solidEarthTides_dC20)) {
                theta_f0[i] <- -drop(asteRiskData::solidEarthTides_dC20[i, 1:5] %*% delaunayVars)
                # m(thetag + pi) is 0 because m=0
            }
            dC20 <- 1e-12*sum(asteRiskData::solidEarthTides_dC20[, 6]*cos(theta_f0) - 
                                  asteRiskData::solidEarthTides_dC20[, 7]*sin(theta_f0))
            # 1e-12 multiplier because units of IERS table 6.5a and b are 1e-12
            # m = 1
            dCnmSET[3, 1] <- dCnmSET[3, 1] + dC20
            theta_f1 <- vector(mode="numeric", length=nrow(asteRiskData::solidEarthTides_dC21dS21))
            for(i in 1:nrow(asteRiskData::solidEarthTides_dC21dS21)) {
                theta_f1[i] <- theta_gPi - drop(asteRiskData::solidEarthTides_dC21dS21[i, 1:5] %*% delaunayVars)
            }
            dC21 <- 1e-12*sum(asteRiskData::solidEarthTides_dC21dS21[, 6]*sin(theta_f1) + 
                                  asteRiskData::solidEarthTides_dC21dS21[, 7]*cos(theta_f1))
            dS21 <- 1e-12*sum(asteRiskData::solidEarthTides_dC21dS21[, 6]*cos(theta_f1) - 
                                  asteRiskData::solidEarthTides_dC21dS21[, 7]*sin(theta_f1))
            dCnmSET[3,2] <- dCnmSET[3,2] + dC21
            dSnmSET[3,2] <- dSnmSET[3,2] + dS21
            # m = 2
            theta_f2 <- vector(mode="numeric", length=nrow(asteRiskData::solidEarthTides_dC22dS22))
            for(i in 1:nrow(asteRiskData::solidEarthTides_dC22dS22)) {
                theta_f2[i] <- 2*theta_gPi - drop(asteRiskData::solidEarthTides_dC21dS21[i, 1:5] %*% delaunayVars)
            }
            dC22 <- 1e-12*sum(asteRiskData::solidEarthTides_dC22dS22[, 6]*sin(theta_f2))
            dS22 <- 1e-12*sum(asteRiskData::solidEarthTides_dC22dS22[, 6]*cos(theta_f2))
            dCnmSET[3,3] <- dCnmSET[3,3] + dC22
            dSnmSET[3,3] <- dCnmSET[3,3] + dS22
            # step 3 of corrections for solid Earth tides: removal of duplicate zero-tide correction
            # necessary because the C20 of GGM0X is a zero-tide coefficient
            dC20perm <- 4.4228e-8*(-0.31460)*0.30190 # 0.30190 is the nominal value of k20
            dCnmSET[3,1] <- dCnmSET[3,1] - dC20perm
            # Solid Earth pole tides
            # Section 7.1.4 of IERS technical note 2018
            # secular x and y pole coordinates in milliarcseconds
            xs <- 55 + 1.677 * Time_years 
            ys <- 320.5 + 3.46 * Time_years
            # m1 and m2 parameters in arcseconds
            m1 <- x_pole * const_Arcs - xs/1000
            m2 <- -(y_pole * const_Arcs - ys/1000)
            dC21solidEarthPoleTide <- -1.333e-9*(m1 + 0.0115*m2)
            dS21solidEarthPoleTide <- -1.333e-9*(m2 - 0.0115*m2)
            dCnmSET[3,2] <- dCnmSET[3,2] + dC21solidEarthPoleTide
            dSnmSET[3,2] <- dSnmSET[3,2] + dS21solidEarthPoleTide
            # Add the corrections so far (i.e. all solid Earth tides)
            C[1:5, 1:5] <- C[1:5, 1:5] + dCnmSET
            S[1:5, 1:5] <- S[1:5, 1:5] + dSnmSET
        }
        # End corrections for solid Earth tides
        # Start corrections for ocean tides
        if(OTcorrections){
            # Calculate Doodson variables from Delaunay's
            if(SETcorrections) {
                doodsonVars <- delaunayToDoodsonVariables(delaunayVars, theta_g)
            } else {
                theta_g <- iauGmst06(iauCal2jd_results$DJMJD0, UT1, iauCal2jd_results$DJMJD0, TT)
                delaunayVars <- delaunayVariables(Time)
                doodsonVars <- delaunayToDoodsonVariables(delaunayVars, theta_g)
                
            }
            #doodsonVars <- doodsonVariables(TT)
            OTcorrected <- parallelOceanTidesCorrections(asteRiskData::oceanTides_fes2014_dCnmdSnm_tideNames,
                                                         asteRiskData::oceanTides_fes2014_dCnmdSnm_data,
                                                         doodsonVars, C, S, m1, m2)
            C <- OTcorrected[[1]]
            S <- OTcorrected[[2]]
        }
        # End of ocean tide correction
        # Ocean pole tide: currently only for C21 and S21, but should extend to degree 10
        # included in parallelOceanTidesCorrections in C++ code
        # End of all corrections
    }
    r_bf <- E %*% r
    d <- sqrt(sum(r_bf^2))
    latgc <- asin(r_bf[3]/d)
    lon <- atan2(r_bf[2], r_bf[1])
    # Define order of Legendre functions
    n <- m <- earthSPHDegree
    legendre_latgc <- legendre(n, m, latgc)
    legendre_latgc_Pnm <- legendre_latgc[[1]]
    legendre_latgc_dPnm <- legendre_latgc[[2]]
    gradient <- gravityGradientSphericalCoords(legendre_latgc_Pnm, legendre_latgc_dPnm,
                                              C, S, latgc, lon, d, earthRadius_EGM96, GM_Earth_TT, n, m)
    # for(n in 0:n) {
    #     b1 <- (-GM_Earth_TT/d^2)*(earthRadius_EGM96/d)^n*(n+1)
    #     b2 <- (GM_Earth_TT/d)*(earthRadius_EGM96/d)^n
    #     b3 <- (GM_Earth_TT/d)*(earthRadius_EGM96/d)^n
    #     q1 <- sum(legendre_latgc_Pnm[n+1,1:(m+1)]*(C[n+1,1:(m+1)]*cos((0:m)*lon)+S[n+1,1:(m+1)]*sin((0:m)*lon)))
    #     q2 <- sum(legendre_latgc_dPnm[n+1,1:(m+1)]*
    #                   (C[n+1,1:(m+1)]*cos((0:m)*lon)+S[n+1,1:(m+1)]*sin((0:m)*lon)))
    #     q3 <- sum((0:m) * legendre_latgc_Pnm[n+1,1:(m+1)]*(S[n+1,1:(m+1)]*cos((0:m)*lon)-C[n+1,1:(m+1)]*sin((0:m)*lon)))
    #     # for(m in 0:m) {
    #     #     q1 <- q1 + legendre_latgc$normLegendreValues[n+1,m+1]*(C[n+1,m+1]*cos(m*lon)+S[n+1,m+1]*sin(m*lon))
    #     #     q2 <- q2 + legendre_latgc$normLegendreDerivativeValues[n+1,m+1]*
    #     #         (C[n+1,m+1]*cos(m*lon)+S[n+1,m+1]*sin(m*lon))
    #     #     q3 <- q3 + m*legendre_latgc$normLegendreValues[n+1,m+1]*(S[n+1,m+1]*cos(m*lon)-C[n+1,m+1]*sin(m*lon))
    #     # }
    #     dUdr <- dUdr + q1*b1
    #     dUdlatgc <- dUdlatgc + q2*b2
    #     dUdlon <- dUdlon + q3*b3
    #     q3 <- 0
    #     q2 <- 0
    #     q1 <- 0
    # }
    dUdr <- gradient[1]
    dUdlatgc <- gradient[2]
    dUdlon <- gradient[3]
    r2xy <- r_bf[1]^2+r_bf[2]^2
    ax <- (1/d*dUdr-r_bf[3]/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf[1]-(1/r2xy*dUdlon)*r_bf[2]
    ay <- (1/d*dUdr-r_bf[3]/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf[2]+(1/r2xy*dUdlon)*r_bf[1]
    az <- 1/d*dUdr*r_bf[3]+sqrt(r2xy)/d^2*dUdlatgc
    a_bf <- c(ax, ay, az)
    a <- t(E) %*% a_bf
    return(as.vector(a))
}

elasticMoonAcceleration <- function(Mjd_UTC, r, r_sun, r_earth, moonLibrations, UT1_UTC,
                                    TT_UTC, moonSPHDegree, SMTcorrections) {
    ## GRAIL gravity fields use the Moon Principal Axes reference frame
    if((moonSPHDegree+1) > nrow(asteRiskData::GRGM1200B_Cnm)) {
        warning(strwrap(paste("Spherical harmonics coefficients are only available
                            up to degree ", nrow(asteRiskData::GRGM1200B_Cnm)-1, ", but ",
                              moonSPHDegree, " was input. Defaulting to the maximum
                            available degree and order of ", nrow(asteRiskData::GRGM1200B_Cnm)-1, ".")))
        moonSPHDegree <- nrow(asteRiskData::GRGM1200B_Cnm)-1
    }
    C <- asteRiskData::GRGM1200B_Cnm[1:(moonSPHDegree+1), 1:(moonSPHDegree+1)]
    S <- asteRiskData::GRGM1200B_Snm[1:(moonSPHDegree+1), 1:(moonSPHDegree+1)]
    ## Rotation from Moon Principal Axes frame to Lunar Celestial Reference Frame
    ## (Lunar cenetered ICRF). See Eq 8 at https://iopscience.iop.org/article/10.3847/1538-3881/abd414
    ## TODO : check if actually conversion from lunar frames ME to PA is required
    # Solid Moon Tides Corrections based on https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JE004755
    # These are anelastic corrections
    MPAtoLCRF <- iauRz(-moonLibrations[1], 
                       iauRx(-moonLibrations[2], 
                             iauRz(-moonLibrations[3], diag(3))
                       )
    )
    if(SMTcorrections) {
        r_earth <- t(MPAtoLCRF) %*% r_earth
        earth_polar <- CartesianToPolar(r_earth)
        r_sun <- t(MPAtoLCRF) %*% r_sun
        sun_polar <- CartesianToPolar(r_sun)
        # Mjd_TT <- Mjd_UTC + TT_UTC/86400
        # Mjd_UT1 <- Mjd_UTC + UT1_UTC/86400
        # Time <- (Mjd_TT-MJD_J2000)/36525
        MJD_TDB <- MJDUTCtoMJDTDB(Mjd_UTC)
        Centuries_TDB <- (MJD_TDB-MJD_J2000)/36525
        delaunayVars <- delaunayVariables(Centuries_TDB)
        delaunayVarsNoOmega <- delaunayVars[1:4]
        legendre_earthTheta <- legendre(2, 2, earth_polar[2])[[1]]
        legendre_sunTheta <- legendre(2, 2, sun_polar[2])[[1]]
        correctionFactorEarth <- GMratioEarthMoon*((moonRadius/earth_polar[3])^3)
        correctionFactorSun <- GMratioSunMoon*((moonRadius/sun_polar[3])^3)
        dC20 <- k20moon * (correctionFactorEarth * legendre_earthTheta[3, 1] +
                             correctionFactorSun * legendre_sunTheta[3, 1]) / 5
        dC21 <- k21moon * (correctionFactorEarth * legendre_earthTheta[3, 2] * cos(earth_polar[1]) +
                             correctionFactorSun * legendre_sunTheta[3, 2] * cos(sun_polar[1]))/5
        dC22 <- k22moon * (correctionFactorEarth * legendre_earthTheta[3, 3] * cos(2*earth_polar[1]) +
                             correctionFactorSun * legendre_sunTheta[3, 3] * cos(2*sun_polar[1]))/5
        dS21 <- k21moon * (correctionFactorEarth * legendre_earthTheta[3, 2] * sin(earth_polar[1]) +
                             correctionFactorSun * legendre_sunTheta[3, 2] * sin(sun_polar[1]))/5
        dS22 <- k22moon * (correctionFactorEarth * legendre_earthTheta[3, 3] * sin(2*earth_polar[1]) +
                             correctionFactorSun * legendre_sunTheta[3, 3] * sin(2*sun_polar[1]))/5
        # technically the following dissipation terms are anelastic 
        # dissTermsC20 <- dissTermsC21 <- dissTermsC22 <- dissTermsS21 <- dissTermsS22 <- numeric(nrow(solidMoonTidesSimple))
        # for(i in 1:length(dissTermsC20)) {
        #     solidMoonTidesRow <- solidMoonTidesSimple[i, ]
        #     argumentZeta <- drop(solidMoonTidesRow[1:4] %*% delaunayVarsNoOmega)
        #     dissTermsC20[i] <- solidMoonTidesRow[7] * sin(argumentZeta)
        #     dissTermsC21[i] <- solidMoonTidesRow[8] * cos(argumentZeta)
        #     dissTermsS21[i] <- solidMoonTidesRow[9] * sin(argumentZeta)
        #     dissTermsC22[i] <- solidMoonTidesRow[10] * sin(argumentZeta)
        #     dissTermsS22[i] <- solidMoonTidesRow[11] * cos(argumentZeta)
        # }
        # C[3,1] <- C[3,1] + dC20 + sum(dissTermsC20)
        # C[3,2] <- C[3,2] + dC21 + sum(dissTermsC21)
        # C[3,3] <- C[3,3] + dC22 + sum(dissTermsC22)
        # S[3,2] <- S[3,2] + dS21 + sum(dissTermsS21)
        # S[3,3] <- S[3,3] + dS22 + sum(dissTermsS22)
        C[3,1] <- C[3,1] + dC20 
        C[3,2] <- C[3,2] + dC21 
        C[3,3] <- C[3,3] + dC22 
        S[3,2] <- S[3,2] + dS21 
        S[3,3] <- S[3,3] + dS22 
    }
    pos_LCRF <- r
    pos_MPA <- t(MPAtoLCRF) %*% pos_LCRF
    d <- sqrt(sum(pos_MPA^2))
    latgc <- asin(pos_MPA[3]/d)
    lon <- atan2(pos_MPA[2], pos_MPA[1])
    # Define order of Legendre functions
    n <- m <- moonSPHDegree
    legendre_latgc <- legendre(n, m, latgc)
    legendre_latgc_Pnm <- legendre_latgc[[1]]
    legendre_latgc_dPnm <- legendre_latgc[[2]]
    gradient <- gravityGradientSphericalCoords(legendre_latgc_Pnm, legendre_latgc_dPnm,
                                               C, S, latgc, lon, d, moonRadius, GM_Moon_GRGM1200B, n, m)
    dUdr <- gradient[1]
    dUdlatgc <- gradient[2]
    dUdlon <- gradient[3]
    r2xy <- pos_MPA[1]^2+pos_MPA[2]^2
    ax <- (1/d*dUdr-pos_MPA[3]/(d^2*sqrt(r2xy))*dUdlatgc)*pos_MPA[1]-(1/r2xy*dUdlon)*pos_MPA[2]
    ay <- (1/d*dUdr-pos_MPA[3]/(d^2*sqrt(r2xy))*dUdlatgc)*pos_MPA[2]+(1/r2xy*dUdlon)*pos_MPA[1]
    az <- 1/d*dUdr*pos_MPA[3]+sqrt(r2xy)/d^2*dUdlatgc
    a_bf <- c(ax, ay, az)
    a <- MPAtoLCRF %*% a_bf
    return(as.vector(a))
}

anelasticMoonAcceleration <- function(Mjd_UTC, r, r_sun, r_earth, moonLibrations, UT1_UTC,
                                    TT_UTC, moonSPHDegree, SMTcorrections) {
    ## GRAIL gravity fields use the Moon Principal Axes reference frame
    if((moonSPHDegree+1) > nrow(asteRiskData::GRGM1200B_Cnm)) {
        warning(strwrap(paste("Spherical harmonics coefficients are only available
                            up to degree ", nrow(asteRiskData::GRGM1200B_Cnm)-1, ", but ",
                              moonSPHDegree, " was input. Defaulting to the maximum
                            available degree and order of ", nrow(asteRiskData::GRGM1200B_Cnm)-1, ".")))
        moonSPHDegree <- nrow(asteRiskData::GRGM1200B_Cnm)-1
    }
    C <- asteRiskData::GRGM1200B_Cnm[1:(moonSPHDegree+1), 1:(moonSPHDegree+1)]
    S <- asteRiskData::GRGM1200B_Snm[1:(moonSPHDegree+1), 1:(moonSPHDegree+1)]
    ## Rotation from Moon Principal Axes frame to Lunar Celestial Reference Frame
    ## (Lunar cenetered ICRF). See Eq 8 at https://iopscience.iop.org/article/10.3847/1538-3881/abd414
    ## TODO : check if actually conversion from lunar frames ME to PA is required
    # Solid Moon Tides Corrections based on https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JE004755
    # These are anelastic corrections
    MPAtoLCRF <- iauRz(-moonLibrations[1], 
                       iauRx(-moonLibrations[2], 
                             iauRz(-moonLibrations[3], diag(3))
                       )
    )
    if(SMTcorrections) {
        dCnmSMT <- matrix(0, nrow=3, ncol=3)
        dSnmSMT <- matrix(0, nrow=3, ncol=3)
        r_earth <- t(MPAtoLCRF) %*% r_earth
        earth_polar <- CartesianToPolar(r_earth)
        r_sun <- t(MPAtoLCRF) %*% r_sun
        sun_polar <- CartesianToPolar(r_sun)
        # Mjd_TT <- Mjd_UTC + TT_UTC/86400
        # Mjd_UT1 <- Mjd_UTC + UT1_UTC/86400
        # Time <- (Mjd_TT-MJD_J2000)/36525
        MJD_TDB <- MJDUTCtoMJDTDB(Mjd_UTC)
        Centuries_TDB <- (MJD_TDB-MJD_J2000)/36525
        delaunayVars <- delaunayVariables(Centuries_TDB)
        delaunayVarsNoOmega <- delaunayVars[1:4]
        legendre_earthTheta <- legendre(2, 2, earth_polar[2])[[1]]
        legendre_sunTheta <- legendre(2, 2, sun_polar[2])[[1]]
        correctionFactorEarth <- GMratioEarthMoon*((moonRadius/earth_polar[3])^3)
        correctionFactorSun <- GMratioSunMoon*((moonRadius/sun_polar[3])^3)
        # dC20 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 1] +
        #                      correctionFactorSun * legendre_sunTheta[3, 1])
        # dC21 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 2] * cos(earth_polar[1]) +
        #                      correctionFactorSun * legendre_sunTheta[3, 2] * cos(sun_polar[1]))/3
        # dC22 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 3] * cos(2*earth_polar[1]) +
        #                      correctionFactorSun * legendre_sunTheta[3, 3] * cos(2*sun_polar[1]))/12
        # dS21 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 2] * sin(earth_polar[1]) +
        #                      correctionFactorSun * legendre_sunTheta[3, 2] * sin(sun_polar[1]))/3
        # dS22 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 3] * sin(2*earth_polar[1]) +
        #                      correctionFactorSun * legendre_sunTheta[3, 3] * sin(2*sun_polar[1]))/12
        dC20 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 1] +
                             correctionFactorSun * legendre_sunTheta[3, 1])/(legendreNormFactor(2,0)^2)
        dC21 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 2] * cos(earth_polar[1]) +
                             correctionFactorSun * legendre_sunTheta[3, 2] * cos(sun_polar[1]))/(3*legendreNormFactor(2,1)^2) 
        dC22 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 3] * cos(2*earth_polar[1]) +
                             correctionFactorSun * legendre_sunTheta[3, 3] * cos(2*sun_polar[1]))/(12*legendreNormFactor(2,2)^2)
        dS21 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 2] * sin(earth_polar[1]) +
                             correctionFactorSun * legendre_sunTheta[3, 2] * sin(sun_polar[1]))/(3*legendreNormFactor(2,1)^2)
        dS22 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 3] * sin(2*earth_polar[1]) +
                             correctionFactorSun * legendre_sunTheta[3, 3] * sin(2*sun_polar[1]))/(12*legendreNormFactor(2,2)^2)
        debugTermsC20 <- debugTermsC21 <- debugTermsC22 <- debugTermsS21 <- debugTermsS22 <- numeric(nrow(asteRiskData::solidMoonTides))
        dissTermsC20 <- dissTermsC21 <- dissTermsC22 <- dissTermsS21 <- dissTermsS22 <- numeric(nrow(asteRiskData::solidMoonTides))
        for(i in 1:length(dissTermsC20)) {
            solidMoonTidesRow <- asteRiskData::solidMoonTides[i, ]
            v <- (2*pi)/solidMoonTidesRow[6]
            gFunctionDebye <- ((vref^2 -v^2) * tau2)/(vref * (1 + v^2 * tau2^2))
            fFunctionDebye <- (v*(1+vref^2*tau2^2))/(vref * (1 + v^2 * tau2^2))
            argumentZeta <- drop(solidMoonTidesRow[1:4] %*% delaunayVarsNoOmega)
            gfFuncsTerm1 <- (gFunctionDebye * cos(argumentZeta) + fFunctionDebye * sin(argumentZeta))
            gfFuncsTerm2 <- (gFunctionDebye * sin(argumentZeta) - fFunctionDebye * cos(argumentZeta))
            dissTermsC20[i] <- solidMoonTidesRow[7] * gfFuncsTerm1
            dissTermsC21[i] <- solidMoonTidesRow[8] * gfFuncsTerm2
            dissTermsS21[i] <- solidMoonTidesRow[9] * gfFuncsTerm1
            dissTermsC22[i] <- solidMoonTidesRow[10] * gfFuncsTerm1
            dissTermsS22[i] <- solidMoonTidesRow[11] * gfFuncsTerm2
            debugTermsC20[i] <- solidMoonTidesRow[7] * cos(argumentZeta)
            debugTermsC21[i] <- solidMoonTidesRow[8] * sin(argumentZeta)
            debugTermsS21[i] <- solidMoonTidesRow[9] * cos(argumentZeta)
            debugTermsC22[i] <- solidMoonTidesRow[10] * cos(argumentZeta)
            debugTermsS22[i] <- solidMoonTidesRow[11] * sin(argumentZeta)
        }
        dC20 <- (dC20 + k2ByQref * sum(dissTermsC20) * 10e-9)
        dC21 <- (dC21 + k2ByQref * sum(dissTermsC21) * 10e-9)
        dS21 <- (dS21 + k2ByQref * sum(dissTermsS21) * 10e-9)
        dC22 <- (dC22 + k2ByQref * sum(dissTermsC22) * 10e-9)
        dS22 <- (dS22 + k2ByQref * sum(dissTermsS22) * 10e-9)
        # dC20 <- (sum(debugTermsC20) * k2ref + k2ByQref * sum(dissTermsC20)) * 10e-9
        # dC21 <- (sum(debugTermsC21) * k2ref + k2ByQref * sum(dissTermsC21)) * 10e-9
        # dS21 <- (sum(debugTermsS21) * k2ref + k2ByQref * sum(dissTermsS21)) * 10e-9
        # dC22 <- (sum(debugTermsC22) * k2ref + k2ByQref * sum(dissTermsC22)) * 10e-9
        # dS22 <- (sum(debugTermsS22) * k2ref + k2ByQref * sum(dissTermsS22)) * 10e-9
        C[3,1] <- C[3,1] + dC20
        C[3,2] <- C[3,2] + dC21
        C[3,3] <- C[3,3] + dC22
        S[3,2] <- S[3,2] + dS21
        S[3,3] <- S[3,3] + dS22
        # dC20 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 1] +
        #                      correctionFactorSun * legendre_sunTheta[3, 1]) / 5
        # dC21 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 2] * cos(earth_polar[1]) +
        #                      correctionFactorSun * legendre_sunTheta[3, 2] * cos(sun_polar[1]))/5
        # dC22 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 3] * cos(2*earth_polar[1]) +
        #                      correctionFactorSun * legendre_sunTheta[3, 3] * cos(2*sun_polar[1]))/5
        # dS21 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 2] * sin(earth_polar[1]) +
        #                      correctionFactorSun * legendre_sunTheta[3, 2] * sin(sun_polar[1]))/5
        # dS22 <- k2ref * (correctionFactorEarth * legendre_earthTheta[3, 3] * sin(2*earth_polar[1]) +
        #                      correctionFactorSun * legendre_sunTheta[3, 3] * sin(2*sun_polar[1]))/5
        # C[3,1] <- C[3,1] + dC20
        # C[3,2] <- C[3,2] + dC21
        # C[3,3] <- C[3,3] + dC22
        # S[3,2] <- S[3,2] + dS21
        # S[3,3] <- S[3,3] + dS22
    }
    pos_LCRF <- r
    pos_MPA <- t(MPAtoLCRF) %*% pos_LCRF
    d <- sqrt(sum(pos_MPA^2))
    latgc <- asin(pos_MPA[3]/d)
    lon <- atan2(pos_MPA[2], pos_MPA[1])
    # Define order of Legendre functions
    n <- m <- moonSPHDegree
    legendre_latgc <- legendre(n, m, latgc)
    legendre_latgc_Pnm <- legendre_latgc[[1]]
    legendre_latgc_dPnm <- legendre_latgc[[2]]
    gradient <- gravityGradientSphericalCoords(legendre_latgc_Pnm, legendre_latgc_dPnm,
                                               C, S, latgc, lon, d, moonRadius, GM_Moon_GRGM1200B, n, m)
    dUdr <- gradient[1]
    dUdlatgc <- gradient[2]
    dUdlon <- gradient[3]
    r2xy <- pos_MPA[1]^2+pos_MPA[2]^2
    ax <- (1/d*dUdr-pos_MPA[3]/(d^2*sqrt(r2xy))*dUdlatgc)*pos_MPA[1]-(1/r2xy*dUdlon)*pos_MPA[2]
    ay <- (1/d*dUdr-pos_MPA[3]/(d^2*sqrt(r2xy))*dUdlatgc)*pos_MPA[2]+(1/r2xy*dUdlon)*pos_MPA[1]
    az <- 1/d*dUdr*pos_MPA[3]+sqrt(r2xy)/d^2*dUdlatgc
    a_bf <- c(ax, ay, az)
    a <- MPAtoLCRF %*% a_bf
    return(as.vector(a))
}

pointMassAcceleration <- function(r, s, GM) {
    # difference position vector
    d <- r - s
    a <- -GM * (d/(sqrt(sum(d^2)))^3 + s/(sqrt(sum(s^2)))^3)
    return(a)
}

solarRadiationAcceleration_deprecated <- function(r, r_earth, r_moon, r_sun, r_sunSSB, 
                                       area, mass, Cr, P0, AU, shm) {
    pccor <- r_moon # position of the moon relative to Earth
    ccor <- r_earth - r_sunSSB # position of the Solar System Barycenter relative to Sun
    pscor <- r_moon - r_sun # Position of Moon relative to Sun
    sbcor <- r # position of Satellite relative to Earth
    bcor <- r - r_sun # Position of Satellite relative to Sun
    sbpcor <- r - r_moon # Position of Satellite relative to Moon
    if(shm == "cylindrical") {
        warning("Using cylindrical shadow model for solar radiation pressure")
        nu <- cylindricalShadow(r, r_sun)
    } else {
        nu <- geometricShadow(pccor,ccor,pscor,sbcor,bcor,sbpcor)$lambda
    }
    a <- nu*Cr*(area/mass)*P0*(AU*AU)*bcor/(norm(bcor, type="2")^3)
    return(a)
}

solarRadiationAcceleration <- function(rSatellite, JPLEphemerides, centralBody, 
                                        area, mass, Cr, solarPressureConstant, AU) {
    # rCentralBody <- JPLEphemerides[[paste("position", centralBody, sep="")]]
    # in principle should be 0 but doing it this way allows using ephemerides with different central body
    rSun <- JPLEphemerides[[paste("position", "Sun", sep="")]]
    # Variables starting with r are position vectors
    # Variables starting with d are distances (norm of position vectors)
    # variables starting with alpha are apparent sizes
    # variables starting with beta are angular separations
    if(centralBody %in% c("Earth", "Moon")) {
        rEarth <- JPLEphemerides[[paste("position", "Earth", sep="")]]
        rEarth_Sun <- rEarth - rSun 
        rMoon <- JPLEphemerides[[paste("position", "Moon", sep="")]]
        rMoon_Sun <- rMoon - rSun
        rSatellite_Sun <- rSatellite - rSun
        rSatellite_Moon <- rSatellite - rMoon
        rSatellite_Earth <- rSatellite - rEarth
        dSatellite_Sun <- sqrt(sum(rSatellite_Sun^2))
        dEarth_Sun <- sqrt(sum(rEarth_Sun^2))
        dMoon_Sun <- sqrt(sum(rMoon_Sun^2))
        dSatellite_Earth <- sqrt(sum(rSatellite_Earth^2))
        alphaSun <- asin(sunRadius/dSatellite_Sun) # apparent radius of Sun from satellite
        if(dSatellite_Sun > dEarth_Sun) {
            # only if distance from Sun to satellite is larger than from Sun
            # to Earth there is a chance that there is occultation of Sun by Earth
            alphaEarth <- asin(earthRadius_EGM96/dSatellite_Earth) # apparent radius of Earth from satellite
            betaEarthSun <- acos((rSatellite_Sun %*% rSatellite_Earth)/(dSatellite_Sun * dSatellite_Earth)) 
            # angular separation Earth-Sun from satellite
            if(betaEarthSun < (alphaEarth - alphaSun)) {
                # Earth completely eclipses Sun, no need to check for Moon eclipse
                eclipseFactor <- 0
                eclipseType <- "EarthTotal"
            } else {
                dSatellite_Moon <- sqrt(sum(rSatellite_Moon^2))
                alphaMoon <- asin(moonRadius/dSatellite_Moon)
                betaMoonSun <- acos((rSatellite_Sun %*% rSatellite_Moon)/(dSatellite_Sun * dSatellite_Moon)) 
                # angular separation Moon-Sun from satellite
                if (betaEarthSun > (alphaEarth + alphaSun)) {
                    # Earth does not eclipse Sun at all, but we need to check for Moon eclipse
                    if(betaMoonSun < (alphaMoon - alphaSun)) {
                        # Moon completely eclipses Sun
                        eclipseFactor <- 0
                        eclipseType <- "MoonTotal"
                    } else {
                        if(betaMoonSun > (alphaMoon + alphaSun)) {
                            # Moon does not eclipse Sun at all either
                            eclipseFactor <- 1
                            eclipseType <- "None"
                        } else {
                            # We are therefore in a situation that:
                            # (alphaMoon - alphaSun) < betaMoonSun < (alphaMoon + alphaSun)
                            # which means there is partial occultation (penumbra) of Sun by Moon
                            a <- min(alphaMoon, alphaSun)
                            b <- max(alphaMoon, alphaSun) # in principle b will be Moon apparent size
                            c <- betaMoonSun
                            luneArea <- calculateLuneArea(a, b, c)
                            if(alphaMoon >= alphaSun) {
                                # Angular diameter of Moon > of Sun, so the fraction of
                                # visible Sun is the lune area
                                eclipseFactor <- luneArea/(pi*alphaSun^2)
                                eclipseType <- "MoonPartial-1"
                            } else {
                                # Angular diameter of Moon < of Sun, so we need
                                # to subtract lune area from Moon apparent area
                                # and then subtract that from Sun apparent area
                                moonOccultingArea <- pi*alphaMoon^2 - luneArea
                                eclipseFactor <- (pi*alphaSun^2 - moonOccultingArea)/(pi*alphaSun^2)
                                eclipseType <- "MoonPartial-2"
                            }
                        }
                    }
                } else {
                    # This is a situation where:
                    # (alphaEarth - alphaSun) < betaEarthSun < (alphaEarth + alphaSun)
                    # Which means there is partial eclipse caused by Earth.
                    # We now need to check if there is also occultation by Moon
                    if(betaMoonSun < (alphaMoon - alphaSun)) {
                        # Moon completely eclipses Sun, therefore partial occultation
                        # by Earth irrelevant. 
                        # TODO: refactor control flow for situations
                        eclipseFactor <- 0
                        eclipseType <- "MoonTotal"
                    } else {
                        if(betaMoonSun > (alphaMoon + alphaSun)) {
                            # Moon does not eclipse Sun at all, so we just calculate
                            # partial eclipse by Earth
                            a <- min(alphaEarth, alphaSun)
                            b <- max(alphaEarth, alphaSun) # in principle b will be Moon apparent size
                            c <- betaEarthSun
                            luneArea <- calculateLuneArea(a, b, c)
                            if(alphaEarth >= alphaSun) {
                                # Angular diameter of Earth > of Sun, so the fraction of
                                # visible Sun is the lune area
                                eclipseFactor <- luneArea/(pi*alphaSun^2)
                                eclipseType <- "EarthPartial-1"
                            } else {
                                # Angular diameter of Earth < of Sun, so we need
                                # to subtract lune area from Earth apparent area
                                # and then subtract that from Sun apparent area
                                earthOccultingArea <- pi*alphaEarth^2 - luneArea
                                eclipseFactor <- (pi*alphaSun^2 - earthOccultingArea)/(pi*alphaSun^2)
                                eclipseType <- "EarthPartial-2"
                            }
                        } else {
                            # This is the most complex situation. There is partial
                            # occultation caused by both the Earth and Moon
                            # Simpler case: Earth and Moon disks do not overlap,
                            # their occultations of Sun are independent. We need to check
                            # eclipse between Earth and Moon
                            betaMoonEarth <- acos((rSatellite_Earth %*% rSatellite_Moon)/(dSatellite_Earth * dSatellite_Moon))
                            # This is the angular separation between Moon and Earth from satellite
                            if(betaMoonEarth > (alphaMoon + alphaEarth)) {
                                # Earth and Moon do not overlap at all. So we calculate
                                # their contributions to Sun occultation independently
                                a1 <- min(alphaEarth, alphaSun)
                                b1 <- max(alphaEarth, alphaSun) 
                                c1 <- betaEarthSun
                                luneArea_Earth <- calculateLuneArea(a1, b1, c1)
                                sunDiskArea <- pi*alphaSun^2
                                if(alphaEarth >= alphaSun) {
                                    # Angular diameter of Earth > of Sun, so the area of
                                    # visible Sun is the lune area (considering only Earth eclipse)
                                    sunVisibleArea <- luneArea_Earth
                                } else {
                                    # Angular diameter of Earth < of Sun, so we need
                                    # to subtract lune area from Earth apparent area
                                    # and then subtract that from Sun apparent area
                                    # to get the area of visible Sun considering
                                    # only Earth eclipse
                                    earthOccultingArea <- pi*alphaEarth^2 - luneArea_Earth
                                    sunVisibleArea <- sunDiskArea - earthOccultingArea
                                }
                                # And now apply Moon contribution
                                a2 <- min(alphaMoon, alphaSun)
                                b2 <- max(alphaMoon, alphaSun) 
                                c2 <- betaMoonSun
                                luneArea_Moon <- calculateLuneArea(a2, b2, c2)
                                if(alphaMoon >= alphaSun) {
                                    # Angular diameter of Moon > of Sun, so the area of
                                    # visible Sun is the lune area (considering only Moon eclipse)
                                    # What we will do is, subtract such lune area from 
                                    # area of total Sun disk, to obtain the area eclipsed by
                                    # Moon, and then subtract this from the previously calculated
                                    # Sun visible area considering Earth eclipse
                                    sunVisibleArea <- sunVisibleArea - (sunDiskArea - luneArea_Moon)
                                } else {
                                    # Angular diameter of Moon < of Sun, so we need
                                    # to subtract lune area from Moon apparent area
                                    # and then subtract that from Sun cisible area
                                    # to get the area of visible Sun considering
                                    # both Earth and Moon eclipses
                                    moonOccultingArea <- pi*alphaMoon^2 - luneArea_Earth
                                    sunVisibleArea <- sunVisibleArea - moonOccultingArea
                                }
                                eclipseFactor <- sunVisibleArea/sunDiskArea
                                eclipseType <- "EarthMoonMixed-1"
                            } else {
                                if(betaMoonEarth < (alphaMoon - alphaEarth)) {
                                    # Moon completely eclipses Earth, so we only 
                                    # need to calculate occultation by Moon
                                    a <- min(alphaMoon, alphaSun)
                                    b <- max(alphaMoon, alphaSun) # in principle b will be Moon apparent size
                                    c <- betaMoonSun
                                    luneArea <- calculateLuneArea(a, b, c)
                                    if(alphaMoon >= alphaSun) {
                                        # Angular diameter of Moon > of Sun, so the fraction of
                                        # visible Sun is the lune area
                                        eclipseFactor <- luneArea/(pi*alphaSun^2)
                                        eclipseType <- "MoonPartial-1"
                                    } else {
                                        # Angular diameter of Moon < of Sun, so we need
                                        # to subtract lune area from Moon apparent area
                                        # and then subtract that from Sun apparent area
                                        moonOccultingArea <- pi*alphaMoon^2 - luneArea
                                        eclipseFactor <- (pi*alphaSun^2 - moonOccultingArea)/(pi*alphaSun^2)
                                        eclipseType <- "MoonPartial-2"
                                    }
                                } else if(betaMoonEarth < (alphaEarth - alphaMoon)) {
                                    # Earth completely eclipses Moon, so we only need
                                    # to calculate occultation by Earth
                                    a <- min(alphaEarth, alphaSun)
                                    b <- max(alphaEarth, alphaSun) # in principle b will be Moon apparent size
                                    c <- betaEarthSun
                                    luneArea <- calculateLuneArea(a, b, c)
                                    if(alphaEarth >= alphaSun) {
                                        # Angular diameter of Earth > of Sun, so the fraction of
                                        # visible Sun is the lune area
                                        eclipseFactor <- luneArea/(pi*alphaSun^2)
                                        eclipseType <- "EarthPartial-1"
                                    } else {
                                        # Angular diameter of Earth < of Sun, so we need
                                        # to subtract lune area from Earth apparent area
                                        # and then subtract that from Sun apparent area
                                        earthOccultingArea <- pi*alphaEarth^2 - luneArea
                                        eclipseFactor <- (pi*alphaSun^2 - earthOccultingArea)/(pi*alphaSun^2)
                                        eclipseType <- "EarthPartial-2"
                                    }
                                } else {
                                    # This is the worst-case scenario. Both Earth and Moon
                                    # are partially eclipsing Sun, and they are partially eclipsing
                                    # each other as well. If we just subtracted their occultation
                                    # of the Sun independently, we would be double-subtracting
                                    # some occulted area. Need to work on getting this correctly
                                    # but temporarily, treat them as independently
                                    # TODO: fix this
                                    a1 <- min(alphaEarth, alphaSun)
                                    b1 <- max(alphaEarth, alphaSun) 
                                    c1 <- betaEarthSun
                                    luneArea_Earth <- calculateLuneArea(a1, b1, c1)
                                    sunDiskArea <- pi*alphaSun^2
                                    if(alphaEarth >= alphaSun) {
                                        # Angular diameter of Earth > of Sun, so the area of
                                        # visible Sun is the lune area (considering only Earth eclipse)
                                        sunVisibleArea <- luneArea_Earth
                                    } else {
                                        # Angular diameter of Earth < of Sun, so we need
                                        # to subtract lune area from Earth apparent area
                                        # and then subtract that from Sun apparent area
                                        # to get the area of visible Sun considering
                                        # only Earth eclipse
                                        earthOccultingArea <- pi*alphaEarth^2 - luneArea_Earth
                                        sunVisibleArea <- sunDiskArea - earthOccultingArea
                                    }
                                    # And now apply Moon contribution
                                    a2 <- min(alphaMoon, alphaSun)
                                    b2 <- max(alphaMoon, alphaSun) 
                                    c2 <- betaMoonSun
                                    luneArea_Moon <- calculateLuneArea(a2, b2, c2)
                                    if(alphaMoon >= alphaSun) {
                                        # Angular diameter of Moon > of Sun, so the area of
                                        # visible Sun is the lune area (considering only Moon eclipse)
                                        # What we will do is, subtract such lune area from 
                                        # area of total Sun disk, to obtain the area eclipsed by
                                        # Moon, and then subtract this from the previously calculated
                                        # Sun visible area considering Earth eclipse
                                        sunVisibleArea <- sunVisibleArea - (sunDiskArea - luneArea_Moon)
                                    } else {
                                        # Angular diameter of Moon < of Sun, so we need
                                        # to subtract lune area from Moon apparent area
                                        # and then subtract that from Sun cisible area
                                        # to get the area of visible Sun considering
                                        # both Earth and Moon eclipses
                                        moonOccultingArea <- pi*alphaMoon^2 - luneArea_Earth
                                        sunVisibleArea <- sunVisibleArea - moonOccultingArea
                                    }
                                    eclipseFactor <- sunVisibleArea/sunDiskArea
                                    eclipseType <- "EarthMoonMixed-2"
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        eclipseFactor <- 1
        eclipseType <- "None"
    }
    eclipseFactor <- as.numeric(eclipseFactor)
    a <- eclipseFactor*Cr*(area/mass)*solarPressureConstant*(AU*AU)*rSatellite_Sun/(sqrt(sum(rSatellite_Sun^2))^3)
    return(a)
}


dragAcceleration <- function(dens, r, v, T, area, mass, CD, Omega) {
    omega = c(0, 0, Omega)
    r_bf <- T %*% r
    v_bf <- T %*% v
    v_rel <- v_bf - c(omega[2]*r_bf[3] - omega[3]*r_bf[2],
                       omega[3]*r_bf[1] - omega[1]*r_bf[3],
                       omega[1]*r_bf[2] - omega[2]*r_bf[1])
    modulo_v_rel <- sqrt(sum(v_rel^2))
    a_bf <- -0.5*CD*(area/mass)*dens*modulo_v_rel*v_rel
    a <- t(T) %*% a_bf
    return(a)
}

relativity <- function(r, v, GM_centralBody) {
    r_Sat <- sqrt(sum(r^2))
    v_Sat <- sqrt(sum(v^2))
    a <- GM_centralBody/(c_light^2*r_Sat^3)*((4*GM_centralBody/r_Sat-v_Sat^2)*r+4*as.vector((r %*% v))*v)
    return(a)
}

accel <- function(t, Y, MJD_UTC, solarArea, satelliteMass, satelliteArea, Cr, Cd, 
                  earthSPHDegree, SETcorrections, OTcorrections, moonSPHDegree, 
                  SMTcorrections, centralBody, autoCentralBodyChange) {
    MJD_UTC <- MJD_UTC + t/86400
    MJD_TT <- MJDUTCtoMJDTT(MJD_UTC)
    MJD_TDB <- MJDUTCtoMJDTDB(MJD_UTC)
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, "l")
    x_pole <- IERS_results$x_pole[[1]]
    y_pole <- IERS_results$y_pole[[1]]
    UT1_UTC <- IERS_results$UT1_UTC[[1]]
    MJD_UT1 <- MJD_UTC + UT1_UTC/86400
    LOD <- IERS_results$LOD[[1]]
    dpsi <- IERS_results$dpsi[[1]]
    deps <- IERS_results$deps[[1]]
    dx_pole <- IERS_results$dx_pole[[1]]
    dy_pole <- IERS_results$dy_pole[[1]]
    TAI_UTC <- IERS_results$TAI_UTC[[1]]
    timeDiffs_results <- timeDiffs(UT1_UTC, TAI_UTC)
    UT1_TAI <- timeDiffs_results$UT1_TAI[[1]]
    UTC_GPS <- timeDiffs_results$UTC_GPS[[1]]
    UT1_GPS <- timeDiffs_results$UT1_GPS[[1]]
    TT_UTC <- timeDiffs_results$TT_UTC[[1]]
    GPS_UTC <- timeDiffs_results$GPS_UTC[[1]]
    # Add 2400000.5 to convert modified Julian date to Julian date
    # invjday_results <- invjday(MJD_UTC+2400000.5)
    # year <- invjday_results$year
    # month <- invjday_results$month
    # day <- invjday_results$day
    # hour <- invjday_results$hour
    # minute <- invjday_results$min
    # sec <- invjday_results$sec
    # iauCal2jd_results <- iauCal2jd(year, month, day)
    # DJMJD0 <-iauCal2jd_results$DJMJD0
    # DATE <- iauCal2jd_results$DATE
    # TIME <- (60*(60*hour + minute) + sec)/86400
    # UTC <- DATE + TIME
    # TT <- UTC + TT_UTC/86400
    # TUT <- TIME + UT1_UTC/86400
    # UT1 <- DATE + TUT
    #PMM <- iauPom00(x_pole, y_pole, iauSp00(DJMJD0, TT))
    PMM <- iauPom00(x_pole, y_pole, iauSp00(2400000.5, MJD_TT))
    # NPB <- iauPnm06a(DJMJD0, TT)
    NPB <- iauPnm06a(2400000.5, MJD_TT)
    # theta <- iauRz(iauGst06(DJMJD0, UT1, DJMJD0, TT, NPB), diag(3))
    theta <- iauRz(iauGst06(2400000.5, MJD_UT1, 2400000.5, MJD_TT, NPB), diag(3))
    E <- PMM %*% theta %*% NPB
    # MJD_TDB <- Mjday_TDB(TT)
    JPL_ephemerides <- JPLephemeridesDE440(MJD_TDB, centralBody, derivativesOrder=2)
    if(autoCentralBodyChange) {
        # realCentralBody <- determineCentralBody(Y[1:3], JPL_ephemerides[-c(1, 2, 12, 13)], JPL_ephemerides[[12]])
        realCentralBody <- determineCentralBody(Y[1:3], JPL_ephemerides[c("positionMercury",
                                                                          "positionVenus",
                                                                          "positionEarth",
                                                                          "positionMars",
                                                                          "positionJupiter",
                                                                          "positionSaturn",
                                                                          "positionUranus",
                                                                          "positionNeptune",
                                                                          "positionPluto")], 
                                                JPL_ephemerides[["positionMoon"]])
    }
    else {
        realCentralBody <- centralBody
    }
    if(centralBody == "Earth") {
        # Acceleration due to Earth, with zonal harmonics
        a <- anelasticEarthAcceleration(MJD_UTC, JPL_ephemerides$positionSun,
                                        JPL_ephemerides$positionMoon, Y[1:3], E,
                                        UT1_UTC, TT_UTC,
                                        x_pole, y_pole, earthSPHDegree, SETcorrections,
                                        OTcorrections)
        # Acceleration due to Moon as point mass only if Moon spherical degree is 1
        if(moonSPHDegree == 1) {
            a <- a + pointMassAcceleration(Y[1:3],JPL_ephemerides$positionMoon,GM_Moon_DE440)
        } else {
            a <- a + elasticMoonAcceleration(MJD_UTC, Y[1:3] - JPL_ephemerides$positionMoon, 
                                               JPL_ephemerides$positionSun - JPL_ephemerides$positionMoon, 
                                               -JPL_ephemerides$positionMoon,
                                               JPL_ephemerides$lunarLibrationAngles, 
                                               UT1_UTC, TT_UTC, moonSPHDegree, 
                                               SMTcorrections) - 
                GM_Moon_GRGM1200B * JPL_ephemerides$positionMoon/(sqrt(sum(JPL_ephemerides$positionMoon^2)))^3
        }
        # Acceleration due to solar radiation pressure
        a <- a + solarRadiationAcceleration(Y[1:3], JPL_ephemerides$positionEarth, 
                                            JPL_ephemerides$positionMoon, 
                                            JPL_ephemerides$positionSun,
                                            JPL_ephemerides$positionSunSSBarycentric, 
                                            solarArea,
                                            satelliteMass,
                                            Cr, solarPressureConst, AU, "geometrical")
        # Acceleration due to atmospheric drag
        # Omega according to section III from https://hpiers.obspm.fr/iers/bul/bulb/explanatory.html
        Omega <- omegaEarth - 8.43994809e-10*LOD
        dens <- NRLMSISE00(MJD_UTC, E%*%Y[1:3], UT1_UTC, TT_UTC)["Total"]
        # a <- a + dragAcceleration(dens, Y[1:3], Y[4:6], NPB, satelliteArea, 
        #                           satelliteMass, Cd, Omega)
        a <- a + dragAcceleration(dens, Y[1:3], Y[4:6], E, satelliteArea, 
                                  satelliteMass, Cd, Omega)
        # Relativistic effects
        a <- a + relativity(Y[1:3], Y[4:6], GM_Earth_TDB)
    } else if(centralBody == "Moon") {
        # Acceleration due to Moon with spherical harmonics
        a <- elasticMoonAcceleration(MJD_UTC, Y[1:3], JPL_ephemerides$positionSun, 
                                        JPL_ephemerides$positionEarth,
                                        JPL_ephemerides$lunarLibrationAngles, 
                                        UT1_UTC, TT_UTC, moonSPHDegree, 
                                        SMTcorrections)
        # Acceleration due to Earth as point mass if spherical harmonics degree is set to 1
        if(earthSPHDegree == 1) {
            a <- a + pointMassAcceleration(Y[1:3],JPL_ephemerides$positionEarth,GM_Earth_DE440)
        }
        else {
            a <- a + anelasticEarthAcceleration(MJD_UTC, JPL_ephemerides$positionSun - JPL_ephemerides$positionEarth,
                                                JPL_ephemerides$positionMoon - JPL_ephemerides$positionEarth,
                                                Y[1:3] - JPL_ephemerides$positionEarth, E,
                                                UT1_UTC, TT_UTC, x_pole, y_pole, earthSPHDegree,
                                                #SETcorrections, OTcorrections) - GM_Earth_TT * JPL_ephemerides$positionEarth/(sqrt(sum(JPL_ephemerides$positionEarth^2)))^3
                                                SETcorrections, OTcorrections) - #JPL_ephemerides$accelerationEarth
                anelasticEarthAcceleration(MJD_UTC, JPL_ephemerides$positionSun - JPL_ephemerides$positionEarth,
                                           JPL_ephemerides$positionMoon - JPL_ephemerides$positionEarth,
                                           0 - JPL_ephemerides$positionEarth, E,
                                           UT1_UTC, TT_UTC, x_pole, y_pole, earthSPHDegree,
                                           #SETcorrections, OTcorrections) - GM_Earth_TT * JPL_ephemerides$positionEarth/(sqrt(sum(JPL_ephemerides$positionEarth^2)))^3
                                           SETcorrections, OTcorrections)
            
        }
        a <- a + solarRadiationAcceleration(Y[1:3], JPL_ephemerides, "Moon", solarArea,
                                             satelliteMass, Cr, solarPressureConst, AU)
        a <- a + relativity(Y[1:3], Y[4:6], GM_Moon_GRGM1200B)
    } else {
        # Acceleration due to Earth and Moon as point masses
        a <- pointMassAcceleration(Y[1:3], JPL_ephemerides$positionEarth, GM_Earth_DE440)
        a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionMoon, GM_Moon_DE440)
    }
    # Acceleration due to Sun
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionSun,GM_Sun_DE440)
    # Accelerations due to planets
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionMercury, GM_Mercury_DE440)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionVenus, GM_Venus_DE440)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionMars, GM_Mars_DE440)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionJupiter, GM_Jupiter_DE440)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionSaturn, GM_Saturn_DE440)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionUranus, GM_Uranus_DE440)    
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionNeptune, GM_Neptune_DE440)
    a <- a + pointMassAcceleration(Y[1:3], JPL_ephemerides$positionPluto, GM_Pluto_DE440)
    dY <- matrix(c(Y[4:6], a), byrow=TRUE, ncol=3, nrow=2)
    colnames(dY) <- c("X", "Y", "Z")
    rownames(dY) <- c("Velocity", "Acceleration")
    return(list(dY, realCentralBody))
}
