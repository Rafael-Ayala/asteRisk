sgp4 <- function(n0, e0, i0, M0, omega0, OMEGA0, Bstar, initialDateTime=NULL, targetTime,
                 keplerAccuracy=10e-12, maxKeplerIterations=10) {
    checkTimeInput(initialDateTime, targetTime, "sgp4")
    if(2*pi/n0 >= 225) {
        deep_space <- TRUE
    }
    t0 <- 0
    if(is.character(initialDateTime) & is.character(targetTime)) {
        t <- as.numeric(difftime(targetTime, initialDateTime, units="secs"))/60
    } else {
        t <- targetTime
    }
    a1 <- (ke/n0)^(2/3) # Note: ke is the inverse of tumin in Vallados implementation
    delta1 <- 1.5 * (k2/(a1^2)) * (((3*(cos(i0))^2) - 1)/((1-e0^2)^(3/2)))
    a0 <- a1 * (1 - (1/3)*delta1 - delta1^2 - (134/81) * delta1^3)
    delta0 <- 1.5 * (k2/(a0^2)) * (((3*(cos(i0))^2) - 1)/((1-e0^2)^(3/2))) #TODO replace common term of delta1 and delta0 by variable definition, divide in each case by a1 and a0
    n0dprime <- n0/(1+delta0) # un_kozai mean motion, called no_unkozai
    #a0dprime <- a0/(1-delta0)
    # use no_unkozai semi major axis
    a0dprime <- (ke/n0dprime)^(2/3) # called ao in Vallado's implementation
    perigee <- (a0dprime*(1-e0)-ae)*earthRadius_SGP4 # equivalent to rp in Vallado's implementation, subtracting 1 and multiplying by Earth Radius
    high_perigee_flag <- TRUE
    if(perigee >= 98 & perigee <= 156) {
        s <- a0dprime*(1-e0) - s + ae # equivalent to: (perigee - 78)/6378.135 + 1
    } else if(perigee < 98) {
        s <- 20/earthRadius_SGP4 + ae # equivalent to 20/6378.135 + 1
    }
    if(perigee < 220) {
        high_perigee_flag <- FALSE # Vallado's implementation uses isimp, which is a low perigee flag (therefore values are inverted)
    }
    theta <- cos(i0)
    xi <- 1/(a0dprime - s) # called tsi in Vallado's implementation
    beta0 <- sqrt(1-e0^2)
    eta <- a0dprime * e0 * xi
    C2 <-(q0 - s)^4 * xi^4 * n0dprime * (abs(1 - eta^2))^(-7/2) * (a0dprime * (1 + 1.5*eta^2 + 4*e0*eta + e0*eta^3) + 1.5 * ((k2*xi)/abs(1 - eta^2)) * (-0.5 * 1.5*theta^2) * (8+ 24*eta^2+ 3*eta^4))
    C1 <- Bstar * C2
    # Introduction of C3 = 0 if eccentricity lower than 1e-4 to match Vallado's implementation
    if(e0 > 1e-4) {
        C3 <- (q0 - s)^4 * xi^5 * A30 * n0dprime * ae * sin(i0) *2/J2 *(1/e0)
    } else {
        C3 <- 0
    }
    # C3 <- (q0 - s)^4 * xi^5 * A30 * n0dprime * ae * sin(i0)
    C4 <- 2 * n0dprime * (q0 - s)^4 * xi^4 * a0dprime * beta0^2 * (1-eta^2)^(-7/2) * ( (2*eta*(1+e0*eta) + 0.5*e0 + 0.5*eta^3) - ((2*k2*xi)/(a0dprime*abs(1-eta^2))) *  (3*(1-3*theta^2) * (1+1.5*eta^2-2*e0*eta-0.5*e0*eta^3) + 0.75*(1-theta^2)*(2*eta^2 - e0*eta - e0*eta^3)*cos(2*omega0)) )
    C5 <- 2*(q0 - s)^4*xi^4*a0dprime*beta0^2*(1-eta^2)^(-7/2)*(1+(11/4)*eta*(eta+e0)+e0*eta^3)
    D2 <- 4*a0dprime * xi * C1^2
    D3 <- (4/3) * a0dprime * xi^2 *(17*a0dprime + s)*C1^3
    D4 <- (2/3) * a0dprime * xi^3 * (221*a0dprime + 31*s)*C1^4
    MDF <- M0 + ( 1 + ( (3*k2*(-1+3*theta^2) ) / (2*a0dprime^2*beta0^3) ) + ( (3*k2^2*(13-78*theta^2+137*theta^4)) / (16*a0dprime^4*beta0^7) ) ) * n0dprime*(t - t0)
    omegaDF <- omega0 + ( ( -(3*k2*(1-5*theta^2)) / (2*a0dprime^2*beta0^4) ) + ( (3*k2^2*(7-114*theta^2+395*theta^4)) / (16*a0dprime^4*beta0^8) ) + ( (5*k4*(3-36*theta^2+49*theta^4)) / (4*a0dprime^4*beta0^8) ) )*n0dprime*(t - t0)
    OMEGADF <- OMEGA0 + ( ( (-3*k2*theta) / (a0dprime^2*beta0^4) ) + ( (3*k2^2*(4*theta-19*theta^3)) / (2*a0dprime^4*beta0^8) ) + ( (5*k4*theta*(3-7*theta^2)) / (2*a0dprime^4*beta0^8) ) ) * n0dprime * (t - t0)
    deltaomega <-Bstar*C3*(cos(omega0))*(t - t0)
    # Introduction of deltaM = 0 if eccentricity lower than 1e-4 to match Vallado's implementation
    if(e0 > 1e-4) {
        deltaM <- (-2/3) * (q0 - s)^4*Bstar*xi^4*(ae/(e0*eta)) * ( (1 + eta*cos(MDF))^3 - (1 + eta*cos(M0))^3 )
    } else {
        deltaM <- 0
    }
    # Mp <- MDF + deltaomega + deltaM # Changed to match Vallado's implementation 
    Mp <- MDF + high_perigee_flag*deltaomega + high_perigee_flag*deltaM 
    omega <- omegaDF - deltaomega - deltaM
    OMEGA <- OMEGADF - (21/2)*( (n0dprime*k2*theta) / (a0dprime^2 * beta0^2) )*C1*(t - t0)^2
    e <- e0 - Bstar*C4*(t-t0) - Bstar*C5*(sin(Mp) - sin(M0))
    a <- a0dprime*((1 - C1*(t - t0) - high_perigee_flag*D2*(t - t0)^2 - high_perigee_flag*D3*(t - t0)^3 - high_perigee_flag*D4 * (t - t0)^4)^2)
    IL <- Mp + omega + OMEGA + n0dprime * (1.5*C1*(t - t0)^2 + high_perigee_flag*(D2 + 2*C1^2)*(t-t0)^3 + high_perigee_flag*0.25*(3*D3+12*C1*D2+10*C1^3)*(t-t0)^4 + high_perigee_flag*0.2*(3*D4 + 12*C1*D3 + 6*D2^2 + 30*C1^2*D2 + 15*C1^4)*(t-t0)^5)
    beta <- sqrt(1-e^2)
    n <- ke/(a^(1.5))
    axN <- e*cos(omega)
    ILL <- ( (A30*sin(i0)) / (8*k2*a*beta^2) ) * e * cos(omega) * ( (3 + 5* theta) / (1 + theta) )
    ayNL <- (A30*sin(i0))/(4*k2*a*beta^2)
    ILT <- IL + ILL
    # ILT <- rem(IL + ILL, 2*pi)
    ayN <- e * sin(omega) + ayNL
    U <- rem(ILT - OMEGA, 2*pi)
    kepler_sol_1 <- U
    kepler_sol_current <- kepler_sol_1
    convergence <- FALSE
    iterations <- 0
    while(!convergence) {
        iterations <- iterations + 1
        delta_kepler_sol <- ( (U-ayN*cos(kepler_sol_current) + axN*sin(kepler_sol_current) - kepler_sol_current) / (-ayN*sin(kepler_sol_current) - axN*cos(kepler_sol_current) + 1) )
        kepler_sol_current <- kepler_sol_current + delta_kepler_sol
        if(iterations > maxKeplerIterations | abs(delta_kepler_sol) < keplerAccuracy) {
            convergence <- TRUE
        }
    }
    kepler_sol <- kepler_sol_current
    ecosE <- axN*cos(kepler_sol) + ayN*sin(kepler_sol)
    esinE <- axN*sin(kepler_sol) - ayN*cos(kepler_sol)
    eL <- sqrt((axN^2 + ayN^2))
    pL <- a*(1 - eL^2)
    r <- a*(1 - ecosE)
    # rderivative <- ke*(sqrt(a)/r)*esinE # Changed to match Vallado's implementation
    rderivative <- (sqrt(a)/r)*esinE
    # rfderivative <- ke*sqrt(pL)/r # Changed to match Vallado's implementation
    rfderivative <- sqrt(pL)/r
    cosu <- (a/r) * (cos(kepler_sol) - axN + ( (ayN*esinE) / (1 + sqrt(1 - eL^2)) ))
    sinu <- (a/r)*(sin(kepler_sol) - ayN - ( (axN*esinE) / (1+ sqrt(1 - eL^2)) ))
    u <- atan2(sinu,cosu)
    Deltar <- (k2/(2*pL)) * (1-theta^2) * cos(2*u)
    Deltau <- (-k2/(4*pL^2)) * (7*theta^2 - 1) * sin(2*u)
    DeltaOMEGA <- ( (3*k2*theta) / (2*pL^2) ) * sin(2*u)
    Deltai <- ( (3*k2*theta) / (2*pL^2) ) * sin(i0) * cos(2*u)
    # Deltarderivative <- ( (-k2*n) / pL ) * (1-theta^2) *sin(2*u) # Changed to match Vallado's implementation
    Deltarderivative <- ( (-k2*n) / pL ) * (1-theta^2) *sin(2*u)/ke
    # Deltarfderivative <- ( (k2*n) / pL ) * ( (1-theta^2) * cos(2*u) - 1.5*(1-3*theta^2) ) # Changed to match Vallado's implementation
    Deltarfderivative <- ( (k2*n) / pL ) * ( (1-theta^2) * cos(2*u) - 1.5*(1-3*theta^2) )/ke
    rk <- r * (1 - 1.5*k2*( ( sqrt(1-eL^2) ) / ( pL^2 ) ) * (3*theta^2-1) ) + Deltar
    uk <- u + Deltau
    OMEGAk <- OMEGA + DeltaOMEGA
    ik <- i0 + Deltai
    rkderivative <- rderivative + Deltarderivative
    rfkderivative <- rfderivative + Deltarfderivative
    Mx <- -sin(OMEGAk) * cos(ik)
    My <- cos(OMEGAk) * cos(ik)
    Mz <- sin(ik)
    Mvector <- c(Mx, My, Mz)
    Nx <- cos(OMEGAk)
    Ny <- sin(OMEGAk)
    Nz <- 0
    Nvector <- c(Nx, Ny, Nz)
    Uvector <- Mvector*sin(uk) + Nvector*cos(uk)
    Vvector <- Mvector*cos(uk) - Nvector*sin(uk)
    rvector <- rk*Uvector
    # rderivativevector <- rkderivative*Uvector -rfkderivative * Vvector # Changed to match Vallado's implementation
    rderivativevector <- rkderivative*Uvector + rfkderivative * Vvector
    position_result <- rvector*earthRadius_SGP4
    # velocity_result <- rderivativevector*earthRadius_SGP4*1440/86400 # Changed to match Vallado's implementation
    velocity_result <- rderivativevector*earthRadius_SGP4*ke/60
    return(list(
        position=position_result,
        velocity=velocity_result,
        algorithm="sgp4"))
}

sdp4 <- function(n0, e0, i0, M0, omega0, OMEGA0, Bstar, initialDateTime, targetTime,
                 keplerAccuracy=10e-8, maxKeplerIterations=100) {
    checkTimeInput(initialDateTime, targetTime, "sdp4")
    t0 <- 0
    if(is.character(initialDateTime) & is.character(targetTime)) {
        t <- as.numeric(difftime(targetTime, initialDateTime, units="secs"))/60
    } else {
        t <- targetTime
    }
    
    # Initialize auxiliary variables
    atime <- xli <- xni <- xfact <- ssl <- ssg <- ssh <- sse <- ssi <- xlamo <-
        gmst <- del1 <- del2 <- del3 <- fasx2 <- fasx4 <- fasx6 <- d2201 <- 
        d2211 <- d3210 <- d3222 <- d4410 <- d4422 <- d5220 <- d5232 <- d5421 <-
        d5433 <- xnddt <- xndot <- xldot <- zmos <- se2 <- se3 <- si2 <- si3 <-
        sl2 <- sl3 <- sl4 <- sgh2 <- sgh3 <- sgh4 <- sh2 <- sh3 <- zmol <-
        ee2 <- e3 <- xi2 <- xi3 <- xl2 <- xl3 <- xl4 <- xgh2 <- xgh3 <- xgh4 <-
        xh2 <- xh3 <- pe <- pinc <- pgh <- ph <- pl <- pgh0 <- ph0 <- pe0 <-
        pinc0 <- pl0 <- se <- si <- sl <- sgh <- shdq <- 0.0
    
    ## start equivalent of sgp4_init
    a1 <- (ke/n0)^(2/3)
    delta1 <- 1.5 * (k2/(a1^2)) * (((3*(cos(i0))^2) - 1)/((1-e0^2)^(3/2)))
    a0 <- a1 * (1 - (1/3)*delta1 - delta1^2 - (134/81) * delta1^3)
    delta0 <- 1.5 * (k2/(a0^2)) * (((3*(cos(i0))^2) - 1)/((1-e0^2)^(3/2)))
    n0dprime <- n0/(1+delta0)
    #a0dprime <- a0/(1-delta0)
    # use no_unkozai semi major axis
    a0dprime <- (ke/n0dprime)^(2/3)
    perigee <- (a0dprime*(1-e0)-ae)*earthRadius_SGP4
    if (perigee < 156) {
        if (perigee < 98) {
            s <- 20/earthRadius_SGP4 + ae
        } else {
            s <- a0dprime*(1-e0) - s + ae
        }
    }
    theta <- cos(i0)
    xi <- 1/(a0dprime - s)
    beta0 <- sqrt(1-e0^2)
    eta <- a0dprime * e0 * xi

    aux0 <- abs(1-eta^2)
    aux1 <- aux0^(-7/2)
    aux2 <- xi^4*a0dprime*beta0^2*aux1

    C2 <- (q0 - s)^4 * xi^4 * n0dprime * aux1 * (a0dprime * (1 + 1.5*eta^2 + 4*e0*eta + e0*eta^3) + 1.5 * ((k2*xi)/aux0) * (-0.5 * 1.5*theta^2) * (8+ 24*eta^2+ 3*eta^4))

    C1 <- Bstar * C2

    C4 <- 2 * n0dprime * (q0 - s)^4 * aux2 * ( (2*eta*(1+e0*eta) + 0.5*e0 + 0.5*eta^3) - ((2*k2*xi)/(a0dprime*aux0)) *  (3*(1-3*theta^2) * (1+1.5*eta^2-2*e0*eta-0.5*e0*eta^3) + 0.75*(1-theta^2)*(2*eta^2 - e0*eta - e0*eta^3)*cos(2*omega0)) )

    Mdot <- ( 1 + ( (3*k2*(-1+3*theta^2) ) / (2*a0dprime^2*beta0^3) ) + ( (3*k2^2*(13-78*theta^2+137*theta^4)) / (16*a0dprime^4*beta0^7) ) ) * n0dprime

    omegadot <- ( ( -(3*k2*(1-5*theta^2)) / (2*a0dprime^2*beta0^4) ) + ( (3*k2^2*(7-114*theta^2+395*theta^4)) / (16*a0dprime^4*beta0^8) ) + ( (5*k4*(3-36*theta^2+49*theta^4)) / (4*a0dprime^4*beta0^8) ) )*n0dprime

    OMEGA1dot <- ( (-3*k2*theta) / (a0dprime^2*beta0^4) )*n0dprime

    OMEGAdot <- OMEGA1dot + (( (3*k2^2*(4*theta-19*theta^3)) / (2*a0dprime^4*beta0^8) ) + ( (5*k4*theta*(3-7*theta^2)) / (2*a0dprime^4*beta0^8) ))*n0dprime

    # Proper deep space initialization algorithm, equivalent to dsinit. Still inside
    # sgp4_init

    e0_2 <- e0^2
    e0_3 <- e0^3
    sqrt_1_e0_2 <- sqrt(1-e0_2)
    inv_a0dprime <- 1/a0dprime
    inv_n0dprime <- 1/n0dprime
    sin_i0 <- sin(i0)
    cos_i0 <- cos(i0)
    sin_OMEGA0 <- sin(OMEGA0)
    cos_OMEGA0 <- cos(OMEGA0)
    sin_omega0 <- sin(omega0)
    cos_omega0 <- cos(omega0)
    sin_i0_2 <- sin_i0^2
    cos_i0_2 <- cos_i0^2

    xpidot <- omegadot + OMEGAdot

    ishq <- if (i0 >= 3*pi/180) TRUE else FALSE

    # Avoid sin(i0) = 0
    if(abs(sin(i0)) < 1e-12 ) {
        sign_i0 <- sign(i0)
        i0 <- sign_i0 * 1e-12
        sin_i0 <- sin(i0)
    }
    # get GMST at TLE epoch
    gmst <- UTCdateTimeToGMST(initialDateTime)
    day <- as.numeric(difftime(initialDateTime, "1900-01-01 12:00:00 UTC", tz="UTC")) + 1

    xnodce <- (4.5236020 - 9.2422029e-4*day) %% (2*pi)
    stem <- sin(xnodce)
    ctem <- cos(xnodce)
    zcosil <- 0.91375164 - 0.03568096*ctem
    zsinil <- sqrt(1 - zcosil^2)
    zsinhl <- 0.089683511*(stem/zsinil)
    zcoshl <- sqrt(1 - zsinhl^2)
    gam <- 5.8351514 + 0.0019443680*day

    zx <- 0.39785416*(stem/zsinil)
    zy <- zcoshl*ctem + 0.91744867*zsinhl*stem
    zx <- atan2(zx, zy)
    zx <- gam + zx - xnodce
    zsingl <- sin(zx)
    zcosgl <- cos(zx)
    zmol <- (4.7199672 + 0.22997150*day - gam) %% (2*pi)
    zmos <- (6.2565837 + 0.017201977*day) %% (2*pi)

    # Solar terms
    zcosg <- ZCOSGS
    zsing <- ZSINGS
    zcosi <- ZCOSIS
    zsini <- ZSINIS
    zcosh <- cos_OMEGA0
    zsinh <- sin_OMEGA0
    cc <- C1SS
    zn <- ZNS
    ze <- ZES
    zmo <- zmos
    lunarSolarDone <- FALSE
    while(TRUE) {
        a1  <- zcosg*zcosh + zsing*zcosi*zsinh
        a3  <- -zsing*zcosh + zcosg*zcosi*zsinh
        a7  <- -zcosg*zsinh + zsing*zcosi*zcosh
        a8  <- zsing*zsini
        a9  <- zsing*zsinh + zcosg*zcosi*zcosh
        a10 <- zcosg*zsini
        a2  <- cos_i0*a7 + sin_i0*a8
        a4  <- cos_i0*a9 + sin_i0*a10
        a5  <- -sin_i0*a7 + cos_i0*a8
        a6  <- -sin_i0*a9 + cos_i0*a10

        x1 <- a1*cos_omega0 + a2*sin_omega0
        x2 <- a3*cos_omega0 + a4*sin_omega0
        x3 <- -a1*sin_omega0 + a2*cos_omega0
        x4 <- -a3*sin_omega0 + a4*cos_omega0
        x5 <- a5*sin_omega0
        x6 <- a6*sin_omega0
        x7 <- a5*cos_omega0
        x8 <- a6*cos_omega0

        z31 <- 12*x1^2 - 3*x3^2
        z32 <- 24*x1 * x2 - 6*x3 * x4
        z33 <- 12*x2^2 - 3*x4^2
        z1 <- 3*(a1^2 + a2^2) + z31 * e0_2
        z2 <- 6*(a1 * a3 + a2 * a4) + z32 * e0_2
        z3 <- 3*(a3^2 + a4^2) + z33 * e0_2
        z11 <- -6*a1 * a5 + e0_2 * (-24*x1 * x7 - 6*x3 * x5)
        z12 <- -6*(a1 * a6 + a3 * a5) + e0_2 * (-24*(x2 * x7 + x1 * x8) - 6*(x3 * x6 + x4 * x5) )
        z13 <- -6*a3 * a6 + e0_2 * (-24*x2 * x8 - 6*x4 * x6)
        z21 <- 6*a2 * a5 + e0_2 * (+24*x1 * x5 - 6*x3 * x7)
        z22 <- 6*(a4 * a5 + a2 * a6) + e0_2 * (24*(x2 * x5 + x1 * x6) - 6*(x4 * x7 + x3 * x8))
        z23 <- 6*a4 * a6 + e0_2 * (24*x2 * x6 - 6*x4 * x8)
        z1 <- 2*z1 + (1 - e0_2) * z31
        z2 <- 2*z2 + (1 - e0_2) * z32
        z3 <- 2*z3 + (1 - e0_2) * z33
        s3 <- cc * inv_n0dprime
        s2 <- -0.5*s3 / sqrt_1_e0_2
        s4 <- s3 * sqrt_1_e0_2
        s1 <- -15*e0 * s4
        s5 <- x1 * x3 + x2 * x4
        s6 <- x2 * x3 + x1 * x4
        s7 <- x2 * x4 - x1 * x3
        se <- s1 * zn * s5
        si <- s2 * zn * (z11 + z13)
        sl <- -zn * s3 * (z1 + z3 - 14 - 6*e0_2)
        sgh <- s4 * zn * (z31 + z33 - 6)
        if(ishq) {
            sh <- -zn * s2 * (z21 + z23)
            shdq = sh / sin_i0
        }
        ee2 <- 2*s1 * s6
        e3 <- 2*s1 * s7
        xi2 <- 2*s2 * z12
        xi3 <- 2*s2 * (z13 - z11)
        xl2 <- -2*s3 * z2
        xl3 <- -2*s3 * (z3 - z1)
        xl4  <- -2*s3 * (-21 - 9*e0_2) * ze
        xgh2 <- 2*s4 * z32
        xgh3 <- 2*s4 * (z33 - z31)
        xgh4 <- -18*s4 * ze
        xh2 <- -2*s2 * z22
        xh3 <- -2*s2 * (z23 - z21)

        if(lunarSolarDone) break

        ## Lunar terms
        sse <- se
        ssi <- si
        ssl <- sl
        ssh <- shdq
        ssg <- sgh - cos_i0 * ssh
        se2 <- ee2
        si2 <- xi2
        sl2 <- xl2
        sgh2 <- xgh2
        sh2 <- xh2
        se3 <- e3
        si3 <- xi3
        sl3 <- xl3
        sgh3 <- xgh3
        sh3 <- xh3
        sl4 <- xl4
        sgh4 <- xgh4
        zcosg <- zcosgl
        zsing <- zsingl
        zcosi <- zcosil
        zsini <- zsinil
        zcosh <- cos_OMEGA0 * zcoshl + sin_OMEGA0 * zsinhl
        zsinh = sin_OMEGA0 * zcoshl - cos_OMEGA0 * zsinhl
        zn <- ZNL
        cc <- C1L
        ze <- ZEL
        zmo <- zmol
        lunarSolarDone <- TRUE
    }
    sse <- sse + se
    ssi <- ssi + si
    ssl <- ssl + sl
    ssg <- ssg + sgh - cos_i0*shdq
    ssh <- ssh + shdq
    if(n0dprime < 0.0052359877 & n0dprime > 0.0034906585) {
        # object is in a 1-day resonant period regime
        iresfl <- TRUE
        isynfl <- TRUE
        g200 <- e0_2 * (0.8125*e0_2 - 2.5) + 1
        g310 <- 2*e0_2 + 1
        g300 <- e0_2 * (6.60937*e0_2 - 6) + 1
        f220 <- 0.75*(cos_i0 + 1)^2
        f311 <- 0.9375*(3.0*cos_i0 + 1)*sin_i0^2 - 0.75 * (cos_i0 + 1)
        f330 <- 1.875 * (cos_i0 + 1)^3
        del1 <- 3*(n0dprime^2 * inv_a0dprime^2)
        del2 <- 2*del1 * f220 * g200 * Q22
        del3 <- 3*del1 * f330 * g300 * Q33 * inv_a0dprime
        del1 <-  del1 * f311 * g310 * Q31 * inv_a0dprime
        fasx2 <- 0.13130908
        fasx4 <- 2.8843198
        fasx6 <- 0.37448087
        xlamo <- (M0 + OMEGA0 + omega0 - gmst) %% (2*pi)
        bfact <- Mdot + xpidot - THDT + ssl + ssg + ssh
    } else if (n0dprime >0.00826 & n0dprime <= 0.00924 & e0 >= 0.5) {
        # object is in a half-day resonant period regime
        iresfl <- TRUE
        isynfl <- FALSE
        g201 <- -0.306 - 0.44*(e0 - 0.64)
        if (e0 <= 0.65) {
            g211 <- 3.6160 - 13.2470*e0 + 16.29000*e0_2
            g310 <- - 19.3020 + 117.3900*e0 - 228.4190*e0_2 + 156.5910*e0_3
            g322 <- -18.9068 +109.7927*e0 -214.6334*e0_2 +146.5816*e0_3
            g410 <- -41.1220 +242.6940*e0 -471.0940*e0_2 +313.9530*e0_3
            g422 <- -146.407 +841.8800*e0 -1629.014*e0_2 +1083.435*e0_3
            g520 <- -532.114 +3017.977*e0 -5740.032*e0_2 +3708.276*e0_3
        } else {
            g211 <- -72.099 +331.8190*e0 -508.7380*e0_2 +266.7240*e0_3
            g310 <- -346.844 +1582.851*e0 -2415.925*e0_2 +1246.113*e0_3
            g322 <- -342.585 +1554.908*e0 -2366.899*e0_2 +1215.972*e0_3
            g410 <- -1052.797 +4758.686*e0 -7193.992*e0_2 +3651.957*e0_3
            g422 <- -3581.690 +16178.11*e0 -24462.77*e0_2 +12422.52*e0_3
            if (e0 <= 0.715) {
                g520 <- +1464.74 -4664.75*e0 +3763.64*e0_2
            }
            else {
                g520 <- -5149.66 +29936.92*e0 -54087.36*e0_2 +31324.56*e0_3
            }
        }
        if (e0 < 0.7) {
            g533 <- -919.22770 +4988.6100*e0 -9064.7700*e0_2 +5542.210*e0_3
            g521 <- -822.71072 +4568.6173*e0 -8491.4146*e0_2 +5337.524*e0_3
            g532 <- -853.66600 +4690.2500*e0 -8624.7700*e0_2 +5341.400*e0_3
        } else {
            g533 <- -37995.780 +161616.52*e0 -229838.20*e0_2 +109377.94*e0_3
            g521 <- -51752.104 +218913.95*e0 -309468.16*e0_2 +146349.42*e0_3
            g532 <- -40023.880 +170470.89*e0 -242699.48*e0_2 +115605.82*e0_3
        }
        f220 <- 0.75*(1 + 2*cos_i0 + cos_i0_2)
        f221 <- 1.5*sin_i0_2
        f321 <- 1.875*sin_i0 * (1 - 2*cos_i0 - 3*cos_i0_2)
        f322 <- -1.875*sin_i0 * (1 + 2*cos_i0 - 3*cos_i0_2)
        f441 <- 35*sin_i0_2 * f220
        f442 <- 39.375*sin_i0_2^2
        f522 <- 9.84375*sin_i0 * (sin_i0_2 * (1 - 2*cos_i0 - 5*cos_i0_2) +
                                      (1/3)*(-2 + 4*cos_i0 + 6*cos_i0_2) )
        f523 <- sin_i0 * (4.92187512*sin_i0_2 * (-2 - 4*cos_i0 + 10*cos_i0_2) +
                              6.56250012*(1 + 2*cos_i0 -  3*cos_i0_2) )
        f542 <- 29.53125*sin_i0 * (2 - 8*cos_i0 + cos_i0_2 * (-12 + 8*cos_i0 + 10*cos_i0_2))
        f543 <- 29.53125*sin_i0 * (-2 - 8*cos_i0 + cos_i0_2 * (12 + 8*cos_i0 - 10*cos_i0_2))
        temp1 <- 3*(n0dprime * inv_a0dprime)^2
        temp0 <- temp1 * ROOT22
        d2201 <- temp0 * f220 * g201
        d2211 <- temp0 * f221 * g211
        temp1 <- temp1 * inv_a0dprime
        temp0 <- temp1 * ROOT32
        d3210 <- temp0 * f321 * g310
        d3222 <- temp0 * f322 * g322
        temp1 <- temp1 * inv_a0dprime
        temp0 <- 2*temp1 * ROOT44
        d4410 <- temp0 * f441 * g410
        d4422 <- temp0 * f442 * g422
        temp1 <- temp1 * inv_a0dprime
        temp0 <- temp1 * ROOT52
        d5220 <- temp0 * f522 * g520
        d5232 <- temp0 * f523 * g532
        temp0 <- 2*temp1 * ROOT54
        d5421 <- temp0 * f542 * g521
        d5433 <- temp0 * f543 * g533
        xlamo <- (M0 + 2*OMEGA0 - 2*gmst) %% (2*pi)
        bfact <- Mdot + 2*OMEGAdot - 2 *THDT + ssl + 2*ssh
    } else {
        # object is not in a resonant period regime
        iresfl <- FALSE
        isynfl <- FALSE
    }
    # Prepare integrator
    if (iresfl) {
        xfact <- bfact - n0dprime
        xli <- xlamo
        atime = 0.0
        xni <- n0dprime #potentially unnecessary
    }
    # Compute dot (derivative) terms
    if (isynfl) {
        sin_1 <- sin(xli-fasx2)
        cos_1 <- cos(xli-fasx2)
        sin_2 <- sin(2*(xli - fasx4))
        cos_2 <- cos(2*(xli - fasx4))
        sin_3 <- sin(3*(xli - fasx6))
        cos_3 <- cos(3*(xli - fasx6))
        xndot <- del1 * sin_1 + del2 * sin_2 + del3 * sin_3
        xnddt <- del1 * cos_1 + 2*del2 * cos_2 + 3*del3 * cos_3
    } else if(iresfl) {
        omega <- omega0 + omegadot + atime
        sin_1 <- sin(2*omega + xli - G22)
        cos_1 <- cos(2*omega + xli - G22)
        sin_2 <- sin(xli  - G22)
        cos_2 <- cos(xli  - G22)
        sin_3 <- sin(omega + xli  - G32)
        cos_3 <- cos(omega + xli  - G32)
        sin_4 <- sin(-omega + xli  - G32)
        cos_4 <- cos(-omega + xli  - G32)
        sin_5 <- sin(omega + xli  - G52)
        cos_5 <- cos(omega + xli  - G52)
        sin_6 <- sin(-omega + xli  - G52)
        cos_6 <- cos(-omega + xli  - G52)
        sin_7 <- sin(2*omega + 2*xli - G44)
        cos_7 <- cos(2*omega + 2*xli - G44)
        sin_8 <- sin(2*xli - G44)
        cos_8 <- cos(2*xli - G44)
        sin_9 <- sin(omega + 2*xli - G54)
        cos_9 <- cos(omega + 2*xli - G54)
        sin_10 <- sin(-omega + 2*xli - G54)
        cos_10 <- cos(-omega + 2*xli - G54)
        xndot <- d2201 * sin_1 + d2211 * sin_2 + d3210 * sin_3 +
            d3222 * sin_4 + d5220 * sin_5 + d5232 * sin_6 +
            d4410 * sin_7 + d4422 * sin_8 + d5421 * sin_9 +
            d5433 * sin_10
        xnddt <- d2201 * cos_1 + d2211 * cos_2 + d3210 * cos_3 +
            d3222 * cos_4 + d5220 * cos_5 + d5232 * cos_6 +
            2*(d4410 * cos_7 + d4422 * cos_8 + d5421 * cos_9 +
                  d5433 * cos_10)
    }
    xldot <- xni + xfact
    xnddt <- xnddt * xldot

    nk <- n0dprime
    ak <- a0dprime
    ek <- e0
    ik <- i0
    OMEGAk <- OMEGA0
    omegak <- omega0
    Mk <- M0

    ## Introduce secular effects of atmospheric drag (dependent on perigee altitude)
    ## and gravity.
    ## Mk, OMEGAk and omegak equivalent to MDF, OMEGA and omega
    Mk <- M0 + Mdot*(t-t0)
    OMEGAk <- OMEGA0 + OMEGAdot*(t-t0) - (21/2)*( (n0dprime*k2*theta) / (a0dprime^2*beta0^2) )*C1*(t-t0)^2
    omegak <- omega0 + omegadot*(t-t0)

    ## introduction of deep space secular effects
    Msec <- Mk + ssl*(t-t0)
    esec <- e0 + sse*(t-t0)
    isec <- i0 + ssi*(t-t0)
    OMEGAsec <- OMEGAk + ssh*(t-t0)
    omegasec <- omegak + ssg*(t-t0)
    # clarify requirement for theta_ds. value equal to original theta
    theta_ds <- (gmst + THDT*(t-t0)) %% (2*pi)
    if(iresfl) {
        if ((atime == 0) | ((t-t0) * atime <= 0) | (abs(t-t0) < abs(atime))) {
            atime <- 0
            xni <- n0dprime
            xli <- xlamo
        }
        ft <- (t - t0) - atime
        delt <- if((t - t0) >= atime) STEP else -STEP
        while(TRUE) {
            if(isynfl) {
                sin_1 <- sin(xli-fasx2)
                cos_1 <- cos(xli-fasx2)
                sin_2 <- sin(2*(xli - fasx4))
                cos_2 <- cos(2*(xli - fasx4))
                sin_3 <- sin(3*(xli - fasx6))
                cos_3 <- cos(3*(xli - fasx6))
                xndot <- del1 * sin_1 + del2 * sin_2 + del3 * sin_3
                xnddt <- del1 * cos_1 + 2*del2 * cos_2 + 3*del3 * cos_3
            } else {
                omega <- omega0 + omegadot + atime
                sin_1 <- sin(2*omega + xli - G22)
                cos_1 <- cos(2*omega + xli - G22)
                sin_2 <- sin(xli  - G22)
                cos_2 <- cos(xli  - G22)
                sin_3 <- sin(omega + xli  - G32)
                cos_3 <- cos(omega + xli  - G32)
                sin_4 <- sin(-omega + xli  - G32)
                cos_4 <- cos(-omega + xli  - G32)
                sin_5 <- sin(omega + xli  - G52)
                cos_5 <- cos(omega + xli  - G52)
                sin_6 <- sin(-omega + xli  - G52)
                cos_6 <- cos(-omega + xli  - G52)
                sin_7 <- sin(2*omega + 2*xli - G44)
                cos_7 <- cos(2*omega + 2*xli - G44)
                sin_8 <- sin(2*xli - G44)
                cos_8 <- cos(2*xli - G44)
                sin_9 <- sin(omega + 2*xli - G54)
                cos_9 <- cos(omega + 2*xli - G54)
                sin_10 <- sin(-omega + 2*xli - G54)
                cos_10 <- cos(-omega + 2*xli - G54)

                xndot <- d2201 * sin_1 + d2211 * sin_2 + d3210 * sin_3 +
                    d3222 * sin_4 + d5220 * sin_5 + d5232 * sin_6 +
                    d4410 * sin_7 + d4422 * sin_8 + d5421 * sin_9 +
                    d5433 * sin_10

                xnddt <- d2201 * cos_1 + d2211 * cos_2 + d3210 * cos_3 +
                    d3222 * cos_4 + d5220 * cos_5 + d5232 * cos_6 +
                    2*(d4410 * cos_7 + d4422 * cos_8 + d5421 * cos_9 +
                           d5433 * cos_10)
            }
            xldot <- xni + xfact
            xnddt <- xnddt * xldot
            ft <- (t - t0) - atime
            if(abs(ft) < STEP) break
            xli <- xli + delt*(xldot + delt*xndot/2)
            xni <- xni + delt*(xndot + delt*xnddt/2)
            atime <- atime + delt
        }
        xl <- xli + ft*(xldot + ft * xndot/2)
        nsec <- xni + ft*(xndot + ft * xnddt/2)
        Msec <- if(!isynfl) (xl - 2*OMEGAsec + 2*theta_ds) else (xl - OMEGAsec - omegasec + theta_ds)
    }
    nk <- if(iresfl) nsec else n0dprime
    ek <- esec
    ik <- isec
    OMEGAk <- OMEGAsec
    omegak <- omegasec
    Mk <- Msec

    ak <- (ke/nk)^(2/3) * (1 - C1*(t - t0))^2
    ek <- ek - Bstar*C4*(t-t0)
    Mk <- Mk + n0dprime*(1.5*C1*(t-t0)^2)

    # normalizations
    Mk_aux <- Mk + omegak + OMEGAk  # related to IL
    OMEGAk <- rem(OMEGAk, 2*pi)
    omegak <- rem(omegak, 2*pi)
    Mk_aux <- rem(Mk_aux, 2*pi)
    Mk <- rem(Mk_aux - omegak - OMEGAk, 2*pi)

    ## calculate lunar-solar periodics effects
    ## update solar terms
    zm <- zmos + ZNS*(t-t0)
    zf <- zm + 2*ZES * sin(zm)
    sinzf <- sin(zf)
    coszf <- cos(zf)
    f2 <- sinzf * sinzf/2 - 0.25
    f3 <- -sinzf * coszf/2
    ses <- se2 * f2 + se3 * f3
    sis <- si2 * f2 + si3 * f3
    sls <- sl2 * f2 + sl3 * f3 + sl4 * sinzf
    sghs <- sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
    shs <- sh2  * f2 + sh3  * f3
    ## update lunar terms
    zm <- zmol + ZNL*(t-t0)
    zf <- zm   + 2*ZEL * sin(zm)
    sinzf <- sin(zf)
    coszf <- cos(zf)
    f2 <- sinzf * sinzf / 2 - 0.25
    f3 <- -sinzf * coszf / 2
    sel <- ee2 * f2 + e3 * f3
    sil <- xi2 * f2 + xi3 * f3
    sll <- xl2 * f2 + xl3 * f3 + xl4 * sinzf
    sghl <- xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
    shl <- xh2  * f2 + xh3  * f3
    # Store computed values
    pgh <- sghs + sghl
    ph <- shs  + shl
    pe <- ses  + sel
    pinc <- sis  + sil
    pl <- sls  + sll
    # Update inclination and eccentricity.
    e_per = ek + pe
    i_per = ik + pinc
    sinis <- sin(i_per)
    cosis <- cos(i_per)
    if (i_per >= 0.2) {
        tmp_ph <- ph / sinis
        omega_per <- omegak + pgh -cosis*tmp_ph
        OMEGA_per <- OMEGAk + tmp_ph
        M_per <- Mk + pl
    } else {
        sinok <- sin(OMEGAk)
        cosok <- cos(OMEGAk)
        # dalf
        alfdp <- sinis * sinok + ph * cosok + pinc * cosis * sinok
        # dbet
        betdp <- sinis * cosok - ph * sinok + pinc * cosis * cosok
        # reduce OMEGA per to 0-2pi interval since it is used outside of trigonometric function
        OMEGA_per <- OMEGAk %% (2*pi)
        # dls
        xls <- Mk + omegak + cosis * OMEGA_per + pl + pgh - pinc * OMEGA_per * sinis
        OMEGA_aux <- OMEGA_per
        OMEGA_per <- atan2(alfdp, betdp) %% (2*pi)
        if(abs(OMEGA_aux-OMEGA_per) > pi) {
            OMEGA_per <- if(OMEGA_per < OMEGA_aux) (OMEGA_per + 2*pi) else (OMEGA_per - 2*pi)
        }
        M_per <- Mk + pl
        omega_per <- xls - M_per - cosis * OMEGA_per
    }

    ek <- e_per
    ik <- i_per
    OMEGAk <- OMEGA_per
    omegak <- omega_per
    Mk <- M_per
    IL <- Mk + omegak + OMEGAk
    # enforce inclination positivity
    if(ik < 0) {
        ik <- -ik
        OMEGAk <- OMEGAk + pi
        omegak <- omegak - pi
    }
    # recalculate theta and sine values affected by updated inclination
    theta <- cos(ik)
    sin_ik <- sin(ik)
    if(ek < 1e-6) {
        ek <- 1e-6
    }
    beta <- sqrt(1-ek^2)
    nk <- ke/(ak^(3/2))
    sin_omegak <- sin(omegak)
    cos_omegak <- cos(omegak)
    axN <- ek*cos_omegak
    ayNL <- A30*sin_ik/(4*k2*ak*beta^2)
    ayN <- ek*sin_omegak + ayNL
    ILL <- 0.5*ayNL*axN*(3+5*theta)/(1+theta)
    ILT <- IL + ILL
    U <- rem(ILT - OMEGAk, 2*pi)
    Eomega <- U
    sinEomega <- 0
    cosEomega <- 0

    # Solve Kepler's equation
    convergence <- FALSE
    iterations <- 0
    while(!convergence) {
        iterations <- iterations + 1
        sinEomega <- sin(Eomega)
        cosEomega <- cos(Eomega)
        DeltaEomega <- (U-ayN*cosEomega+axN*sinEomega-Eomega)/(1-ayN*sinEomega-axN*cosEomega)
        sign_delta <- sign(DeltaEomega)
        if(abs(DeltaEomega) >= 0.95) {
            DeltaEomega <- 0.95 * sign_delta
        }
        Eomega <- Eomega + DeltaEomega
        if(iterations > maxKeplerIterations | abs(DeltaEomega) < keplerAccuracy) {
            convergence <- TRUE
        }
    }
    sinEomega <- sin(Eomega)
    cosEomega <- cos(Eomega)
    ecosE <- axN*cosEomega + ayN*sinEomega
    esinE <- axN*sinEomega - ayN*cosEomega
    eL <- sqrt((axN^2 + ayN^2))
    pL <- ak*(1 - eL^2)
    r <- ak*(1 - ecosE)
    rdot <- ke*(sqrt(ak)/r)*esinE
    rdotf <- ke*sqrt(pL)/r
    aux <- esinE/(1+sqrt(1-eL^2))
    cosu <- (ak/r) * (cosEomega - axN + ayN*aux)
    sinu <- (ak/r) * (sinEomega - ayN - axN*aux)
    cos2u <- 1 - 2*sinu^2
    sin2u <- 2*cosu*sinu
    u <- atan2(sinu, cosu)
    Deltar <- (k2/(2*pL)) * (1 - theta^2) * cos2u
    Deltau <- (-k2/(4*pL^2)) * (7*theta^2 - 1) * sin2u
    DeltaOMEGA <- ((3*k2*theta)/(2*pL^2)) * sin2u
    Deltai <- ((3*k2*theta)/(2*pL^2)) * sin_ik * cos2u
    Deltardot <- ((-k2*nk)/(pL))*(1-theta^2)*sin2u
    Deltardotf <- ( (k2*nk) / pL ) * ( (1-theta^2) * cos2u - 1.5*(1-3*theta^2) )
    rk <- r * (1 - 1.5*k2*( ( sqrt(1-eL^2) ) / ( pL^2 ) ) * (3*theta^2-1) ) + Deltar
    uk <- u + Deltau
    OMEGAk <- OMEGAk + DeltaOMEGA
    ik <- ik + Deltai
    rkdot <- rdot + Deltardot
    rkdotf <- rdotf + Deltardotf
    Mx <- -sin(OMEGAk) * cos(ik)
    My <- cos(OMEGAk) * cos(ik)
    Mz <-sin(ik)
    Mvector <- c(Mx, My, Mz)
    Nx <- cos(OMEGAk)
    Ny <- sin(OMEGAk)
    Nz <- 0
    Nvector <- c(Nx, Ny, Nz)
    Uvector <- Mvector*sin(uk) + Nvector*cos(uk)
    Vvector <- Mvector*cos(uk) - Nvector*sin(uk)
    position_result <- rk*Uvector*earthRadius_SGP4
    velocity_result <- (rkdot*Uvector+rkdotf*Vvector)*earthRadius_SGP4*1440/86400
    return(list(
        position=position_result,
        velocity=velocity_result,
        algorithm="sdp4"))
}

sgdp4 <- function(n0, e0, i0, M0, omega0, OMEGA0, Bstar, initialDateTime=NULL,
                            targetTime, keplerAccuracy=10e-8, maxKeplerIterations=100) {
    if(2*pi/n0 >= 225) {
        deep_space <- TRUE
    } else {
        deep_space <- FALSE
    }
    if(deep_space) {
        return(sdp4(n0, e0, i0, M0, omega0, OMEGA0, Bstar, initialDateTime,
                    targetTime, keplerAccuracy, maxKeplerIterations))
    } else {
        return(sgp4(n0, e0, i0, M0, omega0, OMEGA0, Bstar, initialDateTime,
                    targetTime, keplerAccuracy, maxKeplerIterations))
    }
}
