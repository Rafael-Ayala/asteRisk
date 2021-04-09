glatf <- function(lat) {
    dgtr <- 1.74533e-2
    c2 <- cos(2*dgtr*lat)
    gv <- 980.616 * (1 - 0.0026373 * c2)
    reff <- 2 * gv / (3.085462E-6 + 2.27E-9 * c2) * 1e-5
    return(list(
        gv = gv,
        reff = reff
    ))
}

ccor <- function(alt, r, h1, zh) {
    # chemistry/dissociation correction for msis models
    # ALT - altitude
    # R - target ratio
    # H1 - transition scale length
    # ZH - altitude of 1/2 R
    e <- (alt - zh) / h1
    if (e > 70) {
        corr <- 1
    } else if (e < -70) {
        corr <- exp(r)
    } else {
        ex <- exp(e)
        e <- r/(1 + ex)
        corr <- exp(e)
    }
    return(corr)
}

ccor2 <- function(alt, r, h1, zh, h2) {
    # chemistry/dissociation correction for msis models 2
    # ALT - altitude
    # R - target ratio
    # H1 - transition scale length
    # ZH - altitude of 1/2 R
    # H2 - transition scale length #2 ?
    e1 <- (alt - zh) / h1
    e2 <- (alt - zh) / h2
    if ((e1 > 70) | (e2 > 70)) {
        corr <- 1
    } else if ((e1 < -70) & (e2 < -70)) {
        corr <- exp(r)
    } else {
        ex1 <- exp(e1)
        ex2 <- exp(e2)
        ccor2v <- r / (1 + 0.5 * (ex1 + ex2))
        corr <- exp(ccor2v)
    }
    return(corr)
}

scalh <- function(alt, xm, temp, gsurf) {
    rgas <- 831.4
    g <- gsurf / ((1 + alt/re)^2)
    g <- rgas * temp / (g * xm)
    return(g)
}

dnet <- function(dd, dm, zhm, xmm, xm) {
    # turbopause correction for msis models
    # Root mean density
    # DD - diffusive density
    # DM - full mixed density
    # ZHM - transition scale length
    # XMM - full mixed molecular weight
    # XM  - species molecular weight
    # Outputs DNET - combined density
    a <- zhm / (xmm-xm)
    if (!((dm > 0) & (dd > 0))) {
        warning(paste("dnet log error", dm, dd, xm, sep=" "))
        if (dd == 0 & dm == 0) {
            dd <- 1
        }
        if (dm == 0) {
            dnet <- dd
        } else if (dd == 0) {
            dnet <- dm
        }
    }
    ylog <- a * log(dm/dd)
    if (ylog < -10) {
        dnet <- dd
    } else if (ylog > 10) {
        dnet <- dm
    } else {
        dnet <- dd * (1 + exp(ylog))^(1/a)
    }
    return(dnet)
}

zeta <- function(zz, zl) {
    result <- ((zz - zl) * (re + zl))/ (re + zz)
    return(result)
}

densm <- function(alt, d0, xm, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2) {
    # Calculate Temperature and Density Profiles for lower atmos
    #xs <- rep(0, 10)
    xs <- numeric(10)
    #ys <- rep(0, 10)
    ys <- numeric(10)
    densm_tmp <- d0
    if (alt > zn2[1]) {
        if(xm == 0){
            densm_tmp <- tz
        } else {
            densm_tmp <- d0
        }
    } else {
        z <- max(c(alt, zn2[mn2]))
        mn <- mn2
        z1 <- zn2[1]
        z2 <- zn2[mn]
        t1 <- tn2[1]
        t2 <- tn2[mn]
        zg <- zeta(z, z1)
        zgdif <- zeta(z2, z1)
        # set up spline nodes
        xs <- zeta(zn2[1:mn], z1)/zgdif
        ys <- 1/tn2[1:mn]
        xs <- xs[xs!=0]
        ys <- ys[ys!=0]
        # for (k in 1:mn) {
        #     xs[k] <- zeta(zn2[k], z1)/zgdif
        #     ys[k] <- 1/tn2[k]
        # }
        cubicspline <- splinefun(xs, ys)
        x <- zg/zgdif
        y <- cubicspline(x)
        tz <- 1/y
        if (xm != 0) {
            glb <- gsurf/((1 + z1/re)^2)
            gamm <- xm*glb*zgdif/rgas
            yi <- integrate(cubicspline, xs[1], x)
            expl <- gamm*yi
            expl <- min(expl, 50)
            densm_tmp <- densm_tmp * (t1/tz) * exp(-expl)
        }
        if (alt > zn3[1]) {
            if (xm == 0) {
                densm_tmp <- tz
            } 
        } else {
            z <- alt
            mn <- mn3
            z1 <- zn3[1]
            z2 <- zn3[mn]
            t1 <- tn3[1]
            t2 <- tn3[mn]
            zg <- zeta(z,z1)
            zgdif <- zeta(z2,z1)
            # set up spline nodes
            xs <- zeta(zn3[1:mn], z1)/zgdif
            ys <- 1/tn3[1:mn]
            xs <- xs[xs!=0]
            ys <- ys[ys!=0]
            cubicspline <- splinefun(xs, ys)
            x <- zg/zgdif
            y <- cubicspline(x)
            tz <- 1/y
            if (xm != 0) {
                glb <- gsurf/((1 + z1/re)^2)
                gamm <- xm*glb*zgdif/rgas
                yi <- integrate(cubicspline, xs[1], x)
                expl <- gamm*yi
                expl <- min(expl, 50)
                densm_tmp <- densm_tmp * (t1/tz) * exp(-expl)
            } else {
                densm_tmp <- tz
            }
        }
    }
    return(list(
        densm_tmp = densm_tmp,
        tz = tz
    ))
}

densu <- function(alt, dlb, tinf, tlb, xm, alpha, tz, zlb, s2, mn1, zn1, tn1, tgn1) {
    # Calculate Temperature and Density Profiles for MSIS models 
    # with new lower thermo polynomial
    xs <- numeric(10)
    ys <- numeric(10)
    za <- zn1[1]
    z <- max(alt, za)
    zg2 <- zeta(z, zlb)
    tt <- tinf - (tinf - tlb) * exp(-s2*zg2)
    ta <- tt
    tz <- tt
    densu_temp <- tz
    if (alt < za) { # calculate temperature profile below ZA
        dta <- (tinf - ta) * s2 * ((re + zlb)/(re + za))^2
        tgn1[1] <- dta
        tn1[1] <- ta
        z <- max(alt, zn1[mn1])
        mn <- mn1
        z1 <- zn1[1]
        z2 <- zn1[mn]
        t1 <- tn1[1]
        t2 <- tn1[mn]
        zg <- zeta(z, z1)
        zgdif <- zeta(z2, z1)
        xs <- zeta(zn1[1:mn], z1)/zgdif
        ys <- 1/tn1[1:mn]
        xs <- xs[xs!=0]
        ys <- ys[ys!=0]
        cubicspline <- splinefun(xs, ys)
        x <- zg/zgdif
        y <- cubicspline(x)
        tz <- 1/y
        densu_temp <- tz
    } else if (xm != 0) {
        glb <- gsurf/((1+zlb/re)^2)
        gamma <- xm * glb / (s2 * rgas * tinf)
        expl <- exp(-s2 * gamma * zg2)
        if(expl > 50 | tt <= 0) expl <- 50
        densa <- dlb * (tlb/tt)^(1+alpha+gamma) * expl
        densu_temp <- densa
        if (alt < za) {
            glb <- gsurf/((1+z1/re)^2)
            gamm <- xm * glb * zgdif / rgas
            yi <- integrate(cubicspline, xs[1], x)
            expl <- gamm * yi
            if(expl > 50 | tt <= 0) expl <- 50
            densu_temp <- densu_temp * (t1 / tz)^(1 + alpha) * exp(-expl)
        }
    }
    return(list(
        densu_temp = densu_temp,
        tz = tz))
}

g0 <- function(a, p) {
    result <- (a - 4 + (p[26] - 1) * (a - 4 + (exp(-abs(p[25]) * (a - 4)) - 1) / sqrt(p[25]*p[25])))
    return(result)
}

sumex <- function(ex) {
    result <- (1 + (1 - (ex^19)) / (1 - ex) * (ex^0.5))
    return(result)
}

sg0 <- function(ex, p, ap) {
    result <- (g0(ap[2], p) + (g0(ap[3], p)*ex + g0(ap[4], p)*ex*ex + 
                g0(ap[5], p)*(ex^3) + (g0(ap[6], p)*(ex^4) + 
                g0(ap[7], p)*(ex^12))*(1-(ex^8))/(1-ex)))/sumex(ex)
    return(result)
}

globe7 <- function(p, input, flags) {
    t <- numeric(15)
    sr <- 7.2722e-5
    dgtr <- 1.74533e-2
    dr <- 1.72142e-2
    hr <- 0.2618
    tloc <- input$lst
    c <- sin(input$g_lat * dgtr)
    s <- cos(input$g_lat * dgtr)
    c2 <- c^2
    c4 <- c^4
    s2 <- s^2
    plg[1,2] <<- c
    plg[1,3] <<- 0.5*(3*c2 -1)
    plg[1,4] <<- 0.5*(5*c*c2-3*c)
    plg[1,5] <<- (35*c4 - 30*c2 + 3)/8
    plg[1,6] <<- (63*c2*c2*c - 70*c2*c + 15*c)/8
    plg[1,7] <<- (11*c*plg[1,6] - 5*plg[1,5])/6
    plg[2,2] <<- s
    plg[2,3] <<- 3*c*s
    plg[2,4] <<- 1.5*(5*c2-1)*s
    plg[2,5] <<- 2.5*(7*c2*c-3*c)*s
    plg[2,6] <<- 1.875*(21*c4 - 14*c2 +1)*s
    plg[2,7] <<- (11*c*plg[2,6]-6*plg[2,5])/5
    plg[3,3] <<- 3*s2
    plg[3,4] <<- 15*s2*c
    plg[3,5] <<- 7.5*(7*c2 -1)*s2
    plg[3,6] <<- 3*c*plg[3,5]-2*plg[3,4]
    plg[3,7] <<- (11*c*plg[3,6]-7*plg[3,5])/4
    plg[3,8] <<- (13*c*plg[3,7]-8*plg[3,6])/5
    plg[4,4] <<- 15*s2*s
    plg[4,5] <<- 105*s2*s*c 
    plg[4,6] <<- (9*c*plg[4,5]-7.*plg[4,4])/2
    plg[4,7] <<- (11*c*plg[4,6]-8.*plg[4,5])/3
    if (!(((flags$sw[8] == 0) & (flags$sw[9] == 0)) & (flags$sw[15] == 0))) {
        stloc <- sin(hr*tloc)
        ctloc <- cos(hr*tloc)
        s2tloc <- sin(2*hr*tloc)
        c2tloc <- cos(2*hr*tloc)
        s3tloc <- sin(3*hr*tloc)
        c3tloc <- cos(3*hr*tloc)
    }
    cd32 <- cos(dr*(input$doy-p[32]))
    cd18 <- cos(2*dr*(input$doy-p[18]))
    cd14 <- cos(dr*(input$doy-p[14]))
    cd39 <- cos(2*dr*(input$doy-p[39]))
    # F10.7 effect
    df <- input$f107 - input$f107A
    dfa <- input$f107A - 150
    t[1] <-  p[20]*df*(1+p[60]*dfa) + p[21]*df*df + p[22]*dfa + p[30]*(dfa^2)
    f1 <- 1 + (p[48]*dfa +p[20]*df+p[21]*df*df)*flags$swc[2]
    f2 <- 1 + (p[50]*dfa+p[20]*df+p[21]*df*df)*flags$swc[2]
    # time independent
    t[2] = (p[2]*plg[1,3] + p[3]*plg[1,5] + p[23]*plg[1,7]) + 
        (p[15]*plg[1,3])*dfa*flags$swc[2] + p[27]*plg[1,2]
    # symmetrical annual
    t[3] <- p[19]*cd32
    # symmetrical semiannual
    t[4] <- (p[16] + p[17]*plg[1,3])*cd18
    # asymmetrical annual
    t[5] <- f1*(p[10]*plg[1,2] + p[11]*plg[1,4])*cd14
    # asymmetrical semiannual
    t[6] <- p[38]*plg[1,2]*cd39
    # diurnal
    if (flags$sw[8]) {
        t71 <- (p[12]*plg[2,3])*cd14*flags$swc[6]
        t72 <- (p[13]*plg[2,3])*cd14*flags$swc[6]
        t[7] <- f2*((p[4]*plg[2,2] + p[5]*plg[2,4] + p[28]*plg[2,6] + t71)*ctloc + 
                        (p[7]*plg[2,2] + p[8]*plg[2,4] + p[29]*plg[2,6] + t72)*stloc) 
    }
    # semidiurnal
    if (flags$sw[9]) {
        t81 <- (p[24]*plg[3,4]+p[36]*plg[3,6])*cd14*flags$swc[6]
        t82 <- (p[34]*plg[3,4]+p[37]*plg[3,6])*cd14*flags$swc[6]
        t[8] <- f2*((p[6]*plg[3,3]+ p[42]*plg[3,5] + t81)*c2tloc + 
                        (p[9]*plg[3,3] + p[43]*plg[3,5] + t82)*s2tloc)
    }
    # terdiurnal
    if (flags$sw[9]) {
        t[14] <- f2 * ((p[40]*plg[4,4]+(p[94]*plg[4,5]+p[47]*plg[4,7])*cd14*flags$swc[6])* s3tloc +
                           (p[41]*plg[4,4]+(p[95]*plg[4,5]+p[49]*plg[4,7])*cd14*flags$swc[6])*c3tloc)
    }
    # magnetic activity based on daily ap
    if (flags$sw[10] == -1) {
        ap <- input$ap_a
        if (p[52] != 0) {
            exp1 <- exp(-10800*sqrt(p[52]*p[52])/(1+p[139]*(45-abs(input$g_lat))))
            exp1 <- min(exp1, 0.99999)
            p[25] <- max(p[25], 1e-4)
        }
        apt[1] <<- sg0(exp1, p, ap$a)
        if(flags$sw[10]) {
            t[9] <- apt[1]*(p[51] + p[97]*plg[1,3] + p[55]*plg[1,5] +  
                                (p[126]*plg[1,2] + p[127]*plg[1,4] + p[128]*plg[1,6])*cd14*flags$swc[6] +  
                                (p[129]*plg[2,2] + p[130]*plg[2,4] + p[131]*plg[2,6])*flags$swc[8]* 
                                cos(hr*(tloc - p[132])))
        }
    } else {
        apd <- input$ap - 4
        p44 <- p[44]
        p45 <- p[45]
        if (p44<0)
            p44 <- 1E-5
        end
        apdf <- apd + (p45-1)*(apd + (exp(-p44 * apd) - 1)/p44)
        if (flags$sw[10])
            t[9] <- apdf*(p[33] + p[46]*plg(1,3) + p[35]*plg(1,5) +  ...
                       (p[101]*plg(1,2) + p[102]*plg(1,4) + p[103]*plg(1,6))*cd14*flags$swc[6] +  
                       (p[122]*plg(2,2) + p[123]*plg(2,4) + p[124]*plg(2,6))*flags$swc[8]* 
                       cos(hr*(tloc-p[125])))
    }
    if (((flags$sw[11]) & (input$g_long > -1000))) {
        # longitudinal
        if (flags$sw[12]) {
            t[11] <- (1 + p[81]*dfa*flags$swc[2])*((p[65]*plg[2,3]+p[66]*plg[2,5]+p[67]*plg[2,7] +
                        p[104]*plg[2,2]+p[105]*plg[2,4]+p[106]*plg[2,6] +
                        flags$swc[6]*(p[110]*plg[2,2]+p[111]*plg[2,4] + 
                                          p[112]*plg[2,6])*cd14)* cos(dgtr*flags$input$g_long)
                        + (p[91]*plg[2,3]+p[92]*plg[2,5]+p[93]*plg[2,7] + 
                               p[107]*plg[2,2]+p[108]*plg[2,4]+p[109]*plg[2,6] + 
                               flags$swc[6]*(p[113]*plg[2,2]+p[114]*plg[2,4] + 
                                                 p[115]*plg[2,6])*cd14)*sin(dgtr*flags$input$g_long))
        }
        # ut and mixed ut, longitude
        if (flags$sw[13]) {
            t[12] <- (1+p[96]*plg[1,2])*(1+p[82]*dfa*flags$swc[2])*
                (1+p[120]*plg[1,2]*flags$swc[6]*cd14)*
                ((p[69]*plg[1,2]+p[70]*plg[1,4]+p[71]*plg[1,6])*
                     cos(sr*(input$sec - p[72])))
            t[12] <- t[12] + flags$swc[12]* (p[77]*plg[3,4]+p[78]*
                plg[3,6]+p[79]*plg[3,8])* cos(sr*(input$sec-p[80])+2*dgtr*input$g_long) *
                (1+p[138]*dfa*flags$swc[2])
        }
        # ut, longitude magnetic activity
        if (flags$sw[14]) {
            if (flags$sw[10] == -1) {
                if(p[52]) {
                    t[13] <- apt(1)*flags$swc[12]*(1.+p[133]*plg[1,2])*
                        ((p[53]*plg[2,3]+p[99]*plg[2,5]+p[68]*plg[2,7])*
                        cos(dgtr*(input$g_long-p[98])))+ apt(1)*flags$swc[12]*
                        flags$swc[6]*(p[134]*plg[2,2]+p[135]*plg[2,4]+p[136]*plg[2,6])*
                        cd14*cos(dgtr*(input$g_long-p[137]))+apt(1)*flags$swc[13]*
                        (p[56]*plg[1,2]+p[57]*plg[1,4]+p[58]*plg[1,6])*cos(sr*(input$sec-p[59]))
                }
            } else {
                t[13] <- apdf*flags$swc[12]*(1+p[121]*plg[1,2])*
                    ((p[61]*plg[2,3]+p[62]*plg[2,5]+p[63]*plg[2,7])*
                         cos(dgtr*(input$g_long-p[64]))) + apdf*flags$swc[12]*flags$swc[6]* 
                    (p[116]*plg[2,2]+p[117]*plg[2,4]+p[118]*plg[2,6])* cd14*
                    cos(dgtr*(input$g_long-p[119])) + apdf*flags$swc[13]*
                    (p[84]*plg[1,2]+p[85]*plg[1,4]+p[86]*plg[1,6])*
                    cos(sr*(input$sec - p[76]))
            }
        }
    }
    tinf <- p[31]
    for (i in 1:14) {
        tinf <- tinf + abs(flags$sw[i+1])*t[i]
    }
    return(tinf)
}

glob7s <- function(p, input, flags) {
    # version of globe7 for lower atmosphere
    pset <- 2
    t <- numeric(14)
    dr <- 1.72142e-2
    dgtr <- 1.74533e-2
    if(p[100] == 0) p[100] <- pset
    if(p[100] != pset) {
        stop("Wrong parameter set for GLOB7S")
    }
    cd32 <- cos(dr*(input$doy-p[32]))
    cd18 <- cos(2*dr*(input$doy-p[18]))
    cd14 <- cos(dr*(input$doy-p[14]))
    cd39 <- cos(2*dr*(input$doy-p[39]))
    # effects F10.7
    t[1] <- p[22] * dfa
    # time independent
    t[2] <- p[2]*plg[1,3] + p[3]*plg[1,5] + p[23]*plg[1,7] + p[27]*plg[1,2] + p[15]*plg[1,4] + p[60]*plg[1,6]
    # symmetrical annual
    t[3] <- (p[19]+p[48]*plg[1,3]+p[30]*plg[1,5])*cd32
    # symmetrical semiannual
    t[4] <- (p[16]+p[17]*plg[1,3]+p[31]*plg[1,5])*cd18
    # asymmetrical annual
    t[5] <- (p[10]*plg[1,2]+p[11]*plg[1,4]+p[21]*plg[1,6])*cd14
    # asymmetrical semiannual
    t[6] <- (p[38]*plg[1,2])*cd39
    # diurnal
    if (flags$sw[8]) {
        t71 <- p[12]*plg[2,3]*cd14*flags$swc[6]
        t72 <- p[13]*plg[2,3]*cd14*flags$swc[6]
        t[7] <- ((p[4]*plg[2,2] + p[5]*plg[2,4] + t71) * ctloc + (p[7]*plg[2,2] + p[8]*plg[2,4] + t72) * stloc)
    }
    # semidiurnal
    if (flags$sw[9]) {
        t81 <- (p[24]*plg[3,4]+p[36]*plg[3,6])*cd14*flags$swc[6]
        t82 <- (p[34]*plg[3,4]+p[37]*plg[3,6])*cd14*flags$swc[6]
        t[8] <- ((p[6]*plg[3,3] + p[42]*plg[3,5] + t81) * c2tloc + (p[9]*plg[3,3] + p[43]*plg[3,5] + t82) * s2tloc)
    }
    # terdiurnal
    if (flags$sw[15]) {
        t[14] <- p[40] * plg[4,4] * s3tloc + p[41] * plg[4,4] * c3tloc
    }
    # magnetic activity
    if (flags$sw[10]) {
        if (flags.sw[10] == 1) {
            t[9] <- apdf * (p[33] + p[46] * plg[1, 3] * flags$swc[3])
        } else if (flags.sw[10] == -1) {
            t[9] <- (p[51]*apt[1] + p[97]*plg[1,3] * apt[1]*flags$swc[3])
        }
    }
    # longitudinal
    if (!((flags$sw[11] == 0) | (flags$sw[12] == 0) | (input$g_long <= -1000))) {
        t[11] <- (1 + plg[1,2]*(p[81]*flags$swc[6]*cos(dr*(input$doy-p[82]))
                                + p[86]*flags$swc[7]*cos(2*dr*(input$doy-p[87])))
                  + p[84]*flags$swc[4]*cos(dr*(input$doy-p[85]))
                  + p[88]*flags$swc[5]*cos(2*dr*(input$doy-p[89]))) *
         ((p[65]*plg[2,3]+p[66]*plg[2,5]+p[67]*plg[2,7] +
                p[75]*plg[2,2]+p[76]*plg[2,4]+p[77]*plg[2,6]) *
               cos(dgtr*input$g_long) + (p[91]*plg[2,3]+p[92]*plg[2,5]+p[93]*plg[2,7] +
                                             p[78]*plg[2,2]+p[79]*plg[2,4]+p[80]*plg[2,6])
           * sin(dgtr*input$g_long))
    }
    tt <- 0
    for (i in 1:14) {
        tt <- tt + abs(flags$sw[i+1])*t[i]
    }
    return(tt)
}

gts7 <- function(input, flags) {
    # Thermospheric portion of NRLMSISE-00 (altitude higher than 72.5 km)
    zn1 <- c(120.0, 110.0, 100.0, 90.0, 72.5)
    mn1 <- 5
    dgtr <- 1.74533e-2
    dr <- 1.72142e-2
    alpha <- c(-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0)
    altl <- c(200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0)
    za <- pdl[2, 16]
    zn1[1] <- za
    output <- list(
        d = numeric(9),
        t = numeric(2)
    )
    if (input$alt > zn1[1]) { # tinf variations not important below za or zn1[1]
        tinf <- ptm[1] * pt[1] * (1 + flags$sw[17] * globe7(pt,input,flags))
    } else {
        tinf <- ptm[1] * pt[1]
    }
    output$t[1] <- tinf
    if (input$alt > zn1[5]) { # gradient variations not important below zn1[5]
        g0 <- ptm[4] * ps(1) * (1 + flags$sw[20] * globe7(ps,input,flags))
    } else {
        g0 <- ptm[4] * ps(1)
    }
    tlb <- ptm[2] * (1 + flags$sw[18] * globe7(pd[4, ], input, flags)) * pd[4, 1]
    s <- g0 / (tinf - tlb)
    # Lower thermosphere temp variations not significant for density above 300 km
    if (input$alt < 300) {
        meso_tn1[2] <- ptm[7]*ptl[1,1]/(1-flags$sw[19]*glob7s(ptl[1,], input, flags))
        meso_tn1[3] <- ptm[3]*ptl[2,1]/(1-flags$sw[19]*glob7s(ptl[2,], input, flags))
        meso_tn1[4] <- ptm[8]*ptl[3,1]/(1-flags$sw[19]*glob7s(ptl[3,], input, flags))
        meso_tn1[5] <- ptm[5]*ptl[4,1]/(1-flags$sw[19]*flags$sw[21]*glob7s(ptl[4, ], input, flags))
        meso_tgn1[2] <- ptm[9]*pma[9,1]*(1+flags$sw[19]*flags$sw[21]*glob7s(pma[9, ], input, flags))*
            meso_tn1[5]*meso_tn1[5]/((ptm[5]*ptl[4,1])^2)
    } else {
        meso_tn1[2] <- ptm[7]*ptl[1,1]
        meso_tn1[3] <- ptm[3]*ptl[2,1]
        meso_tn1[4] <- ptm[8]*ptl[3,1]
        meso_tn1[5] <- ptm[5]*ptl[4,1]
        meso_tgn1[2] <- ptm[9] * pma[9,1] * meso_tn1[5]*meso_tn1[5]/((ptm[5]*ptl[4,1])^2)
    }
    # N2 variation factor at Zlb
    g28 <- flags$sw[22] * globe7(pd[3, ], input, flags)
    # variation of turbopause height
    zhf <- pdl[2, 25] * (1 + flags$sw[6] * pdl[1, 25] * sin(dgtr * input$g_lat) * 
                             cos(dr * (input$doy-pt[14])))
    ### output$t[1] <- tinf
    xmm <- pdm[3, 5]
    z <- input$alt
    # N2 density
    db28 <- pdm[3,1] * exp(g28) * pd[3, 1]
    densu_output <- densu(z,db28,tinf,tlb,28,alpha(3),
                          output$t[2],ptm(6),s,mn1,zn1,
                          meso_tn1,meso_tgn1)
    output$d[3] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    dd <<- output$d[3]
    zh28 <- pdm[3, 3] * zhf
    zhm28 <- pdm[3, 4] * pdl[2, 6] 
    xmd <- 28 - xmm
    densu_output <- densu(zh28,db28,tinf,tlb,xmd,(alpha(3)-1),
                          tz,ptm[6],s,mn1, zn1,
                          meso_tn1,meso_tgn1)
    b28 <- densu_output$densu_temp
    tz <- densu_output$tz
    if ((flags$sw[16]) & (z <= altl[3])) {
        dm28 <<- densu(z,b28,tinf,tlb,xmm,alpha[3],
                       tz,ptm[6],s,mn1,zn1,
                       meso_tn1,meso_tgn1)$densu_temp
        output$d[3] <- dnet(output$d[3], dm28, zhm28, xmm, 28)
    }
    # He density
    g4 <- flags$sw[22] * globe7(pd[1, ], input, flags)
    db04 <- pdm[1, 1] * exp(g4) * pd[1,1]
    densu_output <- densu(z,db04,tinf,tlb, 4,alpha[1],
                          output$t[2],ptm[6],s,mn1,zn1,
                          meso_tn1,meso_tgn1)
    output$d[1] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    dd <<- output$d[1]
    if ((flags$sw[16]) & ( z<altl[1]) ) {
        zh04 <- pdm[1,3]
        densu_output <- densu(zh04,db04,tinf,tlb,4-xmm,alpha[1]-1.0,
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        b04 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b04,tinf,tlb,xmm,0,
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        dm04 <<- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        zhm04 <- zhm28
        output$d[1] <- dnet(output$d[1],dm04,zhm04,xmm,4)
        rl <- log(b28*pdm[1,2]/b04)
        zc04 <- pdm[1,5]*pdl[2,1]
        hc04 <- pdm[1,6]*pdl[2,2]
        output$d[1] <- output$d[1]*ccor(z,rl,hc04,zc04)
    }
    # O density
    g16 <- flags.sw[22] * globe7(pd[2, ], input, flags)
    db16 <- pdm[2,1] * exp(g16) * pd[2, 1]
    densu_output <- ensu(z,db16,tinf,tlb, 16,alpha[2],
                         output$t[2],ptm[6],s,mn1,zn1,
                         meso_tn1,meso_tgn1)
    output$d[2] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    dd <<- output$d[2]
    if ((flags$sw[16]) & ( z<=altl[2] )) {
        zh16 <- pdm[2, 3]
        densu_output <- densu(zh16,db16,tinf,tlb,16-xmm,(alpha[2]-1),
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        b16 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b16,tinf,tlb,xmm,0,
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        dm16 <<- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        zhm16 <- zhm28
        output$d[2] <- dnet(output$d[2],dm16,zhm16,xmm,16)
        rl <- pdm[2,2]*pdl[2,17]*(1+flags.sw[2]*pdl[1,24]*(input$f107A-150))
        hc16 <- pdm[2,6]*pdl[2,4]
        zc16 <- pdm[2,5]*pdl[2,3]
        hc216 <- pdm[2,6]*pdl[2,5]
        output$d[2] <- output$d[2]*ccor2(z,rl,hc16,zc16,hc216)
        hcc16 <- pdm[2,8]*pdl[2,14]
        zcc16 <- pdm[2,7]*pdl[2,13]
        rc16 <- pdm[2,4]*pdl[2,15]
        output$d[2] <- output$d[2]*ccor(z,rc16,hcc16,zcc16)
    }
    # O2 density
    g32 <- flags$sw[22]*globe7(pd[5, ], input, flags)
    db32 <- pdm[4,1] * exp(g32) * pd[5,1]
    densu_output <- densu(z,db32,tinf,tlb,32,alpha[4],
                          output$t[2],ptm[6],s,mn1,zn1,
                          meso_tn1,meso_tgn1)
    dd <<- output$d[4]
    if (flags$sw[16]) {
        if (z <= altl[4]) {
            zh32 <- pdm[4,3]
            densu_output <- densu(zh32,db32,tinf,tlb,32-xmm,alpha[4]-1,
                                  output$t[2],ptm[6],s,mn1,zn1,
                                  meso_tn1,meso_tgn1)
            b32 <- densu_output$densu_temp
            output$t[2] <- densu_output$tz
            densu_output <- densu(z,b32,tinf,tlb,xmm,0,output$t[2],
                                  ptm[6],s,mn1,zn1,
                                  meso_tn1,meso_tgn1)
            dm32 <<- densu_output$densu_temp
            output$t[2] <- densu_output$tz
            zhm32 <- zhm28
            output$d[4] <- dnet(output$d[4],dm32,zhm32,xmm,32)
            rl <- log(b28*pdm[4,2]/b32)
            hc32 <- pdm[4,6]*pdl[2,8]
            zc32 <- pdm[4,5]*pdl[2,7]
            output$d[4] <- output$d[4]*ccor(z,rl,hc32,zc32)
        }
        hcc32 <- pdm[4,8]*pdl[2,23]
        hcc232 <- pdm[4,8]*pdl[1,23]
        zcc32 <- pdm[4,7]*pdl[2,22]
        rc32 <- pdm[4,4]*pdl[2,24]*(1+flags$sw(2)*pdl[1,24]*(input$f107A - 150))
        output$d[4] <- output$d[4]*ccor2(z,rc32,hcc32,zcc32,hcc232)
    }
    # Ar density
    g40 <- flags$sw[22] * globe7(pd[6, ],input,flags)
    db40 <- pdm[5,1]*exp(g40)*pd[6, 1]
    densu_output <- densu(z,db40,tinf,tlb,40,alpha[5],
                          output$t[2],ptm[6],s,mn1,zn1,
                          meso_tn1,meso_tgn1)
    output$d[5] <- densu_output$densu_temp
    outpu$t[2] <- densu_output$tz
    dd <<- output$d[5]
    if ((flags$sw[16]) & (z <= altl[5])) {
        zh40 <- pdm[5, 3]
        densu_output <- densu(zh40,db40,tinf,tlb,40-xmm,alpha[5]-1,
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        b40 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b40,tinf,tlb,xmm,0,
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        dm40 <<- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        zhm40 <- zhm28
        output$d[5] <- dnet(output$d[5], dm40, zhm40, xmm, 40)
        rl <- log(b28*pdm[5,2]/b40)
        hc40 <- pdm[5, 6] * pdl[2, 10]
        zc40 <- pdm[5, 5] * pdl[2, 9]
        output$d[5] <- output$d[5] * ccor(z,rl,hc40,zc40)
    }
    # H density
    g1 <- flags$sw[22] * globe7(pd[7, ], input, flags)
    db01 <- pdm[6, 1] * exp(g1) * pd[7, 1]
    densu_output <- densu(z,db01,tinf,tlb,1,alpha[7],
                          output$t[2],ptm[6],s,mn1,zn1,meso_tn1,meso_tgn1)
    output$d[7] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    dd <<- output$d[7]
    if ((flags$sw[16]) & (z <= altl[7])) {
        zh01 <- pdm[6, 3]
        densu_output <- densu(zh01,db01,tinf,tlb,1-xmm,alpha[7]-1,
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        b01 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b01,tinf,tlb,xmm,0,output$t[2],
                              ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        dm01 <<- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        zhm01 <- zhm28
        output$d[7] <- dnet(output$d[7],dm01,zhm01,xmm,1)
        rl <- log(b28*pdm[6,2]*sqrt(pdl[2,18]*pdl[2,18])/b01)
        hc01 <- pdm[6,6]*pdl[2,12]
        zc01 <- pdm[6,5]*pdl[2,11]
        output$d[7] <- output$d[7]*ccor(z,rl,hc01,zc01)
        hcc01 <- pdm[6,8]*pdl[1,20]
        zcc01 <- pdm[6,7]*pdl[2,19]
        rc01 <- pdm[6,4]*pdl[2,21]
        output$d[7] <- output$d[7]*ccor(z,rc01,hcc01,zcc01)
    }
    # Atomic nitrogen (N) density
    g14 <- flags$sw[22]*globe7(pd[8, ],input,flags)
    db14 <- pdm[7,1]*exp(g14)*pd[8,1]
    densu_output <- densu(z,db14,tinf,tlb,14,alpha[8],
                          output$t[2],ptm[6],s,mn1,zn1,
                          meso_tn1,meso_tgn1)
    output$d[8] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    dd <<- output$d[8]
    if ((flags$sw[16]) & (z <= altl[8])) {
        zh14 <- pdm[7, 3]
        densu_output <- densu(zh14,db14,tinf,tlb,14-xmm,alpha[8]-1,
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        b14 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b14,tinf,tlb,xmm,0,
                              output$t[2],ptm[6],s,mn1,zn1,
                              meso_tn1,meso_tgn1)
        dm14 <<- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        zhm14 <- zhm28
        output$d[8] <- dnet(output$d[8],dm14,zhm14,xmm,14)
        rl <- log(b28*pdm[7,2]*sqrt(pdl[1,3]*pdl[1,3])/b14)
        hc14 <- pdm[7,6]*pdl[1,2]
        zc14 <- pdm[7,5]*pdl[1,1]
        output$d[8] <- output$d[8]*ccor(z,rl,hc14,zc14)
        hcc14 <- pdm[7,8]*pdl[1,5]
        zcc14 <- pdm[7,7]*pdl[1,4]
        rc14 <- pdm[7,4]*pdl[1,6]
        output$d[8] <- output$d[8]*ccor(z,rc14,hcc14,zcc14)
    }
    # Anomalous oxygen density
    g16h <- flags$sw[22]*globe7(pd[9, ],input,flags)
    db16h <- pdm[8,1]*exp(g16h)*pd[9,1]
    tho <- pdm[8,10]*pdl[1,7]
    densu_output <- densu(z,db16h,tho,tho,16,alpha[9],
                          output$t[2],ptm[6],s,mn1,zn1,
                          meso_tn1,meso_tgn1)
    dd <<- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    zsht <- pdm[8,6]
    zmho <- pdm[8,5]
    zsho <- scalh(zmho,16,tho,gsurf)
    output$d[9] <- dd*exp(-zsht/zsho*(exp(-(z-zmho)/zsht)-1))
    # Total density (atomic mass-weighted sums)
    output$d[6] <- 1.66e-24*(4*output$d[1]+16*output$d[2]+28*output$d[3]+32*output$d[4]+40*output$d[5]+ output$d[7]+14*output$d[8])
    # Temperature
    z <- abs(input$alt)
    densu_output <- densu(z,1,tinf,tlb,0,0,
                          output$t[2],ptm[6],s,mn1,zn1, 
                          meso_tn1,meso_tgn1)
    output$t[2] <- densu_output$tz
    if (flags$sw[1]) {
        for(i in 1:9) {
            output$d[i] <- output$d[i]*1e6;
        }
        output$d[6] <- output$d[6]/1000
    }
    return(output)
}

gtd7 <- function(input, flags) {
    for (i in 1:24) {
        if (i == 10) {
            flags$sw[i] <<- flags$switches[i]
            flags$swc[i] <<- flags$switches[i]
        } else {
            flags$sw[i] <<- if(flags$switches[i] == 1) 1 else 0
            flags$swc[i] <<- if(flags$switches[i] > 0) 1 else 0
        }
    }
    mn3 <- 5
    zn3 <- c(32.5, 20.0, 15.0, 10.0, 0.0)
    mn2 <- 4
    zn2 <- c(72.5, 55.0, 45.0, 32.5)
    zmix <- 62.5
    output <- list(
        d = numeric(9),
        t = numeric(2)
    )
    tz <- numeric(1)
    xlat <- if(flags$sw[3] == 0) 45 else input$g_lat
    glatf_results <- glatf(xlat)
    gsurf <<- glatf_results[[1]]
    re <<- glatf_results[[2]]
    xmm <- pdm[3,5]
    altt <- max(input$alt, zn2[1])
    tmp <- input$alt
    input$alt <- altt
    soutput <- gts7(input, flags)
    input$alt <- tmp
    if (flags$sw[1]) {
        dm28m <- dm28 * 1e6
    } else {
        dm28m <- dm28
    }
    output$t <- soutput$t
    if (input$alt >= zn2[1]) {
        output$d <- soutput$d
    } else { # lower mesosphere/upper stratosphere (between zn3[1] and zn2[1])
        meso_tgn2[1] <- meso_tgn1[2];
        meso_tn2[1] <- meso_tn1[5];
        meso_tn2[2] <- pma[1,1]*pavgm[1]/(1-flags$sw[21]*glob7s(pma[1, ], input, flags))
        meso_tn2[3] <- pma[2,1]*pavgm[2]/(1-flags$sw[21]*glob7s(pma[2, ], input, flags))
        meso_tn2[4] <- pma[3,1]*pavgm[3]/(1-flags$sw[21]*flags$sw[23]*glob7s(pma[3, ], input, flags))
        meso_tgn2[2] <- pavgm[9]*pma[10,1]*(1+flags$sw[21]*flags$sw[23]*
                                            glob7s(pma[10, ], input, flags))*
                        meso_tn2[4]*meso_tn2[4]/((pma[3,1]*pavgm[3])^2)
        meso_tn3[1] <- meso_tn2[4]
        if (input.alt<zn3(1)) { # lower stratosphere/troposphere (below zn3[1])
            meso_tgn3[1] <- meso_tgn2[2]
            meso_tn3[2] <- pma[4,1]*pavgm[4]/(1-flags $sw[23]*glob7s(pma[4, ], input, flags))
            meso_tn3[3] <- pma[5,1]*pavgm[5]/(1-flags $sw[23]*glob7s(pma[5, ], input, flags))
            meso_tn3[4] <- pma[6,1]*pavgm[6]/(1-flags $sw[23]*glob7s(pma[6, ], input, flags))
            meso_tn3[5] <- pma[7,1]*pavgm[7]/(1-flags $sw[23]*glob7s(pma[7, ], input, flags))
            meso_tgn3[2] <- pma[8,1]*pavgm[8]*(1+flags $sw[23]*
                                               glob7s(pma[8, ], input, flags)) *
                            meso_tn3[5]*meso_tn3[5]/((pma[7,1]*pavgm[7])^2)
        } 
        # linear transition to full mixing below zn2[1]
        dmc <- 0
        if (input$alt > zmix) {
            dmc <- 1 - (zn2[1] - input$alt)/(zn2[1] - zmix)
        }
        dz28 <- soutput$d[3]
        ## N2 density 
        dmr <- soutput$d[3] / dm28m - 1
        densm_results <- densm(input$alt, dm28m, xmm, tz, mn3, zn3, 
                               meso_tn3, meso_tgn3, mn2, zn2, 
                               meso_tn2, meso_tgn2)
        output$d[3] <- densm_results$densm_tmp
        tz < densm_results$tz
        output$d[3] <- output$d[3] * (1 + dmr*dmc)
        ## He density
        dmr <- soutput$d[1] / (dz28 * pdm[1,2]) - 1
        output$d[1] <- output$d[3] * pdm[1,2] * (1 + dmr*dmc)
        ## O density
        output$d[c(2, 9)] <- 0
        ## O2 density
        dmr <- soutput$d[4] / (dz28 * pdm[4,2]) - 1
        output$d[4] <- output$d[3] * pdm[4,2] * (1 + dmr*dmc)
        # Ar density
        dmr <- soutput$d[5] / (dz28 * pdm[5,2]) - 1
        output$d[5] <- output$d[3] * pdm[5,2] * (1 + dmr*dmc)
        # H density
        output$d[7] <- 0
        # Atomic nitrogen (N) density
        output$d[8] <- 0
        # Total mass density
        output$d[6] <- 1.66e-24 * (4 * output$d[1] + 16 * output$d[2] + 
                                   28 * output$d[3] + 32 * output$d[4] + 40 * 
                                   output$d[5] + output$d[7] + 14 * output$d[8])
        if (flags$sw[1]) {
            output$d[6] <- output$d[6]/1000
        }
        densm_results <- densm(input$alt, 1, 0, tz, mn3, zn3, 
                               meso_tn3, meso_tgn3, mn2, zn2, 
                               meso_tn2, meso_tgn2);
        dd <<- densm_results$densm_temp
        tz <- densm_results$tz
        output$t[2] <- tz
    }
    return(output)
}

gtd7d <- function(input, flags) {
    output <- gtd7(input, flags)
    output$d[6] <- 1.66e-24 * (4 * output$d[1] + 16 * output$d[2] + 
                               28 * output$d[3] + 32 * output$d[4] + 
                               40 * output$d[5] + output$d[7] + 
                               14 * output$d[8] + 16 * output$d[9])
    if (flags$sw[1]) {
        output$d[6] <- output$d[6]/1000
    }
    return(output)
}

NRLMSISE00 <- function(Mjd_UTC, r_ECEF, UT1_UTC, TT_UTC) {
    gsurf <<- 0
    # re <- 6367.088132098377173
    dd <<- 0
    dm04 <<- 0
    dm16 <<- 0
    dm28 <<- 0
    dm32 <<- 0
    dm40 <<- 0
    dm01 <<- 0
    dm14 <<- 0
    meso_tn1 <<- rep(0, 5)
    meso_tn2 <<- rep(0, 4)
    meso_tn3 <<- rep(0, 5)
    meso_tgn1 <<- rep(0, 2)
    meso_tgn2 <<- rep(0, 2)
    meso_tgn3 <<- rep(0, 2)
    dfa <<- 0
    plg <<- matrix(0, nrow=4, ncol=9)
    ctloc <<- 0
    stloc <<- 0
    c2tloc <<- 0
    s2tloc <<- 0
    s3tloc <<- 0
    c3tloc <<- 0
    apdf <<- 0
    apt <<- rep(0, 4)
    flags <- list(
        switches = c(0, rep(1, 8), -1, rep(1, 14)),
        sw = rep(0, 24),
        swc = rep(0, 24)
    )
    ap_a <- numeric(7)
    input <- list(
        year = 0,
        doy = 0,
        sec = 0,
        alt = 0,
        g_lat = 0,
        g_long = 0,
        lst = 0,
        f107A = 0,
        f107 = 0,
        ap = 0,
        ap_a = ap_a
    )
    invjday_results <- invjday(Mjd_UTC+2400000.5)
    days <- fractionalDays(invjday_results$year,
                           invjday_results$month,
                           invjday_results$day,
                           invjday_results$hour,
                           invjday_results$min,
                           invjday_results$sec)
    input$doy <- floor(days)
    input$year <- 0
    input$sec <- invjday_results$hour*3600 + invjday_results$min*60 + invjday_results$sec
    geodetic_results <- ECEFtoGeodetic(r_ECEF)
    input$alt <- geodetic_results$height/1000
    input$g_lat <- geodetic_results$lat*(180/pi)
    input$g_long <- geodetic_results$lon*(180/pi)
    iauCal2jd_results <- iauCal2jd(invjday_results$year, invjday_results$month, invjday_results$day)

    TIME <- (60*(60*invjday_results$hour + invjday_results$minute) + invjday_results$sec)/86400
    UTC <- iauCal2jd_results$DATE + TIME
    TT <- UTC + TT_UTC/86400
    TUT <- TIME + UT1_UTC/86400
    UT1 <- iauCal2jd_results$DATE + TUT
    rnpb <- iauPnm06a(iauCal2jd_results$DJMJD0, TT)
    lst <- geodetic_results$lon + iauGst06(iauCal2jd_results$DJMJD0, UT1, 
                                           iauCal2jd_results$DJMJD0, TT, rnpb)
    lst <- lst %% 2*pi
    lst <- (lst*24)/(2*pi)
    input$lst <- lst
    i <- which(((invjday_results$year == spaceWeather[, 1]) & 
                    (invjday_results$mon == spaceWeather[, 2]) & 
                    (invjday_results$day == spaceWeather[, 3])))[1]
    ## BORRAR ESTE COMENTARIO: variable spaceWeather equivalente a global swdata
    ## del HPOP de matlab, que se define en archivo test_HPOP.m (en versiones an
    ## teriores, conocido simplemente como HPOP.m)
    sw <- spaceWeather[i, ]
    input$ap <- sw[23]
    input$ap_a[1] <- sw[23] 
    input$ap_a[2] <- sw[15] 
    sw_1 <- spaceWeather[i-1, ]
    input$ap_a[3] <- sw_1[22] 
    input$ap_a[4] <- sw_1[21] 
    input$ap_a[5] <- sw_1[20] 
    sum <- sw_1[19]+sw_1[18]+sw_1[17]+sw_1[16]+sw_1[15]
    sw_2 <- spaceWeather[i-2, ]
    sum <- sum+sw_2[22]+sw_2[21]+sw_2[20]
    input$ap_a[6] <- sum/8 
    sw_3 <- spaceWeather[i-3, ]
    sum <- sw_2[19] + sw_2[18] + sw_2[17] + sw_2[16] + 
          sw_2[15] + sw_3[22] + sw_3[21] + sw_3[20]
    input$ap_a[7] <- sum/8
    input$f107 <- sw[31]  
    input.f107A <- sw[33]   
    output <- gtd7d(input, flags)
    dens <- 1e3*output$d(6)
    return(dens)
}
