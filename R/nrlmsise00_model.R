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
    rgas <- asteRiskData::rgas
    g <- rgas * temp / (gsurf / ((1 + alt/NRLMSISE00.env$re_)^2) * xm)
    return(g)
}

dnet <- function(dd, dm, zhm, xmm, xm) {
    # turbopause correction for msis models
    # Root mean density
    # NRLMSISE00.env$dd - diffusive density
    # dm - full mixed density
    # zhm - transition scale length
    # xmm - full mixed molecular weight
    # xm  - species molecular weight
    # Outputs DNET - combined density
    a <- zhm / (xmm-xm)
    if ((dm <= 0) | (NRLMSISE00.env$dd <= 0)) {
        warning(paste("dnet log error", dm, NRLMSISE00.env$dd, xm, sep=" "))
        if (dm == 0) {
            dnet <- NRLMSISE00.env$dd
            if (NRLMSISE00.env$dd == 0){
                assign("dd", 1, NRLMSISE00.env)
            }
        } else if (NRLMSISE00.env$dd == 0) {
            dnet <- dm
        }
    }
    ylog <- a * log(dm/NRLMSISE00.env$dd)
    if (ylog < -10) {
        dnet <- NRLMSISE00.env$dd
    } else if (ylog > 10) {
        dnet <- dm
    } else {
        dnet <- NRLMSISE00.env$dd * (1 + exp(ylog))^(1/a)
    }
    return(dnet)
}

zeta <- function(zz, zl) {
    result <- ((zz - zl) * (NRLMSISE00.env$re_ + zl))/ (NRLMSISE00.env$re_ + zz)
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
        pos <- xs != 0
        xs <- xs[pos]
        ys <- ys[pos]
        # for (k in 1:mn) {
        #     xs[k] <- zeta(zn2[k], z1)/zgdif
        #     ys[k] <- 1/tn2[k]
        # }
        #yd1 <- -tgn2[1] / (t1*t1) * zgdif
        #yd2 <- -tgn2[2] / (t2*t2) * zgdif *  (((NRLMSISE00.env$re_+z2)/(NRLMSISE00.env$re_+z1))^2)
        # Calculate spline coefficients
        #y2out <- spline(xs, ys, mn, yd1, yd2)
        cubicspline <- splinefun(xs, ys)
        x <- zg/zgdif
        #y <- splint(xs, ys, y2out, mn, x)
        y <- cubicspline(x)
        # Temperature at altitude
        tz <- 1/y
        if (xm != 0) {
            glb <- NRLMSISE00.env$gsurf/((1 + z1/NRLMSISE00.env$re_)^2)
            gamm <- xm*glb*zgdif/asteRiskData::rgas
            # yi <- splini(xs, ys, y2out, mn, x)
            yi <- integrate(cubicspline, xs[1], x)
            expl <- gamm*yi$value
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
            pos <- xs != 0
            xs <- xs[pos]
            ys <- ys[pos]
            cubicspline <- splinefun(xs, ys)
            x <- zg/zgdif
            y <- cubicspline(x)
            tz <- 1/y
            if (xm != 0) {
                glb <- NRLMSISE00.env$gsurf/((1 + z1/NRLMSISE00.env$re_)^2)
                gamm <- xm*glb*zgdif/asteRiskData::rgas
                yi <- integrate(cubicspline, xs[1], x)
                expl <- gamm*yi$value
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
        dta <- (tinf - ta) * s2 * ((NRLMSISE00.env$re_ + zlb)/(NRLMSISE00.env$re_ + za))^2
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
        pos <- xs != 0
        xs <- xs[pos]
        ys <- ys[pos]
        cubicspline <- splinefun(xs, ys)
        # yd1 <- -tgn1[1] / (t1*t1) * zgdif
        # yd2 <- -tgn1[2] / (t2*t2) * zgdif * ((NRLMSISE00.env$re_+z2)/(NRLMSISE00.env$re_+z1))^2
        # y2out <- spline(xs, ys, mn, yd1, yd2)
        x <- zg/zgdif
        # y <- splint(xs, ys, y2out, mn, x)
        y <- cubicspline(x)
        tz <- 1/y
        densu_temp <- tz
    }  
    if (xm != 0) {
        glb <- NRLMSISE00.env$gsurf/((1+zlb/NRLMSISE00.env$re_)^2)
        gamma <- xm * glb / (s2 * asteRiskData::rgas * tinf)
        expl <- exp(-s2 * gamma * zg2)
        if(expl > 50 | tt <= 0) expl <- 50
        densa <- dlb * (tlb/tt)^(1+alpha+gamma) * expl
        densu_temp <- densa
        if (alt < za) {
            glb <- NRLMSISE00.env$gsurf/((1+z1/NRLMSISE00.env$re_)^2)
            gamm <- xm * glb * zgdif / asteRiskData::rgas
            # yi <- splini(xs, ys, y2out, mn, x)
            yi <- integrate(cubicspline, xs[1], x)
            expl <- gamm * yi$value
            if(expl > 50 | tt <= 0) expl <- 50
            densu_temp <- densu_temp * (t1 / tz)^(1 + alpha) * exp(-expl)
        }
    }
    return(list(
        densu_temp = densu_temp,
        tz = tz))
}

g0 <- function(a, p) {
    a <- a[[1]]
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
    new_plg_values <- c(0, c, 0.5*(3*c2 -1), 0.5*(5*c*c2-3*c), (35*c4 - 30*c2 + 3)/8,
                        (63*c2*c2*c - 70*c2*c + 15*c)/8, (11*c*((63*c2*c2*c - 70*c2*c + 15*c)/8) - 
                                                          5*((35*c4 - 30*c2 + 3)/8))/6, 
                        0, 0, # end of row 1
                        0, s, 3*c*s, 1.5*(5*c2-1)*s, 2.5*(7*c2*c-3*c)*s,
                        1.875*(21*c4 - 14*c2 +1)*s, (11*c*(1.875*(21*c4 - 14*c2 +1)*s)-
                                                         6*(2.5*(7*c2*c-3*c)*s))/5, 
                        0, 0, # end of row 2
                        0, 0, 3*s2, 15*s2*c, 7.5*(7*c2 -1)*s2, 
                        3*c*(7.5*(7*c2 -1)*s2)-2*(15*s2*c),
                        (11*c*(3*c*(7.5*(7*c2 -1)*s2)-2*(15*s2*c))-7*(7.5*(7*c2 -1)*s2))/4,
                        (13*c*(11*c*(3*c*(7.5*(7*c2 -1)*s2)-2*(15*s2*c))-7*(7.5*(7*c2 -1)*s2))/4 -
                             8*(3*c*(7.5*(7*c2 -1)*s2)-2*(15*s2*c)))/5, 0, # end of row 3
                        0, 0, 0, 15*s2*s, 105*s2*s*c, 
                        (9*c*(105*s2*s*c)-7*(15*s2*s))/2, 
                        (11*c*((9*c*(105*s2*s*c)-7*(15*s2*s))/2)-8*(105*s2*s*c))/3, 0, 0) # end of row 4
    assign("plg", matrix(new_plg_values, nrow=4, ncol=9, byrow=TRUE), NRLMSISE00.env)
    if (!(((flags$sw[8] == 0) & (flags$sw[9] == 0)) & (flags$sw[15] == 0))) {
        assign("stloc", sin(hr*tloc), NRLMSISE00.env)
        assign("ctloc", cos(hr*tloc), NRLMSISE00.env)
        assign("s2tloc", sin(2*hr*tloc), NRLMSISE00.env)
        assign("c2tloc", cos(2*hr*tloc), NRLMSISE00.env)
        assign("s3tloc", sin(3*hr*tloc), NRLMSISE00.env)
        assign("c3tloc", cos(3*hr*tloc), NRLMSISE00.env)
    }
    cd32 <- cos(dr*(input$doy-p[32]))
    cd18 <- cos(2*dr*(input$doy-p[18]))
    cd14 <- cos(dr*(input$doy-p[14]))
    cd39 <- cos(2*dr*(input$doy-p[39]))
    # F10.7 effect
    df <- as.numeric(input$f107 - input$f107A)
    assign("dfa", as.numeric(input$f107A - 150), NRLMSISE00.env)
    t[1] <-  p[20]*df*(1+p[60]*NRLMSISE00.env$dfa) + p[21]*df*df + p[22]*NRLMSISE00.env$dfa + p[30]*(NRLMSISE00.env$dfa^2)
    f1 <- 1 + (p[48]*NRLMSISE00.env$dfa+p[20]*df+p[21]*df*df)*flags$swc[2]
    f2 <- 1 + (p[50]*NRLMSISE00.env$dfa+p[20]*df+p[21]*df*df)*flags$swc[2]
    # time independent
    t[2] = (p[2]*NRLMSISE00.env$plg[1,3] + p[3]*NRLMSISE00.env$plg[1,5] + p[23]*NRLMSISE00.env$plg[1,7]) + 
        (p[15]*NRLMSISE00.env$plg[1,3])*NRLMSISE00.env$dfa*flags$swc[2] + p[27]*NRLMSISE00.env$plg[1,2]
    # symmetrical annual
    t[3] <- p[19]*cd32
    # symmetrical semiannual
    t[4] <- (p[16] + p[17]*NRLMSISE00.env$plg[1,3])*cd18
    # asymmetrical annual
    t[5] <- f1*(p[10]*NRLMSISE00.env$plg[1,2] + p[11]*NRLMSISE00.env$plg[1,4])*cd14
    # asymmetrical semiannual
    t[6] <- p[38]*NRLMSISE00.env$plg[1,2]*cd39
    # diurnal
    if (flags$sw[8]) {
        t71 <- (p[12]*NRLMSISE00.env$plg[2,3])*cd14*flags$swc[6]
        t72 <- (p[13]*NRLMSISE00.env$plg[2,3])*cd14*flags$swc[6]
        t[7] <- f2*((p[4]*NRLMSISE00.env$plg[2,2] + p[5]*NRLMSISE00.env$plg[2,4] + p[28]*NRLMSISE00.env$plg[2,6] + t71)*NRLMSISE00.env$ctloc + 
                        (p[7]*NRLMSISE00.env$plg[2,2] + p[8]*NRLMSISE00.env$plg[2,4] + p[29]*NRLMSISE00.env$plg[2,6] + t72)*NRLMSISE00.env$stloc) 
    }
    # semidiurnal
    if (flags$sw[9]) {
        t81 <- (p[24]*NRLMSISE00.env$plg[3,4]+p[36]*NRLMSISE00.env$plg[3,6])*cd14*flags$swc[6]
        t82 <- (p[34]*NRLMSISE00.env$plg[3,4]+p[37]*NRLMSISE00.env$plg[3,6])*cd14*flags$swc[6]
        t[8] <- f2*((p[6]*NRLMSISE00.env$plg[3,3]+ p[42]*NRLMSISE00.env$plg[3,5] + t81)*NRLMSISE00.env$c2tloc + 
                        (p[9]*NRLMSISE00.env$plg[3,3] + p[43]*NRLMSISE00.env$plg[3,5] + t82)*NRLMSISE00.env$s2tloc)
    }
    # terdiurnal
    if (flags$sw[9]) {
        t[14] <- f2 * ((p[40]*NRLMSISE00.env$plg[4,4]+(p[94]*NRLMSISE00.env$plg[4,5]+p[47]*NRLMSISE00.env$plg[4,7])*cd14*flags$swc[6])* NRLMSISE00.env$s3tloc +
                           (p[41]*NRLMSISE00.env$plg[4,4]+(p[95]*NRLMSISE00.env$plg[4,5]+p[49]*NRLMSISE00.env$plg[4,7])*cd14*flags$swc[6])*NRLMSISE00.env$c3tloc)
    }
    # magnetic activity based on daily ap
    if (flags$sw[10] == -1) {
        ap <- input$ap_a
        if (p[52] != 0) {
            exp1 <- exp(-10800*sqrt(p[52]*p[52])/(1+p[139]*(45-abs(input$g_lat))))
            exp1 <- min(exp1, 0.99999)
            p[25] <- max(p[25], 1e-4)
            assign("apt", sg0(exp1, p, ap), NRLMSISE00.env)
            if(flags$sw[10]) {
                t[9] <- NRLMSISE00.env$apt*(p[51] + p[97]*NRLMSISE00.env$plg[1,3] + p[55]*NRLMSISE00.env$plg[1,5] +  
                                    (p[126]*NRLMSISE00.env$plg[1,2] + p[127]*NRLMSISE00.env$plg[1,4] + p[128]*NRLMSISE00.env$plg[1,6])*cd14*flags$swc[6] +  
                                    (p[129]*NRLMSISE00.env$plg[2,2] + p[130]*NRLMSISE00.env$plg[2,4] + p[131]*NRLMSISE00.env$plg[2,6])*flags$swc[8]* 
                                    cos(hr*(tloc - p[132])))
            }
        }
    } else {
        apd <- input$ap - 4
        p44 <- p[44]
        p45 <- p[45]
        if (p44<0) {
            p44 <- 1E-5
        }
        assign("apdf", apd + (p45-1)*(apd + (exp(-p44 * apd) - 1)/p44), NRLMSISE00.env)
        if (flags$sw[10])
            t[9] <- NRLMSISE00.env$apdf*(p[33] + p[46]*NRLMSISE00.env$plg[1,3] + p[35]*NRLMSISE00.env$plg[1,5] +  
                       (p[101]*NRLMSISE00.env$plg[1,2] + p[102]*NRLMSISE00.env$plg[1,4] + p[103]*NRLMSISE00.env$plg[1,6])*cd14*flags$swc[6] +  
                       (p[122]*NRLMSISE00.env$plg[2,2] + p[123]*NRLMSISE00.env$plg[2,4] + p[124]*NRLMSISE00.env$plg[2,6])*flags$swc[8]* 
                       cos(hr*(tloc-p[125])))
    }
    if (((flags$sw[11]) & (input$g_long > -1000))) {
        # longitudinal
        if (flags$sw[12]) {
            t[11] <- (1 + p[81]*NRLMSISE00.env$dfa*flags$swc[2])*((p[65]*NRLMSISE00.env$plg[2,3]+p[66]*NRLMSISE00.env$plg[2,5]+p[67]*NRLMSISE00.env$plg[2,7] +
                        p[104]*NRLMSISE00.env$plg[2,2]+p[105]*NRLMSISE00.env$plg[2,4]+p[106]*NRLMSISE00.env$plg[2,6] +
                        flags$swc[6]*(p[110]*NRLMSISE00.env$plg[2,2]+p[111]*NRLMSISE00.env$plg[2,4] + 
                                          p[112]*NRLMSISE00.env$plg[2,6])*cd14)* cos(dgtr*input$g_long)
                        + (p[91]*NRLMSISE00.env$plg[2,3]+p[92]*NRLMSISE00.env$plg[2,5]+p[93]*NRLMSISE00.env$plg[2,7] + 
                               p[107]*NRLMSISE00.env$plg[2,2]+p[108]*NRLMSISE00.env$plg[2,4]+p[109]*NRLMSISE00.env$plg[2,6] + 
                               flags$swc[6]*(p[113]*NRLMSISE00.env$plg[2,2]+p[114]*NRLMSISE00.env$plg[2,4] + 
                                                 p[115]*NRLMSISE00.env$plg[2,6])*cd14)*sin(dgtr*input$g_long))
        }
        # ut and mixed ut, longitude
        if (flags$sw[13]) {
            t[12] <- (1+p[96]*NRLMSISE00.env$plg[1,2])*(1+p[82]*NRLMSISE00.env$dfa*flags$swc[2])*
                (1+p[120]*NRLMSISE00.env$plg[1,2]*flags$swc[6]*cd14)*
                ((p[69]*NRLMSISE00.env$plg[1,2]+p[70]*NRLMSISE00.env$plg[1,4]+p[71]*NRLMSISE00.env$plg[1,6])*
                     cos(sr*(input$sec - p[72])))
            t[12] <- t[12] + flags$swc[12]* (p[77]*NRLMSISE00.env$plg[3,4]+p[78]*
                NRLMSISE00.env$plg[3,6]+p[79]*NRLMSISE00.env$plg[3,8])* cos(sr*(input$sec-p[80])+2*dgtr*input$g_long) *
                (1+p[138]*NRLMSISE00.env$dfa*flags$swc[2])
        }
        # ut, longitude magnetic activity
        if (flags$sw[14]) {
            if (flags$sw[10] == -1) {
                if(p[52]) {
                    t[13] <- NRLMSISE00.env$apt*flags$swc[12]*(1.+p[133]*NRLMSISE00.env$plg[1,2])*
                        ((p[53]*NRLMSISE00.env$plg[2,3]+p[99]*NRLMSISE00.env$plg[2,5]+p[68]*NRLMSISE00.env$plg[2,7])*
                        cos(dgtr*(input$g_long-p[98])))+ NRLMSISE00.env$apt*flags$swc[12]*
                        flags$swc[6]*(p[134]*NRLMSISE00.env$plg[2,2]+p[135]*NRLMSISE00.env$plg[2,4]+p[136]*NRLMSISE00.env$plg[2,6])*
                        cd14*cos(dgtr*(input$g_long-p[137]))+NRLMSISE00.env$apt*flags$swc[13]*
                        (p[56]*NRLMSISE00.env$plg[1,2]+p[57]*NRLMSISE00.env$plg[1,4]+p[58]*NRLMSISE00.env$plg[1,6])*cos(sr*(input$sec-p[59]))
                }
            } else {
                t[13] <- NRLMSISE00.env$apdf*flags$swc[12]*(1+p[121]*NRLMSISE00.env$plg[1,2])*
                    ((p[61]*NRLMSISE00.env$plg[2,3]+p[62]*NRLMSISE00.env$plg[2,5]+p[63]*NRLMSISE00.env$plg[2,7])*
                         cos(dgtr*(input$g_long-p[64]))) + NRLMSISE00.env$apdf*flags$swc[12]*flags$swc[6]* 
                    (p[116]*NRLMSISE00.env$plg[2,2]+p[117]*NRLMSISE00.env$plg[2,4]+p[118]*NRLMSISE00.env$plg[2,6])* cd14*
                    cos(dgtr*(input$g_long-p[119])) + NRLMSISE00.env$apdf*flags$swc[13]*
                    (p[84]*NRLMSISE00.env$plg[1,2]+p[85]*NRLMSISE00.env$plg[1,4]+p[86]*NRLMSISE00.env$plg[1,6])*
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
    t[1] <- p[22] * NRLMSISE00.env$dfa
    # time independent
    t[2] <- p[2]*NRLMSISE00.env$plg[1,3] + p[3]*NRLMSISE00.env$plg[1,5] + p[23]*NRLMSISE00.env$plg[1,7] + p[27]*NRLMSISE00.env$plg[1,2] + p[15]*NRLMSISE00.env$plg[1,4] + p[60]*NRLMSISE00.env$plg[1,6]
    # symmetrical annual
    t[3] <- (p[19]+p[48]*NRLMSISE00.env$plg[1,3]+p[30]*NRLMSISE00.env$plg[1,5])*cd32
    # symmetrical semiannual
    t[4] <- (p[16]+p[17]*NRLMSISE00.env$plg[1,3]+p[31]*NRLMSISE00.env$plg[1,5])*cd18
    # asymmetrical annual
    t[5] <- (p[10]*NRLMSISE00.env$plg[1,2]+p[11]*NRLMSISE00.env$plg[1,4]+p[21]*NRLMSISE00.env$plg[1,6])*cd14
    # asymmetrical semiannual
    t[6] <- (p[38]*NRLMSISE00.env$plg[1,2])*cd39
    # diurnal
    if (flags$sw[8]) {
        t71 <- p[12]*NRLMSISE00.env$plg[2,3]*cd14*flags$swc[6]
        t72 <- p[13]*NRLMSISE00.env$plg[2,3]*cd14*flags$swc[6]
        t[7] <- ((p[4]*NRLMSISE00.env$plg[2,2] + p[5]*NRLMSISE00.env$plg[2,4] + t71) * NRLMSISE00.env$ctloc + (p[7]*NRLMSISE00.env$plg[2,2] + p[8]*NRLMSISE00.env$plg[2,4] + t72) * NRLMSISE00.env$stloc)
    }
    # semidiurnal
    if (flags$sw[9]) {
        t81 <- (p[24]*NRLMSISE00.env$plg[3,4]+p[36]*NRLMSISE00.env$plg[3,6])*cd14*flags$swc[6]
        t82 <- (p[34]*NRLMSISE00.env$plg[3,4]+p[37]*NRLMSISE00.env$plg[3,6])*cd14*flags$swc[6]
        t[8] <- ((p[6]*NRLMSISE00.env$plg[3,3] + p[42]*NRLMSISE00.env$plg[3,5] + t81) * NRLMSISE00.env$c2tloc + (p[9]*NRLMSISE00.env$plg[3,3] + p[43]*NRLMSISE00.env$plg[3,5] + t82) * NRLMSISE00.env$s2tloc)
    }
    # terdiurnal
    if (flags$sw[15]) {
        t[14] <- p[40] * NRLMSISE00.env$plg[4,4] * NRLMSISE00.env$s3tloc + p[41] * NRLMSISE00.env$plg[4,4] * NRLMSISE00.env$c3tloc
    }
    # magnetic activity
    if (flags$sw[10]) {
        if (flags$sw[10] == 1) {
            t[9] <- NRLMSISE00.env$apdf * (p[33] + p[46] * NRLMSISE00.env$plg[1, 3] * flags$swc[3])
        } else if (flags$sw[10] == -1) {
            t[9] <- (p[51]*NRLMSISE00.env$apt + p[97]*NRLMSISE00.env$plg[1,3] * NRLMSISE00.env$apt*flags$swc[3])
        }
    }
    # longitudinal
    if (!((flags$sw[11] == 0) | (flags$sw[12] == 0) | (input$g_long <= -1000))) {
        t[11] <- (1 + NRLMSISE00.env$plg[1,2]*(p[81]*flags$swc[6]*cos(dr*(input$doy-p[82]))
                                + p[86]*flags$swc[7]*cos(2*dr*(input$doy-p[87])))
                  + p[84]*flags$swc[4]*cos(dr*(input$doy-p[85]))
                  + p[88]*flags$swc[5]*cos(2*dr*(input$doy-p[89]))) *
         ((p[65]*NRLMSISE00.env$plg[2,3]+p[66]*NRLMSISE00.env$plg[2,5]+p[67]*NRLMSISE00.env$plg[2,7] +
                p[75]*NRLMSISE00.env$plg[2,2]+p[76]*NRLMSISE00.env$plg[2,4]+p[77]*NRLMSISE00.env$plg[2,6]) *
               cos(dgtr*input$g_long) + (p[91]*NRLMSISE00.env$plg[2,3]+p[92]*NRLMSISE00.env$plg[2,5]+p[93]*NRLMSISE00.env$plg[2,7] +
                                             p[78]*NRLMSISE00.env$plg[2,2]+p[79]*NRLMSISE00.env$plg[2,4]+p[80]*NRLMSISE00.env$plg[2,6])
           * sin(dgtr*input$g_long))
    }
    tt <- 0
    for (i in 1:14) {
        tt <- tt + abs(flags$sw[i+1])*t[[i]]
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
    za <- asteRiskData::pdl[2, 16]
    zn1[1] <- za
    output <- list(
        d = numeric(9),
        t = numeric(2)
    )
    if (input$alt > zn1[1]) { # tinf variations not important below za or zn1[1]
        tinf <- asteRiskData::ptm[1] * asteRiskData::pt[1] * (1 + flags$sw[17] * globe7(asteRiskData::pt,input,flags))
    } else {
        tinf <- asteRiskData::ptm[1] * asteRiskData::pt[1]
    }
    output$t[1] <- tinf
    if (input$alt > zn1[5]) { # gradient variations not important below zn1[5]
        g0 <- asteRiskData::ptm[4] * asteRiskData::ps[1] * (1 + flags$sw[20] * globe7(asteRiskData::ps,input,flags))
    } else {
        g0 <- asteRiskData::ptm[4] * asteRiskData::ps[1]
    }
    tlb <- asteRiskData::ptm[2] * (1 + flags$sw[18] * globe7(asteRiskData::pd[4, ], input, flags)) * asteRiskData::pd[4, 1]
    s <- g0 / (tinf - tlb)
    # Lower thermosphere temp variations not significant for density above 300 km
    if (input$alt < 300) {
        new_meso_tn1 <- c(0,
                          asteRiskData::ptm[7]*asteRiskData::ptl[1,1]/(1-flags$sw[19]*glob7s(asteRiskData::ptl[1,], input, flags)),
                          asteRiskData::ptm[3]*asteRiskData::ptl[2,1]/(1-flags$sw[19]*glob7s(asteRiskData::ptl[2,], input, flags)),
                          asteRiskData::ptm[8]*asteRiskData::ptl[3,1]/(1-flags$sw[19]*glob7s(asteRiskData::ptl[3,], input, flags)),
                          asteRiskData::ptm[5]*asteRiskData::ptl[4,1]/(1-flags$sw[19]*flags$sw[21]*glob7s(asteRiskData::ptl[4, ], input, flags)))
        assign("meso_tn1", new_meso_tn1, NRLMSISE00.env)
        assign("meso_tgn1", c(0, asteRiskData::ptm[9]*asteRiskData::pma[9,1]*(1+flags$sw[19]*flags$sw[21]*
                                                 glob7s(asteRiskData::pma[9, ], input, flags))* 
                   NRLMSISE00.env$meso_tn1[5]*NRLMSISE00.env$meso_tn1[5]/((asteRiskData::ptm[5]*asteRiskData::ptl[4,1])^2)),
               NRLMSISE00.env)
    } else {
        new_meso_tn1 <- c(0,
                          asteRiskData::ptm[7]*asteRiskData::ptl[1,1],
                          asteRiskData::ptm[3]*asteRiskData::ptl[2,1],
                          asteRiskData::ptm[8]*asteRiskData::ptl[3,1],
                          asteRiskData::ptm[5]*asteRiskData::ptl[4,1])
        assign("meso_tn1", new_meso_tn1, NRLMSISE00.env)
        assign("meso_tgn1", c(0, asteRiskData::ptm[9] * asteRiskData::pma[9,1] * NRLMSISE00.env$meso_tn1[5]*
                                  NRLMSISE00.env$meso_tn1[5]/((asteRiskData::ptm[5]*asteRiskData::ptl[4,1])^2)), 
               NRLMSISE00.env)
    }
    # N2 variation factor at Zlb
    g28 <- flags$sw[22] * globe7(asteRiskData::pd[3, ], input, flags)
    # variation of turbopause height
    zhf <- asteRiskData::pdl[2, 25] * (1 + flags$sw[6] * asteRiskData::pdl[1, 25] * sin(dgtr * input$g_lat) * 
                             cos(dr * (input$doy-asteRiskData::pt[14])))
    ### output$t[1] <- tinf
    xmm <- asteRiskData::pdm[3, 5]
    z <- input$alt
    # N2 density
    db28 <- asteRiskData::pdm[3,1] * exp(g28) * asteRiskData::pd[3, 1]
    densu_output <- densu(z,db28,tinf,tlb,28,alpha[3],
                          output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                          NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    output$d[3] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    assign("dd", output$d[3], NRLMSISE00.env)
    zh28 <- asteRiskData::pdm[3, 3] * zhf
    zhm28 <- asteRiskData::pdm[3, 4] * asteRiskData::pdl[2, 6] 
    xmd <- 28 - xmm
    densu_output <- densu(zh28,db28,tinf,tlb,xmd,(alpha[3]-1),
                          tz,asteRiskData::ptm[6],s,mn1, zn1,
                          NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    b28 <- densu_output$densu_temp
    tz <- densu_output$tz
    if ((flags$sw[16]) & (z <= altl[3])) {
        assign("dm28", densu(z,b28,tinf,tlb,xmm,alpha[3],
                             tz,asteRiskData::ptm[6],s,mn1,zn1,
                             NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)$densu_temp, 
               NRLMSISE00.env)
        output$d[3] <- dnet(output$d[3], NRLMSISE00.env$dm28, zhm28, xmm, 28)
    }
    # He density
    g4 <- flags$sw[22] * globe7(asteRiskData::pd[1, ], input, flags)
    db04 <- asteRiskData::pdm[1, 1] * exp(g4) * asteRiskData::pd[1,1]
    densu_output <- densu(z,db04,tinf,tlb, 4,alpha[1],
                          output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                          NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    output$d[1] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    assign("dd", output$d[1], NRLMSISE00.env)
    if ((flags$sw[16]) & ( z<altl[1]) ) {
        zh04 <- asteRiskData::pdm[1,3]
        densu_output <- densu(zh04,db04,tinf,tlb,4-xmm,alpha[1]-1.0,
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        b04 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b04,tinf,tlb,xmm,0,
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        assign("dm04", densu_output$densu_temp, NRLMSISE00.env)
        output$t[2] <- densu_output$tz
        zhm04 <- zhm28
        output$d[1] <- dnet(output$d[1],NRLMSISE00.env$dm04,zhm04,xmm,4)
        rl <- log(b28*asteRiskData::pdm[1,2]/b04)
        zc04 <- asteRiskData::pdm[1,5]*asteRiskData::pdl[2,1]
        hc04 <- asteRiskData::pdm[1,6]*asteRiskData::pdl[2,2]
        output$d[1] <- output$d[[1]]*ccor(z,rl,hc04,zc04)
    }
    # O density
    g16 <- flags$sw[22] * globe7(asteRiskData::pd[2, ], input, flags)
    db16 <- asteRiskData::pdm[2,1] * exp(g16) * asteRiskData::pd[2, 1]
    densu_output <- densu(z,db16,tinf,tlb, 16,alpha[2],
                         output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                         NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    output$d[2] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    assign("dd", output$d[2], NRLMSISE00.env)
    if ((flags$sw[16]) & ( z<=altl[2] )) {
        zh16 <- asteRiskData::pdm[2, 3]
        densu_output <- densu(zh16,db16,tinf,tlb,16-xmm,(alpha[2]-1),
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        b16 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b16,tinf,tlb,xmm,0,
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        assign("dm16", densu_output$densu_temp, NRLMSISE00.env)
        output$t[2] <- densu_output$tz
        zhm16 <- zhm28
        output$d[2] <- dnet(output$d[2],NRLMSISE00.env$dm16,zhm16,xmm,16)
        rl <- asteRiskData::pdm[2,2]*asteRiskData::pdl[2,17]*(1+flags$sw[2]*asteRiskData::pdl[1,24]*(input$f107A-150))
        hc16 <- asteRiskData::pdm[2,6]*asteRiskData::pdl[2,4]
        zc16 <- asteRiskData::pdm[2,5]*asteRiskData::pdl[2,3]
        hc216 <- asteRiskData::pdm[2,6]*asteRiskData::pdl[2,5]
        output$d[2] <- output$d[[2]]*ccor2(z,rl,hc16,zc16,hc216)
        hcc16 <- asteRiskData::pdm[2,8]*asteRiskData::pdl[2,14]
        zcc16 <- asteRiskData::pdm[2,7]*asteRiskData::pdl[2,13]
        rc16 <- asteRiskData::pdm[2,4]*asteRiskData::pdl[2,15]
        output$d[2] <- output$d[[2]]*ccor(z,rc16,hcc16,zcc16)
    }
    # O2 density
    g32 <- flags$sw[22]*globe7(asteRiskData::pd[5, ], input, flags)
    db32 <- asteRiskData::pdm[4,1] * exp(g32) * asteRiskData::pd[5,1]
    densu_output <- densu(z,db32,tinf,tlb,32,alpha[4],
                          output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                          NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    output$d[4] <- densu_output$densu_temp
    output$t[4] <- densu_output$tz
    assign("dd", output$d[4], NRLMSISE00.env)
    if (flags$sw[16]) {
        if (z <= altl[4]) {
            zh32 <- asteRiskData::pdm[4,3]
            densu_output <- densu(zh32,db32,tinf,tlb,32-xmm,alpha[4]-1,
                                  output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                                  NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
            b32 <- densu_output$densu_temp
            output$t[2] <- densu_output$tz
            densu_output <- densu(z,b32,tinf,tlb,xmm,0,output$t[2],
                                  asteRiskData::ptm[6],s,mn1,zn1,
                                  NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
            assign("dm32", densu_output$densu_temp, NRLMSISE00.env)
            output$t[2] <- densu_output$tz
            zhm32 <- zhm28
            output$d[4] <- dnet(output$d[4],NRLMSISE00.env$dm32,zhm32,xmm,32)
            rl <- log(b28*asteRiskData::pdm[4,2]/b32)
            hc32 <- asteRiskData::pdm[4,6]*asteRiskData::pdl[2,8]
            zc32 <- asteRiskData::pdm[4,5]*asteRiskData::pdl[2,7]
            output$d[4] <- output$d[[4]]*ccor(z,rl,hc32,zc32)
        }
        hcc32 <- asteRiskData::pdm[4,8]*asteRiskData::pdl[2,23]
        hcc232 <- asteRiskData::pdm[4,8]*asteRiskData::pdl[1,23]
        zcc32 <- asteRiskData::pdm[4,7]*asteRiskData::pdl[2,22]
        rc32 <- asteRiskData::pdm[4,4]*asteRiskData::pdl[2,24]*(1+flags$sw[2]*asteRiskData::pdl[1,24]*(input$f107A - 150))
        output$d[4] <- output$d[[4]]*ccor2(z,rc32,hcc32,zcc32,hcc232)
    }
    # Ar density
    g40 <- flags$sw[22] * globe7(asteRiskData::pd[6, ],input,flags)
    db40 <- asteRiskData::pdm[5,1]*exp(g40)*asteRiskData::pd[6, 1]
    densu_output <- densu(z,db40,tinf,tlb,40,alpha[5],
                          output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                          NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    output$d[5] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    assign("dd", output$d[5], NRLMSISE00.env)
    if ((flags$sw[16]) & (z <= altl[5])) {
        zh40 <- asteRiskData::pdm[5, 3]
        densu_output <- densu(zh40,db40,tinf,tlb,40-xmm,alpha[5]-1,
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        b40 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b40,tinf,tlb,xmm,0,
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        assign("dm40", densu_output$densu_temp, NRLMSISE00.env)
        output$t[2] <- densu_output$tz
        zhm40 <- zhm28
        output$d[5] <- dnet(output$d[5], NRLMSISE00.env$dm40, zhm40, xmm, 40)
        rl <- log(b28*asteRiskData::pdm[5,2]/b40)
        hc40 <- asteRiskData::pdm[5, 6] * asteRiskData::pdl[2, 10]
        zc40 <- asteRiskData::pdm[5, 5] * asteRiskData::pdl[2, 9]
        output$d[5] <- output$d[[5]] * ccor(z,rl,hc40,zc40)
    }
    # H density
    g1 <- flags$sw[22] * globe7(asteRiskData::pd[7, ], input, flags)
    db01 <- asteRiskData::pdm[6, 1] * exp(g1) * asteRiskData::pd[7, 1]
    densu_output <- densu(z,db01,tinf,tlb,1,alpha[7],
                          output$t[2],asteRiskData::ptm[6],s,mn1,zn1,NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    output$d[7] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    assign("dd", output$d[7], NRLMSISE00.env)
    if ((flags$sw[16]) & (z <= altl[7])) {
        zh01 <- asteRiskData::pdm[6, 3]
        densu_output <- densu(zh01,db01,tinf,tlb,1-xmm,alpha[7]-1,
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        b01 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b01,tinf,tlb,xmm,0,output$t[2],
                              asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        assign("dm01", densu_output$densu_temp, NRLMSISE00.env)
        output$t[2] <- densu_output$tz
        zhm01 <- zhm28
        output$d[7] <- dnet(output$d[7],NRLMSISE00.env$dm01,zhm01,xmm,1)
        rl <- log(b28*asteRiskData::pdm[6,2]*sqrt(asteRiskData::pdl[2,18]*asteRiskData::pdl[2,18])/b01)
        hc01 <- asteRiskData::pdm[6,6]*asteRiskData::pdl[2,12]
        zc01 <- asteRiskData::pdm[6,5]*asteRiskData::pdl[2,11]
        output$d[7] <- output$d[[7]]*ccor(z,rl,hc01,zc01)
        hcc01 <- asteRiskData::pdm[6,8]*asteRiskData::pdl[1,20]
        zcc01 <- asteRiskData::pdm[6,7]*asteRiskData::pdl[2,19]
        rc01 <- asteRiskData::pdm[6,4]*asteRiskData::pdl[2,21]
        output$d[7] <- output$d[[7]]*ccor(z,rc01,hcc01,zcc01)
    }
    # Atomic nitrogen (N) density
    g14 <- flags$sw[22]*globe7(asteRiskData::pd[8, ],input,flags)
    db14 <- asteRiskData::pdm[7,1]*exp(g14)*asteRiskData::pd[8,1]
    densu_output <- densu(z,db14,tinf,tlb,14,alpha[8],
                          output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                          NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    output$d[8] <- densu_output$densu_temp
    output$t[2] <- densu_output$tz
    assign("dd", output$d[8], NRLMSISE00.env)
    if ((flags$sw[16]) & (z <= altl[8])) {
        zh14 <- asteRiskData::pdm[7, 3]
        densu_output <- densu(zh14,db14,tinf,tlb,14-xmm,alpha[8]-1,
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        b14 <- densu_output$densu_temp
        output$t[2] <- densu_output$tz
        densu_output <- densu(z,b14,tinf,tlb,xmm,0,
                              output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                              NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
        assign("dm14", densu_output$densu_temp, NRLMSISE00.env)
        output$t[2] <- densu_output$tz
        zhm14 <- zhm28
        output$d[8] <- dnet(output$d[8],NRLMSISE00.env$dm14,zhm14,xmm,14)
        rl <- log(b28*asteRiskData::pdm[7,2]*sqrt(asteRiskData::pdl[1,3]*asteRiskData::pdl[1,3])/b14)
        hc14 <- asteRiskData::pdm[7,6]*asteRiskData::pdl[1,2]
        zc14 <- asteRiskData::pdm[7,5]*asteRiskData::pdl[1,1]
        output$d[8] <- output$d[[8]]*ccor(z,rl,hc14,zc14)
        hcc14 <- asteRiskData::pdm[7,8]*asteRiskData::pdl[1,5]
        zcc14 <- asteRiskData::pdm[7,7]*asteRiskData::pdl[1,4]
        rc14 <- asteRiskData::pdm[7,4]*asteRiskData::pdl[1,6]
        output$d[8] <- output$d[[8]]*ccor(z,rc14,hcc14,zcc14)
    }
    # Anomalous oxygen density
    g16h <- flags$sw[22]*globe7(asteRiskData::pd[9, ],input,flags)
    db16h <- asteRiskData::pdm[8,1]*exp(g16h)*asteRiskData::pd[9,1]
    tho <- asteRiskData::pdm[8,10]*asteRiskData::pdl[1,7]
    densu_output <- densu(z,db16h,tho,tho,16,alpha[9],
                          output$t[2],asteRiskData::ptm[6],s,mn1,zn1,
                          NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    assign("dd", densu_output$densu_temp, NRLMSISE00.env)
    output$t[2] <- densu_output$tz
    zsht <- asteRiskData::pdm[8,6]
    zmho <- asteRiskData::pdm[8,5]
    zsho <- scalh(zmho,16,tho,NRLMSISE00.env$gsurf)
    output$d[9] <- NRLMSISE00.env$dd*exp(-zsht/zsho*(exp(-(z-zmho)/zsht)-1))
    # Total density (atomic mass-weighted sums)
    output$d[6] <- 1.66e-24*(4*output$d[[1]]+16*output$d[[2]]+28*output$d[[3]]+32*output$d[[4]]+40*output$d[[5]]+ output$d[[7]]+14*output$d[[8]])
    # Temperature
    z <- abs(input$alt)
    densu_output <- densu(z,1,tinf,tlb,0,0,
                          output$t[2],asteRiskData::ptm[6],s,mn1,zn1, 
                          NRLMSISE00.env$meso_tn1,NRLMSISE00.env$meso_tgn1)
    output$t[2] <- densu_output$tz
    if (flags$sw[1]) {
        for(i in 1:9) {
            output$d[i] <- output$d[[i]]*1e6;
        }
        output$d[6] <- output$d[[6]]/1000
    }
    return(output)
}

gtd7 <- function(input, flags) {
    for (i in 1:24) {
        if (i == 10) {
            flags$sw[i] <- flags$switches[i]
            flags$swc[i] <- flags$switches[i]
        } else {
            flags$sw[i] <- if(flags$switches[i] == 1) 1 else 0
            flags$swc[i] <- if(flags$switches[i] > 0) 1 else 0
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
    assign("gsurf", glatf_results[[1]], NRLMSISE00.env)
    assign("re_", glatf_results[[2]], NRLMSISE00.env)
    xmm <- asteRiskData::pdm[3,5]
    altt <- max(input$alt, zn2[1])
    tmp <- input$alt
    input$alt <- altt
    soutput <- gts7(input, flags)
    input$alt <- tmp
    if (flags$sw[1]) {
        dm28m <- NRLMSISE00.env$dm28 * 1e6
    } else {
        dm28m <- NRLMSISE00.env$dm28
    }
    output$t <- soutput$t
    if (input$alt >= zn2[1]) {
        output$d <- soutput$d
    } else { # lower mesosphere/upper stratosphere (between zn3[1] and zn2[1])
        new_meso_tn2 <- c(NRLMSISE00.env$meso_tn1[5],
                          asteRiskData::pma[1,1]*asteRiskData::pavgm[1]/(1-flags$sw[21]*glob7s(asteRiskData::pma[1, ], input, flags)),
                          asteRiskData::pma[2,1]*asteRiskData::pavgm[2]/(1-flags$sw[21]*glob7s(asteRiskData::pma[2, ], input, flags)),
                          asteRiskData::pma[3,1]*asteRiskData::pavgm[3]/(1-flags$sw[21]*flags$sw[23]*glob7s(asteRiskData::pma[3, ], input, flags)))
        assign("meso_tn2", new_meso_tn2, NRLMSISE00.env)
        new_meso_tgn2 <-
            c(
                NRLMSISE00.env$meso_tgn1[2],
                asteRiskData::pavgm[9] * asteRiskData::pma[10, 1] * (1 + flags$sw[21] * flags$sw[23] *
                                             glob7s(asteRiskData::pma[10,], input, flags)) *
                    NRLMSISE00.env$meso_tn2[4] * NRLMSISE00.env$meso_tn2[4] /
                    ((asteRiskData::pma[3, 1] * asteRiskData::pavgm[3]) ^ 2)
            )
        assign("meso_tgn2", new_meso_tgn2, NRLMSISE00.env)
        new_meso_tn3 <- c(NRLMSISE00.env$meso_tn2[4], NRLMSISE00.env$meso_tn3[2:5])
        if (input$alt < zn3[1]) { # lower stratosphere/troposphere (below zn3[1])
            new_meso_tn3[2:5] <- c(asteRiskData::pma[4,1]*asteRiskData::pavgm[4]/(1-flags $sw[23]*glob7s(asteRiskData::pma[4, ], input, flags)),
                                   asteRiskData::pma[5,1]*asteRiskData::pavgm[5]/(1-flags $sw[23]*glob7s(asteRiskData::pma[5, ], input, flags)),
                                   asteRiskData::pma[6,1]*asteRiskData::pavgm[6]/(1-flags $sw[23]*glob7s(asteRiskData::pma[6, ], input, flags)),
                                   asteRiskData::pma[7,1]*asteRiskData::pavgm[7]/(1-flags $sw[23]*glob7s(asteRiskData::pma[7, ], input, flags)))
            assign("meso_tn3", new_meso_tn3, NRLMSISE00.env)
            new_meso_tgn3 <- c(NRLMSISE00.env$meso_tgn2[2],
                               asteRiskData::pma[8,1]*asteRiskData::pavgm[8]*(1+flags $sw[23]* glob7s(asteRiskData::pma[8, ], input, flags)) *
                                   NRLMSISE00.env$meso_tn3[5]*NRLMSISE00.env$meso_tn3[5]/
                                   ((asteRiskData::pma[7,1]*asteRiskData::pavgm[7])^2))
            assign("meso_tgn3", new_meso_tgn3, NRLMSISE00.env)
        } else {
            assign("meso_tn3", new_meso_tn3, NRLMSISE00.env)
        }
        # linear transition to full mixing below zn2[1]
        dmc <- 0
        if (input$alt > zmix) {
            dmc <- 1 - (zn2[1] - input$alt)/(zn2[1] - zmix)
        }
        dz28 <- soutput$d[3]
        ## N2 density 
        dmr <- soutput$d[[3]] / NRLMSISE00.env$dm28m - 1
        densm_results <- densm(input$alt, NRLMSISE00.env$dm28m, xmm, tz, mn3, zn3, 
                               NRLMSISE00.env$meso_tn3, NRLMSISE00.env$meso_tgn3, mn2, zn2, 
                               NRLMSISE00.env$meso_tn2, NRLMSISE00.env$meso_tgn2)
        output$d[3] <- densm_results$densm_tmp
        tz < densm_results$tz
        output$d[3] <- output$d[[3]] * (1 + dmr*dmc)
        ## He density
        dmr <- soutput$d[1] / (dz28 * asteRiskData::pdm[1,2]) - 1
        output$d[1] <- output$d[[3]] * asteRiskData::pdm[1,2] * (1 + dmr*dmc)
        ## O density
        output$d[c(2, 9)] <- 0
        ## O2 density
        dmr <- soutput$d[[4]] / (dz28 * asteRiskData::pdm[4,2]) - 1
        output$d[4] <- output$d[[3]] * asteRiskData::pdm[4,2] * (1 + dmr*dmc)
        # Ar density
        dmr <- soutput$d[5] / (dz28 * asteRiskData::pdm[5,2]) - 1
        output$d[5] <- output$d[[3]] * asteRiskData::pdm[5,2] * (1 + dmr*dmc)
        # H density
        output$d[7] <- 0
        # Atomic nitrogen (N) density
        output$d[8] <- 0
        # Total mass density
        output$d[6] <- 1.66e-24 * (4 * output$d[[1]] + 16 * output$d[[2]] + 
                                   28 * output$d[[3]] + 32 * output$d[[4]] + 40 * 
                                   output$d[[5]] + output$d[[7]] + 14 * output$d[[8]])
        if (flags$sw[1]) {
            output$d[6] <- output$d[[6]]/1000
        }
        densm_results <- densm(input$alt, 1, 0, tz, mn3, zn3, 
                               NRLMSISE00.env$meso_tn3, NRLMSISE00.env$meso_tgn3, mn2, zn2, 
                               NRLMSISE00.env$meso_tn2, NRLMSISE00.env$meso_tgn2)
        assign("dd", densm_results$densm_temp, NRLMSISE00.env)
        tz <- densm_results$tz
        output$t[2] <- tz
    }
    return(output)
}

gtd7d <- function(input, flags) {
    output <- gtd7(input, flags)
    output$d[6] <- 1.66e-24 * (4 * output$d[[1]] + 16 * output$d[[2]] + 
                               28 * output$d[[3]] + 32 * output$d[[4]] + 
                               40 * output$d[[5]] + output$d[[7]] + 
                               14 * output$d[[8]] + 16 * output$d[[9]])
    if (flags$sw[1]) {
        output$d[6] <- output$d[[6]]/1000
    }
    return(output)
}

NRLMSISE00 <- function(Mjd_UTC, r_ECEF, UT1_UTC, TT_UTC) {
    hasData()
    assign("gsurf", 0, NRLMSISE00.env)
    assign("re_", asteRiskData::re, NRLMSISE00.env)
    assign("dd", 0, NRLMSISE00.env)
    assign("dm04", 0, NRLMSISE00.env)
    assign("dm16", 0, NRLMSISE00.env)
    assign("dm28", 0, NRLMSISE00.env)
    assign("dm32", 0, NRLMSISE00.env)
    assign("dm40", 0, NRLMSISE00.env)
    assign("dm01", 0, NRLMSISE00.env)
    assign("dm14", 0, NRLMSISE00.env)
    assign("meso_tn1", rep(0, 5), NRLMSISE00.env)
    assign("meso_tn2", rep(0, 4), NRLMSISE00.env)
    assign("meso_tn3", rep(0, 5), NRLMSISE00.env)
    assign("meso_tgn1", rep(0, 2), NRLMSISE00.env)
    assign("meso_tgn2", rep(0, 2), NRLMSISE00.env)
    assign("meso_tgn3", rep(0, 2), NRLMSISE00.env)
    assign("dfa", 0, NRLMSISE00.env)
    assign("plg", matrix(0, nrow=4, ncol=9), NRLMSISE00.env)
    assign("ctloc", 0, NRLMSISE00.env)
    assign("stloc", 0, NRLMSISE00.env)
    assign("c2tloc", 0, NRLMSISE00.env)
    assign("s2tloc", 0, NRLMSISE00.env)
    assign("c3tloc", 0, NRLMSISE00.env)
    assign("s3tloc", 0, NRLMSISE00.env)
    assign("apdf", 0, NRLMSISE00.env)
    assign("apt", 0, NRLMSISE00.env)
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
    geodetic_results <- ECEFtoLATLON(r_ECEF, degreesOutput=FALSE)
    input$alt <- geodetic_results["altitude"]/1000
    input$g_lat <- geodetic_results["latitude"]*(180/pi)
    input$g_long <- geodetic_results["longitude"]*(180/pi)
    iauCal2jd_results <- iauCal2jd(invjday_results$year, invjday_results$month, invjday_results$day)
    TIME <- (60*(60*invjday_results$hour + invjday_results$min) + invjday_results$sec)/86400
    UTC <- iauCal2jd_results$DATE + TIME
    TT <- UTC + TT_UTC/86400
    TUT <- TIME + UT1_UTC/86400
    UT1 <- iauCal2jd_results$DATE + TUT
    rnpb <- iauPnm06a(iauCal2jd_results$DJMJD0, TT)
    lst <- geodetic_results["longitude"] + iauGst06(iauCal2jd_results$DJMJD0, UT1, 
                                           iauCal2jd_results$DJMJD0, TT, rnpb)
    lst <- lst %% (2*pi)
    lst <- (lst*24)/(2*pi)
    print("ESTO es LST")
    print(list)
    input$lst <- lst
    i <- which(((invjday_results$year == asteRiskData::spaceWeather[, 1]) & 
                    (invjday_results$mon == asteRiskData::spaceWeather[, 2]) & 
                    (invjday_results$day == asteRiskData::spaceWeather[, 3])))[1]
    sw <- as.numeric(asteRiskData::spaceWeather[i, ])
    input$ap <- sw[23]
    input$ap_a[1] <- sw[23] 
    input$ap_a[2] <- sw[15] 
    sw_1 <- as.numeric(asteRiskData::spaceWeather[i-1, ])
    input$ap_a[3] <- sw_1[22] 
    input$ap_a[4] <- sw_1[21] 
    input$ap_a[5] <- sw_1[20] 
    sum <- sw_1[19]+sw_1[18]+sw_1[17]+sw_1[16]+sw_1[15]
    sw_2 <- as.numeric(asteRiskData::spaceWeather[i-2, ])
    sum <- sum + sw_2[22] + sw_2[21]  +sw_2[20]
    input$ap_a[6] <- sum/8 
    sw_3 <- as.numeric(asteRiskData::spaceWeather[i-3, ])
    sum <- sw_2[19] + sw_2[18] + sw_2[17] + sw_2[16] + 
          sw_2[15] + sw_3[22] + sw_3[21] + sw_3[20]
    input$ap_a[7] <- sum/8
    input$f107 <- sw[31]  
    input$f107A <- sw[33]   
    output <- gtd7d(input, flags)
    output_densities <- output$d
    #dens <- 1e3*output_densities[[6]]
    output_densities[6] <- 1e3*output_densities[6]
    names(output_densities) <- c("He", "O", "N2", "O2", "Ar", "Total", "H",
                                 "N", "AnomalousO")
    return(output_densities)
}

NRLMSISE00model <- function(position_ECEF, dateTime) {
    date <- strptime(dateTime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
    year <- date$year + 1900
    month <- date$mon + 1
    day <- date$mday
    hour <- date$hour
    minute <- date$min
    second <- date$sec
    Mjd_UTC <- iauCal2jd(year, month, day, hour, minute, second)$DATE
    IERS_results <- IERS(asteRiskData::earthPositions, Mjd_UTC, "l")
    UT1_UTC <- IERS_results$UT1_UTC[[1]]
    TAI_UTC <- IERS_results$TAI_UTC[[1]]
    timeDiffs_results <- timeDiffs(UT1_UTC, TAI_UTC)
    TT_UTC <- timeDiffs_results$TT_UTC[[1]]
    densities <- NRLMSISE00(Mjd_UTC, position_ECEF, UT1_UTC, TT_UTC)
    return(densities)
}