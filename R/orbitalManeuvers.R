lambert2F0elliptic <- function(x, y, N, tau) {
    acot1 <- acot(x/sqrt(1-x^2))
    if(acot1 < 0) {
        acot1 <- acot1 + pi
    }
    acot2 <- acot(y/sqrt(1-y^2))
    return((acot1 - acot2 - x*sqrt(1-x^2) + y*sqrt(1-y^2) + N*pi)/(sqrt((1-x^2)^3)) - tau)
}

lambert2F1elliptic <- function(x, y, N, tau, sigma) {
    return((3*x*lambert2F0elliptic(x, y, N, 0) - 2*(1-sigma^3*x/abs(y)))/(1-x^2))
}

lambert2F2elliptic <- function(x, y, N, tau, sigma) {
    return(((1+4*x^2)*lambert2F1elliptic(x, y, N, 0, sigma) + 2*(1-sigma^5*x^3/abs(y)^3))/(x*(1-x^2)))
}

lambert2F0hyperbolic <- function(x, y, tau) {
    acoth1 <- acoth(x/sqrt(x^2-1))
    acoth2 <- acoth(y/sqrt(y^2-1))
    return((-acoth1 + acoth2 + x*sqrt(x^2 - 1) - y*sqrt(y^2 - 1))/(sqrt((x^2 - 1)^3)) - tau)
}

lambert2F1hyperbolic <- function(x, y, tau, sigma) {
    return((3*x*lambert2F0hyperbolic(x, y, 0) - 2*(1-sigma^3*x/abs(y)))/(1-x^2))
}

lambert2F2hyperbolic <- function(x, y, tau, sigma) {
    return(((1+4*x^2)*lambert2F1hyperbolic(x, y, 0, sigma) + 2*(1-sigma^5*x^3/abs(y)^3))/(x*(1-x^2)))
    # return((3*lambert2F0hyperbolic(x, y, 0) + 5*x*lambert2F1hyperbolic(x, y, 0, sigma) + 
    #         2*(1-sigma^2)*(sigma^3/y^3))/(1-x^2))
}

lambert2Phi0 <- function(x, y, N) {
    acotX <- acot(x/sqrt(1-x^2))
    if(acotX < 0) {
        acotX <- acotX + pi
    }
    acotY <- acot(y/sqrt(1-y^2))
    phiX <- acotX - (1/(3*x)*(2+x^2)*sqrt(1-x^2))
    phiY <- acotY - (1/(3*y)*(2+y^2)*sqrt(1-y^2))
    #return(phiY - phiX - N*pi)
    return(phiX - phiY + N*pi)
}

lambert2Phi1 <- function(x, y, sigma) {
    return(
        (2/3) * (((1-x^2)^(3/2))/(x^2)) * (1-(sigma^5*x^3)/abs(y)^3)
    )
}

lambert2Phi2 <- function(x, y, sigma) {
    return(
        (2/3) * sqrt(1-x^2) * (
            (-2-x^2)/x^3 * (
                1 - (sigma^5*x^3)/abs(y)^3
            ) - (3*sigma^5*(1-sigma^2)*(1-x^2))/abs(y)^5
        )
    )
}

lambert2 <- function(initialPosition, finalPosition, initialTime, finalTime,
                     retrogradeTransfer=FALSE, centralBody="Earth", maxIterations=100) {
    checkTimeInput(initialTime, finalTime, "sdp4")
    t0 <- 0
    if(is.character(initialTime) & is.character(finalTime)) {
        t <- as.numeric(difftime(initialTime, finalTime, units="secs"))
    } else {
        t <- finalTime
    }
    centralBody <- tolower(centralBody)
    if(!(centralBody %in% c("sun", "mercury", "venus", "earth", "moon", "mars", "jupiter",
                            "saturn", "uranus", "neptune"))) {
        stop("Only the following central bodies are currently supported:
             Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, 
             Neptune and Pluto")
    }
    GM <- switch(centralBody, sun=GM_Sun_DE440, mercury=GM_Mercury_DE440,
                 venus=GM_Venus_DE440, earth=GM_Earth_DE440, moon=GM_Moon_DE440,
                 mars=GM_Mars_DE440, jupiter=GM_Jupiter_DE440, saturn=GM_Saturn_DE440,
                 uranus=GM_Uranus_DE440, neptune=GM_Neptune_DE440, pluto=GM_Pluto_DE440)
    initialDistance <- sqrt(sum(initialPosition^2))
    finalDistance <- sqrt(sum(finalPosition^2))
    diffVector <- initialPosition - finalPosition
    diffDistance <- sqrt(sum(diffVector^2))
    transferAngle <- acos(sum(initialPosition*finalPosition)/
                              (initialDistance*finalDistance))
    alpha <- vectorCrossProduct3D(initialPosition, finalPosition)[3]
    if((!retrogradeTransfer & alpha<0) | (retrogradeTransfer & alpha>0)) {
        transferAngle <- 2*pi - transferAngle
    }
    m <- initialDistance + finalDistance + diffDistance
    n <- initialDistance + finalDistance - diffDistance
    sigma2 <- (4*initialDistance*finalDistance*cos(transferAngle/2)^2)/m^2
    if(transferAngle < pi) {
        sigma <- sqrt(sigma2)
    } else {
        sigma <- -sqrt(sigma2)
    }
    tau <- 4*t*sqrt(GM/m^3)
    tauParabolic <- 2/3*(1-sigma^3)
    diffTauParabolic <- tau - tauParabolic
    if(abs(diffTauParabolic) <= .Machine$double.eps) {
        orbitType <- "parabolic"
    } else if(diffTauParabolic > 0) {
        orbitType <- "elliptic"
    } else if(diffTauParabolic < 0) {
        orbitType <- "hyperbolic"
    }
    if(orbitType =="elliptic") {
        Nmax <- floor(tau/pi)
        solutionsList <- vector(mode="list", length=1 + Nmax*2)
        tauME <- acos(sigma) + sigma*sqrt(1-sigma^2)
        if(abs(tauME - tau) <= .Machine$double.eps) {
            x1 <- 0
            path <- "Minimum Energy"
        } else if(tau < tauME) {
            x1 <- 0.5
            path <- "Low"
        } else {
            x1 <- -0.5
            path <- "High"
        }
        x <- x1
        y <- sqrt(1-sigma2*(1-x^2))
        if(sigma < 0) {
            y <- -y
        }
        F0 <- lambert2F0elliptic(x, y, 0, tau) - tau
        F1 <- lambert2F1elliptic(x, y, 0, tau, sigma)
        F2 <- lambert2F2elliptic(x, y, 0, tau, sigma)
        signF1 <- sign(F1)
        degree <- 2
        convergence <- FALSE
        numberIterations <- 0
        while(!convergence) {
            xNew <- x - degree*F0/(F1 + F1/abs(F1)*sqrt((degree-1)^2*F1^2 - degree*(degree-1)*F0*F2))
            difference <- xNew - x
            if(abs(difference) < .Machine$double.eps/2) {
                convergence <- TRUE
            }
            x <- xNew
            y <- sqrt(1-sigma2*(1-x^2)) * sign(sigma)
            F0 <- lambert2F0elliptic(x, y, 0, tau)
            F1 <- lambert2F1elliptic(x, y, tau, 0, sigma)
            F2 <- lambert2F2elliptic(x, y, tau, 0, sigma)
            numberIterations <- numberIterations + 1
            if(numberIterations > maxIterations) {
                degree <- degree + 1
                numberIterations <- 0
            }
        }
        vc <- sqrt(GM) * (y/sqrt(n) + x/sqrt(m))
        vr <- sqrt(GM) * (y/sqrt(n) - x/sqrt(m))
        ec <- (finalPosition - initialPosition)/diffDistance
        er1 <- initialPosition/sqrt(sum(initialPosition^2))
        er2 <- finalPosition/sqrt(sum(finalPosition^2))
        v1 <- vc*ec + vr*er1
        v2 <- vc*ec - vr*er2
        solutionsList[[1]] <- list(numberRevs = 0, path=path, 
                                   orbitType=orbitType, v1=v1, v2=v2)
        if(Nmax >= 1) {
            for(nrev in 1:Nmax) {
                #tauMT <- 2/3*()
                tauME <- nrev*pi + acos(sigma) + sigma*sqrt(1-sigma^2)
                # first get tauMT
                convergence <- FALSE
                numberIterations <- 0
                xMT <- 0.14
                yMT <- sqrt(1-sigma2*(1-xMT^2)) * sign(sigma)
                phi0 <- lambert2Phi0(xMT, yMT, nrev)
                phi1 <- lambert2Phi1(xMT, yMT, sigma)
                phi2 <- lambert2Phi2(xMT, yMT, sigma)
                while(!convergence) {
                    xMTNew <- xMT - degree*phi0/(phi1 + phi1/abs(phi1)*sqrt((degree-1)^2*phi1^2 - degree*(degree-1)*phi0*phi2))
                    difference <- xMTNew - xMT
                    if(abs(difference) < .Machine$double.eps/2) {
                        convergence <- TRUE
                    }
                    xMT <- xMTNew
                    yMT <- sqrt(1-sigma2*(1-xMT^2)) * sign(sigma)
                    phi0 <- lambert2Phi0(xMT, yMT, nrev)
                    phi1 <- lambert2Phi1(xMT, yMT, sigma)
                    phi2 <- lambert2Phi2(xMT, yMT, sigma)
                    numberIterations <- numberIterations + 1
                    if(numberIterations > maxIterations) {
                        degree <- degree + 1
                        numberIterations <- 0
                    }
                }
                # Currently using Laguerre's method for getting xMT, but the below block uses Newton
                # method. Seems to work well as well
                # numberIterations <- 0
                # while(!convergence) {
                #     xMTNew <- xMT - phi0/phi1
                #     difference <- xMTNew - xMT
                #     xMT <- xMTNew
                #     yMT <- sqrt(1-sigma2*(1-xMT^2)) * sign(sigma)
                #     phi0 <- lambert2Phi0(xMT, yMT, nrev)
                #     phi1 <- lambert2Phi1(xMT, yMT, sigma)
                #     #phi2 <- lambert2Phi2(xMT, yMT, sigma)
                #     numberIterations <- numberIterations + 1
                #     if(numberIterations > 99999 | abs(difference) < .Machine$double.eps/2) {
                #         convergence <- TRUE
                #     }
                # }
                tauMT <- (2/3)*(1/xMT - sigma^3/abs(yMT))
                convergence <- FALSE
                numberIterations <- 0
                # First the high path
                x <- -abs(x)/(nrev+1)
                y <- sqrt(1-sigma2*(1-x^2)) * sign(sigma)
                F0 <- lambert2F0elliptic(x, y, nrev, tau)
                F1 <- lambert2F1elliptic(x, y, tau, nrev, sigma)
                F2 <- lambert2F2elliptic(x, y, tau, nrev, sigma)
                while(!convergence) {
                    xNew <- x - degree*F0/(F1 + F1/abs(F1)*sqrt((degree-1)^2*F1^2 - degree*(degree-1)*F0*F2))
                    difference <- xNew - x
                    if(abs(difference) < .Machine$double.eps/2) {
                        convergence <- TRUE
                    }
                    x <- xNew
                    y <- sqrt(1-sigma2*(1-x^2)) * sign(sigma)
                    F0 <- lambert2F0elliptic(x, y, nrev, tau)
                    F1 <- lambert2F1elliptic(x, y, tau, nrev, sigma)
                    F2 <- lambert2F2elliptic(x, y, tau, nrev, sigma)
                    numberIterations <- numberIterations + 1
                    if(numberIterations > maxIterations) {
                        degree <- degree + 1
                        numberIterations <- 0
                    }
                }
                vc <- sqrt(GM) * (y/sqrt(n) + x/sqrt(m))
                vr <- sqrt(GM) * (y/sqrt(n) - x/sqrt(m))
                ec <- (finalPosition - initialPosition)/diffDistance
                er1 <- initialPosition/sqrt(sum(initialPosition^2))
                er2 <- finalPosition/sqrt(sum(finalPosition^2))
                v1 <- vc*ec + vr*er1
                v2 <- vc*ec - vr*er2
                solutionsList[[2*nrev]] <- list(numberRevs = nrev, path="high", v1=v1, v2=v2)
                # Then the low path
                x <- (xMT + 0.75)/2
                y <- sqrt(1-sigma2*(1-x^2)) * sign(sigma)
                F0 <- lambert2F0elliptic(x, y, nrev, tau)
                F1 <- lambert2F1elliptic(x, y, tau, nrev, sigma)
                F2 <- lambert2F2elliptic(x, y, tau, nrev, sigma)
                convergence <- FALSE
                while(!convergence) {
                    xNew <- x - degree*F0/(F1 + F1/abs(F1)*sqrt((degree-1)^2*F1^2 - degree*(degree-1)*F0*F2))
                    difference <- xNew - x
                    if(abs(difference) < .Machine$double.eps/2) {
                        convergence <- TRUE
                    }
                    x <- xNew
                    y <- sqrt(1-sigma2*(1-x^2)) * sign(sigma)
                    F0 <- lambert2F0elliptic(x, y, nrev, tau)
                    F1 <- lambert2F1elliptic(x, y, tau, nrev, sigma)
                    F2 <- lambert2F2elliptic(x, y, tau, nrev, sigma)
                    numberIterations <- numberIterations + 1
                    if(numberIterations > maxIterations) {
                        degree <- degree + 1
                        numberIterations <- 0
                    }
                }
                vc <- sqrt(GM) * (y/sqrt(n) + x/sqrt(m))
                vr <- sqrt(GM) * (y/sqrt(n) - x/sqrt(m))
                ec <- (finalPosition - initialPosition)/diffDistance
                er1 <- initialPosition/sqrt(sum(initialPosition^2))
                er2 <- finalPosition/sqrt(sum(finalPosition^2))
                v1 <- vc*ec + vr*er1
                v2 <- vc*ec - vr*er2
                solutionsList[[2*nrev+1]] <- list(numberRevs = nrev, path="low", 
                                                  orbitType=orbitType, v1=v1, v2=v2)
            }
        }
    } else if(orbitType == "parabolic") {
        
        x <- 1
        y <- sign(sigma)
        vc <- sqrt(GM) * (y/sqrt(n) + x/sqrt(m))
        vr <- sqrt(GM) * (y/sqrt(n) - x/sqrt(m))
        ec <- (finalPosition - initialPosition)/diffDistance
        er1 <- initialPosition/sqrt(sum(initialPosition^2))
        er2 <- finalPosition/sqrt(sum(finalPosition^2))
        v1 <- vc*ec + vr*er1
        v2 <- vc*ec - vr*er2
        solutionsList <- vector(mode="list", length=1)
        solutionsList[[1]] <- list(numberRevs = 0, path="low", 
                                   orbitType=orbitType, v1=v1, v2=v2)
    } else if (orbitType == "hyperbolic") {
        x <- 10
        y <- sqrt(1-sigma2*(1-x^2))
        if(sigma < 0) {
            y <- -y
        }
        F0 <- lambert2F0hyperbolic(x, y, tau)
        F1 <- lambert2F1hyperbolic(x, y, tau, sigma)
        F2 <- lambert2F2hyperbolic(x, y, tau, sigma)
        signF1 <- sign(F1)
        degree <- 2
        convergence <- FALSE
        numberIterations <- 0
        while(!convergence) {
            xNew <- x - degree*F0/(F1 + F1/abs(F1)*sqrt((degree-1)^2*F1^2 - degree*(degree-1)*F0*F2))
            difference <- xNew - x
            if(abs(difference) < .Machine$double.eps/2) {
                convergence <- TRUE
            }
            x <- xNew
            y <- sqrt(1-sigma2*(1-x^2)) * sign(sigma)
            F0 <- lambert2F0hyperbolic(x, y, tau)
            F1 <- lambert2F1hyperbolic(x, y, tau, sigma)
            F2 <- lambert2F2hyperbolic(x, y, tau, sigma)
            numberIterations <- numberIterations + 1
            if(numberIterations > maxIterations) {
                degree <- degree + 1
                numberIterations <- 0
            }
        }
        vc <- sqrt(GM) * (y/sqrt(n) + x/sqrt(m))
        vr <- sqrt(GM) * (y/sqrt(n) - x/sqrt(m))
        ec <- (finalPosition - initialPosition)/diffDistance
        er1 <- initialPosition/sqrt(sum(initialPosition^2))
        er2 <- finalPosition/sqrt(sum(finalPosition^2))
        v1 <- vc*ec + vr*er1
        v2 <- vc*ec - vr*er2
        solutionsList[[1]] <- list(numberRevs = 0, path=NA, 
                                   orbitType=orbitType, v1=v1, v2=v2)
    }
    return(solutionsList)
}


gooding <- function(initialPosition, finalPosition, initialTime, finalTime,
                    retrogradeTransfer=FALSE, centralBody="Earth", maxIterations=100) {
    checkTimeInput(initialTime, finalTime, "lambert")
    t0 <- 0
    if(is.character(initialTime) & is.character(finalTime)) {
        t <- as.numeric(difftime(initialTime, finalTime, units="secs"))
    } else {
        t <- finalTime
    }
    centralBody <- tolower(centralBody)
    if(!(centralBody %in% c("sun", "mercury", "venus", "earth", "moon", "mars", "jupiter",
                            "saturn", "uranus", "neptune"))) {
        stop("Only the following central bodies are currently supported:
             Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, 
             Neptune and Pluto")
    }
    GM <- switch(centralBody, sun=GM_Sun_DE440, mercury=GM_Mercury_DE440,
                 venus=GM_Venus_DE440, earth=GM_Earth_DE440, moon=GM_Moon_DE440,
                 mars=GM_Mars_DE440, jupiter=GM_Jupiter_DE440, saturn=GM_Saturn_DE440,
                 uranus=GM_Uranus_DE440, neptune=GM_Neptune_DE440, pluto=GM_Pluto_DE440)
}