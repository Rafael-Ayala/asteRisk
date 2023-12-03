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
    #return((3*x*lambert2F0elliptic(x, y, N, 0) - 2*(1-sigma^3*x/y))/(1-x^2))
}

lambert2F2elliptic <- function(x, y, N, tau, sigma) {
    return(((1+4*x^2)*lambert2F1elliptic(x, y, N, 0, sigma) + 2*(1-sigma^5*x^3/abs(y)^3))/(x*(1-x^2)))
    #return((3*lambert2F0elliptic(x, y, N, 0) + 5 * x * lambert2F1elliptic(x, y, N, 0, sigma) + 
    #           2*(1-sigma^2)*sigma^3/y^3)/(1-x^2))
}

lambert2F0hyperbolic <- function(x, y, tau) {
    acoth1 <- acoth(x/sqrt(x^2-1))
    acoth2 <- acoth(y/sqrt(y^2-1))
    return((-acoth1 + acoth2 + x*sqrt(x^2 - 1) - y*sqrt(y^2 - 1))/(sqrt((x^2 - 1)^3)) - tau)
}

lambert2F1hyperbolic <- function(x, y, tau, sigma) {
    # return((3*x*lambert2F0hyperbolic(x, y, 0) - 2*(1-sigma^3*x/abs(y)))/(1-x^2))
    return((3*x*lambert2F0hyperbolic(x, y, 0) - 2 * (1-sigma^3*x/abs(y)))/(1-x^2))
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
        t <- as.numeric(difftime(finalTime, initialTime, units="secs"))
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
        F0 <- lambert2F0elliptic(x, y, 0, tau)
        #F0 <- lambert2F0elliptic(x, y, 0, tau) - tau
        F1 <- lambert2F1elliptic(x, y, 0, tau, sigma)
        F2 <- lambert2F2elliptic(x, y, 0, tau, sigma)
        signF1 <- sign(F1)
        degree <- 5
        convergence <- FALSE
        numberIterations <- 0
        while(!convergence) {
            xNew <- x - degree*F0/(F1 + F1/abs(F1)*sqrt((degree-1)^2*F1^2 - degree*(degree-1)*F0*F2))
            difference <- xNew - x
            if(abs(difference) < .Machine$double.eps/2) {
                convergence <- TRUE
            }
            x <- xNew
            if(x <= -1) {
                x <- -1 + .Machine$double.eps
            } else if (x >= 1) {
                x <- 1 - .Machine$double.eps
            }
            y <- sqrt(1-sigma2*(1-x^2)) * sign(sigma)
            F0 <- lambert2F0elliptic(x, y, 0, tau)
            #F0 <- lambert2F0elliptic(x, y, 0, tau) - tau
            F1 <- lambert2F1elliptic(x, y, tau, 0, sigma)
            F2 <- lambert2F2elliptic(x, y, tau, 0, sigma)
            numberIterations <- numberIterations + 1
            if(numberIterations > maxIterations) {
                degree <- degree + 1
                numberIterations <- 0
            }
        }
        degree <- 5
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
                degree <- 1
                while(!convergence) {
                    xMTNew <- xMT - degree*phi0/(phi1 + phi1/abs(phi1)*sqrt((degree-1)^2*phi1^2 - degree*(degree-1)*phi0*phi2))
                    difference <- xMTNew - xMT
                    if(abs(difference) < .Machine$double.eps/2) {
                        convergence <- TRUE
                    }
                    xMT <- xMTNew
                    if(xMT <= -1) {
                        xMT <- -1 + .Machine$double.eps
                    } else if (xMT >= 1) {
                        xMT <- 1 - .Machine$double.eps
                    }
                    yMT <- sqrt(1-sigma2*(1-xMT^2)) * sign(sigma)
                    phi0 <- lambert2Phi0(xMT, yMT, nrev)
                    phi1 <- lambert2Phi1(xMT, yMT, sigma)
                    phi2 <- lambert2Phi2(xMT, yMT, sigma)
                    numberIterations <- numberIterations + 1
                    # if(numberIterations > maxIterations) {
                    #     degree <- degree + 1
                    #     numberIterations <- 0
                    # }
                }
                degree <- 5
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
                degree <- 5
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
                degree <- 5
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
        x <- 1.0001
        y <- sqrt(1-sigma2*(1-x^2))
        if(sigma < 0) {
            y <- -y
        }
        F0 <- lambert2F0hyperbolic(x, y, tau)
        F1 <- lambert2F1hyperbolic(x, y, tau, sigma)
        F2 <- lambert2F2hyperbolic(x, y, tau, sigma)
        signF1 <- sign(F1)
        degree <- 5
        convergence <- FALSE
        numberIterations <- 0
        while(!convergence) {
            xNew <- x - degree*F0/(F1 + F1/abs(F1)*sqrt((degree - 1)^2*F1^2 - degree*(degree-1)*F0*F2))
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
        degree <- 5
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

hybridY <- function(x, sigma) {
    return(sqrt(1 - sigma^2 * (1 - x^2)))
}

hybridPsi <- function(x, y, sigma) {
    if((x >= -1) && (x < 1)) {
        return(acos(x*y + sigma*(1-x^2)))
    } else if (x > 1) {
        return(asinh((y - x*sigma) * sqrt(x^2-1)))
    } else {
        return(0)
    }
}

hybridF0 <- function(x, y, N, tau, sigma) {
    if(N == 0 && x > sqrt(0.6) && x < sqrt(1.4)) {
        eta <- y - sigma*x
        S1 <- (1-sigma-x*eta)/2
        Q <- (4/3) * hyperg_2F1(3, 1, 2.5, S1)
        F0 <- (eta^3*Q + 4*sigma*eta)/2 - tau
    } else {
        psi <- hybridPsi(x, y, sigma)
        F0 <- ( (psi + N*pi) / (sqrt(abs(1 - x^2))) - x + sigma*y ) / (1-x^2) - tau
    }
    return(F0)
}

hybridF1 <- function(x, y, F0, sigma) {
    return((3*F0*x - 2 + 2*sigma^3*(x/y))/(1-x^2))
}

hybridF2 <- function(x, y, F0, F1, sigma) {
    return((3*F0 + 5 * x * F1 + 2*(1-sigma^2)*sigma^3/y^3)/(1-x^2))
}

hybridF3 <- function(x, y, F1, F2, sigma) {
    return((7*x*F2 + 8*F1 - 6*(1 - sigma^2)*sigma^5*x/y^5 )/(1-x^2))
}

hybridMT <- function(sigma, N, maxIterations=100, atol=1e-5, rtol=1e-7) {
    if(sigma==1) {
        xMT <- 0
        yMT <- hybridY(xMT, sigma)
        tauMT <- hybridF0(xMT, yMT, N, 0, sigma)
    } else {
        if(N == 0) {
            xMT <- Inf
            tauMT <- 0
        } else {
            xMT0 <- 0.1
            yMT0 <- hybridY(xMT0, sigma)
            tauMT0 <- hybridF0(xMT0, yMT0, N, 0, sigma)
            xMT <- halleySolver(xMT0, tauMT0, sigma, maxIterations=maxIterations, 
                                atol=atol, rtol=rtol)
            yMT <- hybridY(xMT, sigma)
            tauMT <- hybridF0(xMT, yMT, N, 0, sigma)
        }
    }
    return(list(xMT=xMT, tauMT=tauMT))
}

initialGuessX <- function(tau, sigma, N, lowPath) {
    if(N == 0) {
        tau0 <- acos(sigma) + sigma * sqrt(1 - sigma^2)
        tau1 <- 2 * (1 - sigma^3)/3
        if(tau >= tau0) {
            x0 <- (tau0/tau)^(2/3) - 1
        } else if(tau < tau1) {
            x0 <- 2.5 * tau1/tau * (tau1 - tau)/(1 - sigma^5) + 1
        } else {
            x0 <- (tau0/tau)^(log2(tau1/tau0)) - 1
        }
    } else {
        x01 <- (((N*pi + pi) / (8*tau))^(2/3) - 1) / (((N*pi + pi) / (8*tau)) ^ (2/3) + 1)
        x0r <- (((8*tau) / (N*pi))^(2/3) - 1) / (((8*tau) / (N*pi))^(2/3) + 1)
        if(lowPath) {
            x0 <- max(x01, x0r)
        } else {
            x0 <- min(x01, x0r)
        }
    }
    return(x0)
}

halleySolver <- function(xMT0, tauMT0, sigma, maxIterations=100, atol=1e-5, rtol=1e-7) {
    iteration <- 0
    xMT <- xMT0
    while(iteration < maxIterations) {
        yMT <- hybridY(xMT, sigma)
        F1 <- hybridF1(xMT, yMT, tauMT0, sigma)
        F2 <- hybridF2(xMT, yMT, tauMT0, F1, sigma)
        F3 <- hybridF3(xMT, yMT, F1, F2, sigma)
        xMTNew <- xMT - 2*F1*F2/(2*F2^2-F1*F3)
        if(abs(xMTNew - xMT) < (rtol*abs(xMT) + atol) ) {
            return(xMTNew)
        }
        xMT <- xMTNew
    }
    stop("Failed to converge when solving for tauMT")
}

householderSolver <- function(x0, tau0, sigma, N, maxIterations=100, atol=1e-5, rtol=1e-7) {
    iteration <- 0
    x <- x0
    while(iteration < maxIterations) {
        y <- hybridY(x, sigma)
        F0 <- hybridF0(x, y, N, tau0, sigma)
        F0plustau0 <- F0 + tau0
        F1 <- hybridF1(x, y, F0plustau0, sigma)
        F2 <- hybridF2(x, y, F0plustau0, F1, sigma)
        F3 <- hybridF3(x, y, F1, F2, sigma)
        xNew <- x - F0 * ( (F1^2 - F0*F2/2) / (F1 * (F1^2 - F0 * F2) + F3 * F1^2/6) )
        if(abs(xNew - x) < (rtol*abs(x) + atol) ) {
            return(xNew)
        }
        iteration <- iteration + 1
        x <- xNew
    }
    stop("Failed to converge when solving for x")
}

hybridReconstruct <- function(x, y, r1, r2, sigma, gamma, rho, zeta) {
    Vr1 <- gamma*((sigma*y - x) - rho * (sigma * y +x))/r1
    Vr2 <- -gamma*((sigma*y - x) + rho * (sigma * y +x))/r2
    Vt1 <- gamma * zeta * (y + sigma * x)/r1
    Vt2 <- gamma * zeta * (y + sigma * x)/r2
    return(list(
        Vr1=Vr1,
        Vr2=Vr2,
        Vt1=Vt1,
        Vt2=Vt2
    ))
}

lambert <- function(initialPosition, finalPosition, initialTime, finalTime,
                    retrogradeTransfer=FALSE, centralBody="Earth", maxIterations=2000,
                    atol=0.00001, rtol=0.000001) {
    checkTimeInput(initialTime, finalTime, "lambert")
    t0 <- 0
    if(is.character(initialTime) & is.character(finalTime)) {
        t <- as.numeric(difftime(finalTime, initialTime, units="secs"))
    } else {
        t <- finalTime
    }
    if(is.numeric(centralBody)) {
        GM <- centralBody
    } else if(is.character(centralBody)){
        centralBody <- tolower(centralBody)
        if(!(centralBody %in% c("sun", "mercury", "venus", "earth", "moon", "mars", 
                                "jupiter", "saturn", "uranus", "neptune"))) {
            stop("Only the following central bodies are currently supported:
             Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, 
             Neptune and Pluto")
        }
    }

    GM <- switch(centralBody, sun=GM_Sun_DE440, mercury=GM_Mercury_DE440,
                 venus=GM_Venus_DE440, earth=GM_Earth_DE440, moon=GM_Moon_DE440,
                 mars=GM_Mars_DE440, jupiter=GM_Jupiter_DE440, saturn=GM_Saturn_DE440,
                 uranus=GM_Uranus_DE440, neptune=GM_Neptune_DE440, pluto=GM_Pluto_DE440)
    initialDistance <- sqrt(sum(initialPosition^2))
    finalDistance <- sqrt(sum(finalPosition^2))
    diffVector <- finalPosition - initialPosition
    diffDistance <- sqrt(sum(diffVector^2))
    s <- (initialDistance + finalDistance + diffDistance)/2
    ir1 <- initialPosition/initialDistance
    ir2 <- finalPosition/finalDistance
    ih <- vectorCrossProduct3D(ir1, ir2)
    ih <- ih/sqrt(sum(ih^2))
    sigma <- sqrt(1 - min(1, diffDistance/s))
    if(ih[3] < 0) {
        sigma <- -sigma
        it1 <- vectorCrossProduct3D(ir1, ih)
        it2 <- vectorCrossProduct3D(ir2, ih)
    } else {
        it1 <- vectorCrossProduct3D(ih, ir1)
        it2 <- vectorCrossProduct3D(ih, ir2)
    }
    if(retrogradeTransfer) {
        sigma <- -sigma
        it1 <- -it1
        it2 <- -it2
    }
    tau <- sqrt(2*GM/s^3)*t
    tauParabolic <- 2/3*(1-sigma^3)
    diffTauParabolic <- tau - tauParabolic
    if(abs(diffTauParabolic) <= .Machine$double.eps) {
        orbitType <- "parabolic"
    } else if(diffTauParabolic > 0) {
        orbitType <- "elliptic"
    } else if(diffTauParabolic < 0) {
        orbitType <- "hyperbolic"
    }
    if(orbitType == "elliptic") {
        Nmax <- floor(tau/pi)
    } else {
        Nmax <- 0
    }
    solutionsList <- vector(mode="list", length=1 + Nmax*2)
    tauME <- acos(sigma) + sigma*sqrt(1 - sigma^2)
    if(abs(tauME - tau) <= .Machine$double.eps) {
        path <- "Minimum Energy"
        lowPath <- FALSE
    } else if(tau < tauME) {
        path <- "Low"
        lowPath <- TRUE
    } else {
        path <- "High"
        lowPath <- FALSE
    }
    # First for N=0/ or non-elliptical cases
    x0 <- initialGuessX(tau, sigma, 0, lowPath)
    y0 <- hybridY(x0, sigma)
    x <- householderSolver(x0, tau, sigma, 0, maxIterations = maxIterations,
                           atol = atol, rtol=rtol)
    y <- hybridY(x, sigma)
    gamma <- sqrt(GM*s/2)
    rho <- (initialDistance - finalDistance)/diffDistance
    # Zeta originally named sigma in Izzy's
    zeta <- sqrt(1-rho^2)
    auxV <- hybridReconstruct(x, y, initialDistance, finalDistance,
                              sigma, gamma, rho, zeta)
    v1 <- auxV$Vr1*(initialPosition/initialDistance) + auxV$Vt1*it1
    v2 <- auxV$Vr2*(finalPosition/finalDistance) + auxV$Vt2*it2
    solutionsList[[1]] <- list(numberRevs = 0, path=path, 
                               orbitType=orbitType, v1=v1, v2=v2)
    if(Nmax >= 1) {
        for(nrev in 1:Nmax) {
            # tauME <- nrev*pi + acos(sigma) + sigma*sqrt(1-sigma^2)
            # first get tauMT
            minimumT <- hybridMT(sigma, nrev)
            xMT <- minimumT$xMT
            yMT <- hybridY(xMT, sigma)
            tauMT <- minimumT$tauMT
            # First the high path
            lowPath <- FALSE
            x0 <- initialGuessX(tau, sigma, nrev, lowPath)
            y0 <- hybridY(x0, sigma)
            x <- householderSolver(x0, tau, sigma, nrev, maxIterations = maxIterations,
                                   atol = atol, rtol=rtol)
            y <- hybridY(x, sigma)
            auxV <- hybridReconstruct(x, y, initialDistance, finalDistance,
                                      sigma, gamma, rho, zeta)
            v1 <- auxV$Vr1*(initialPosition/initialDistance) + auxV$Vt1*it1
            v2 <- auxV$Vr2*(finalPosition/finalDistance) + auxV$Vt2*it2
            solutionsList[[2*nrev]] <- list(numberRevs = nrev, path="high", 
                                            orbitType=orbitType, v1=v1, v2=v2)
            # then for the low path
            lowPath <- TRUE
            x0 <- initialGuessX(tau, sigma, nrev, lowPath)
            y0 <- hybridY(x0, sigma)
            x <- householderSolver(x0, tau, sigma, nrev, maxIterations = maxIterations,
                                   atol = atol, rtol=rtol)
            y <- hybridY(x, sigma)
            auxV <- hybridReconstruct(x, y, initialDistance, finalDistance,
                                      sigma, gamma, rho, zeta)
            v1 <- auxV$Vr1*(initialPosition/initialDistance) + auxV$Vt1*it1
            v2 <- auxV$Vr2*(finalPosition/finalDistance) + auxV$Vt2*it2
            solutionsList[[2*nrev + 1]] <- list(numberRevs = nrev, path="low",
                                                orbitType=orbitType, v1=v1, v2=v2)
        }
    }
    return(solutionsList)
}
