# TTSecondsToTDBSeconds <- function(TTSeconds) {
#     # We need to do a double nested Newton integration, one for eccentric anomaly
#     # and another one for obtaining TDB seconds from TT
#     current_TDBSeconds <- TTSeconds
#     maxKeplerIterations_outer <- 25
#     keplerAccuracy_outer=10e-9
#     convergence_outer <- FALSE
#     iterations_outer <- 0
#     eccentricity <- 0.0167
#     keplerAccuracy=10e-8
#     maxKeplerIterations=6
#     while(!convergence_outer) {
#         iterations_outer <- iterations_outer + 1
#         meanAnomaly <- TDB_to_TT_M0 + TDB_to_TT_M1*current_TDBSeconds
#         kepler_sol_current <- meanAnomaly
#         convergence <- FALSE
#         iterations <- 0
#         while(!convergence) {
#             iterations <- iterations + 1
#             delta_kepler_sol <- ( (meanAnomaly - kepler_sol_current + eccentricity*sin(kepler_sol_current)) / (1 - eccentricity*cos(kepler_sol_current)) )
#             kepler_sol_current <- kepler_sol_current + delta_kepler_sol
#             if(iterations > maxKeplerIterations | delta_kepler_sol < keplerAccuracy) {
#                 convergence <- TRUE
#             }
#         }
#         eccentricAnomaly <- kepler_sol_current
#         TDB_minus_TT <- TDB_to_TT_K*sin(eccentricAnomaly)
#         current_TDBSeconds <- current_TDBSeconds + TDB_minus_TT
#         if(iterations_outer > maxKeplerIterations_outer | abs(current_TDBSeconds - TTSeconds) < keplerAccuracy_outer) {
#             convergence_outer <- TRUE
#         }
#     }
#     return(current_TDBSeconds)
# }

TTSecondsToTDBSeconds <- function(TTSeconds) {
    meanAnomaly <- TT_to_TDB_M0 + TT_to_TDB_M1*TTSeconds
    kepler_sol_current <- meanAnomaly
    keplerAccuracy=10e-8
    maxKeplerIterations=6
    eccentricity <- TT_to_TDB_EB
    convergence <- FALSE
    iterations <- 0
    while(!convergence) {
        iterations <- iterations + 1
        delta_kepler_sol <- ( (meanAnomaly - kepler_sol_current + eccentricity*sin(kepler_sol_current)) / (1 - eccentricity*cos(kepler_sol_current)) )
        kepler_sol_current <- kepler_sol_current + delta_kepler_sol
        if(iterations > maxKeplerIterations | delta_kepler_sol < keplerAccuracy) {
            convergence <- TRUE
        }
    }
    eccentricAnomaly <- kepler_sol_current
    TDB_minus_TT <- TT_to_TDB_K*sin(eccentricAnomaly)
    return(TTSeconds + TDB_minus_TT)
    return(current_TDBSeconds)
}

TDBSecondsToTTSeconds <- function(TDBSeconds) {
    meanAnomaly <- TDB_to_TT_M0 + TDB_to_TT_M1*TDBSeconds
    kepler_sol_current <- meanAnomaly
    keplerAccuracy=10e-8
    maxKeplerIterations=6
    eccentricity <- TDB_to_TT_EB
    convergence <- FALSE
    iterations <- 0
    while(!convergence) {
        iterations <- iterations + 1
        delta_kepler_sol <- ( (meanAnomaly - kepler_sol_current + eccentricity*sin(kepler_sol_current)) / (1 - eccentricity*cos(kepler_sol_current)) )
        kepler_sol_current <- kepler_sol_current + delta_kepler_sol
        if(iterations > maxKeplerIterations | delta_kepler_sol < keplerAccuracy) {
            convergence <- TRUE
        }
    }
    eccentricAnomaly <- kepler_sol_current
    TDB_minus_TT <- TDB_to_TT_K*sin(eccentricAnomaly)
    return(TDBSeconds - TDB_minus_TT)
}

TTSecondsToTAISeconds <- function(TTSeconds) {
    return(TTSeconds - 32.184)
}

TAISecondsToTTSeconds <- function(TAISeconds) {
    return(TAISeconds + 32.184)
}

TTSecondsToUTCSeconds_J2000 <- function(TTSeconds_J2000) {
    TT_MJD <- TTSeconds_J2000/86400 + MJD_J2000
    IERS_results <- IERS(asteRiskData::earthPositions, TT_MJD, interp = "n")
    leapSeconds <- IERS_results$TAI_UTC
    UTCSeconds_J2000 <- TTSeconds_J2000 - (leapSeconds + 32.184)
    return(UTCSeconds_J2000)
}

UTCSecondsToTTSeconds_J2000 <- function(UTCSeconds_J2000) {
    UTC_MJD <- UTCSeconds_J2000/86400 + MJD_J2000
    IERS_results <- IERS(asteRiskData::earthPositions, UTC_MJD, interp = "n")
    leapSeconds <- IERS_results$TAI_UTC
    TTSeconds_J2000 <- UTCSeconds_J2000 + leapSeconds + 32.184
    return(TTSeconds_J2000)
}

TDBSecondsToUTCSeconds_J2000 <- function(TDBSeconds_J2000) {
    TTSeconds_J2000 <- TDBSecondsToTTSeconds(TDBSeconds_J2000)
    UTCSeconds_J2000 <- TTSecondsToUTCSeconds_J2000(TTSeconds_J2000)
    return(UTCSeconds_J2000)
}

UTCSecondsToTDBSeconds_J2000 <- function(UTCSeconds_J2000) {
    TTSeconds_J2000 <- UTCSecondsToTTSeconds_J2000(UTCSeconds_J2000)
    TDBSeconds_J2000 <- TTSecondsToTDBSeconds(TTSeconds_J2000)
    return(TDBSeconds_J2000)
}


