MJDUTCtoMJDTT <- function(MJD_UTC) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "n")
    leapSeconds <- IERS_results$TAI_UTC
    MJD_TT <- MJD_UTC + (leapSeconds + 32.184)/86400
    return(MJD_TT)
}

MJDTTtoMJDUTC <- function(MJD_TT) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_TT, interp = "n")
    leapSeconds <- IERS_results$TAI_UTC
    MJD_UTC <- MJD_TT - (leapSeconds + 32.184)/86400
    return(MJD_TT)
}

MJDUTCtoMJDUT1 <- function(MJD_UTC) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UTC, interp = "n")
    MJD_UT1 <- MJD_UTC + IERS_results$UT1_UTC/86400
    return(MJD_UT1)
}

MJDUT1toMJDUTC <- function(MJD_UT1) {
    IERS_results <- IERS(asteRiskData::earthPositions, MJD_UT1, interp = "n")
    MJD_UT1 <- MJD_UT1 - IERS_results$UT1_UTC/86400
    return(MJD_UT1)
}
