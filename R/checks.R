checkTimeInput <- function(initialTime, targetTime, algorithm) {
    if(algorithm == "sgp4") {
        if(!(is.numeric(targetTime) | is.character(targetTime)) | is.nan(targetTime)) {
            stop("Please provide a target time as an amount of minutes or as a
                 date-time string in UTC")
        }
        if(!is.numeric(targetTime)) {
            attempt_target_conversion <- tryCatch({
                as.POSIXct(targetTime, tz="UTC")
            }, error = function(e) return(NA))
            if(is.na(attempt_target_conversion)) {
                stop("Please provide a target time as an amount of minutes or as a
                     date-time string in UTC")
            }
            attempt_initial_conversion <- tryCatch({
                as.POSIXct(initialTime, tz="UTC")
            }, error = function(e) return(NA))
            if(is.na(attempt_initial_conversion)) {
                stop("Target time provided as a date-time string. A valid date-time
                      string must also be provided as initialTime")
            }
        }
    }
    if(algorithm == "sdp4") {
        attempt_initial_conversion <- tryCatch({
            as.POSIXct(initialTime, tz="UTC")
        }, error = function(e) return(NA))
        if(is.na(attempt_initial_conversion)) {
            stop("Please provide a valid date-time string in UTC specifying
                 the epoch time")
        }
        if(!(is.numeric(targetTime) | is.character(targetTime)) | is.nan(targetTime)) {
            stop("Please provide a target time as an amount of minutes or as a
                 date-time string in UTC")
        }
        if(!is.numeric(targetTime)) {
            attempt_target_conversion <- tryCatch({
                as.POSIXct(targetTime, tz="UTC")
            }, error = function(e) return(NA))
            if(is.na(attempt_target_conversion)) {
                stop("Please provide a target time as an amount of minutes or as a
                 date-time string in UTC")
            }
        }
    }
    if(algorithm == "lambert") {
        if(!(is.numeric(targetTime) | is.character(targetTime)) | is.nan(targetTime)) {
            stop("Please provide a target time as an amount of seconds or as a
                 date-time string in UTC")
        }
        if(!is.numeric(targetTime)) {
            attempt_target_conversion <- tryCatch({
                as.POSIXct(targetTime, tz="UTC")
            }, error = function(e) return(NA))
            if(is.na(attempt_target_conversion)) {
                stop("Please provide a target time as an amount of seconds or as a
                     date-time string in UTC")
            }
            attempt_initial_conversion <- tryCatch({
                as.POSIXct(initialTime, tz="UTC")
            }, error = function(e) return(NA))
            if(is.na(attempt_initial_conversion)) {
                stop("Target time provided as a date-time string. A valid date-time
                      string must also be provided as initialTime")
            }
        }
    }
}

checkTLElist <- function(readTLE) {
    if(!all(c("NORADcatalogNumber", "classificationLevel", "internationalDesignator",
              "launchYear", "launchNumber", "launchPiece", "dateTime", "elementNumber",
              "inclination", "ascension", "eccentricity", "perigeeArgument", "meanAnomaly",
              "meanMotion", "meanMotionDerivative", "meanMotionSecondDerivative",
              "Bstar", "ephemerisType", "epochRevolutionNumber", "objectName") %in% names(readTLE))) {
        stop("Please provide a valid TLE read by either the readTLE or parseTLElines functions")
    }
}
