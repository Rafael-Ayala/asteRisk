.pkgenv <- new.env(parent=emptyenv())

.onLoad  <- function(libname, pkgname) {
    has_data <- requireNamespace("asteRiskData", quietly = TRUE)
    .pkgenv[["has_data"]] <- has_data
}

.onAttach <- function(libname, pkgname) {
    if (!.pkgenv$has_data) {
        msg <- strwrap("To use all functionalities available in this package, 
                       including the high-precision numerical orbit propagator,
                       you must install the asteRiskData package. To install it,
                       run `install.packages('asteRiskData', repos='https://rafael-ayala.github.io/drat/')`.")
        msg <- paste(msg, collapse="\n")
        packageStartupMessage(msg)
    }
}

hasData <- function(has_data = .pkgenv$has_data) {                                
    if (!has_data) {
        msg <- strwrap("To use this function, you must have the
                       `asteriskData` package installed.")
        msg <- paste(msg, collapse="\n")
        stop(msg) 
    }
}

globalVariables(c("earthPositions", "DE436coeffs", "gsurf", "re_", "re", "spaceWeather",
                  "rgas", "Cnm", "Snm", "solidEarthTides_dC21dS21", "dfa", "plg",
                  "ctloc", "stloc", "c2tloc", "s2tloc", "c3tloc", "s3tloc", "apdf",
                  "apt", "end", "pdm", "dm28", "meso_tgn1", "meso_tn1", "pma",
                  "pavgm", "pdl", "ptm", "pt", "ps", "pd", "ptl", "solidEarthTides_dC22dS22",
                  "dm04", "dm16", "dm32", "dm40", "dm01", "dm14", "dd", "xls", 
                  "NLS", "xpl", "NPL", "s0", "s1", "s2", "s3", "s4", "w5"))
