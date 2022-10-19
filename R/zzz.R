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
                       `asteRiskData` package installed.")
        msg <- paste(msg, collapse="\n")
        stop(msg) 
    }
}

globalVariables(c("gsurf", "re_", "re",
                  "rgas", "dfa", "plg",
                  "ctloc", "stloc", "c2tloc", "s2tloc", "c3tloc", "s3tloc", "apdf",
                  "apt", "end", "pdm", "dm28", "meso_tgn1", "meso_tn1", "pma",
                  "pavgm", "pdl", "ptm", "pt", "ps", "pd", "ptl", "solidEarthTides_dC22dS22",
                  "dm04", "dm16", "dm32", "dm40", "dm01", "dm14", "dd"))
NRLMSISE00.env <- new.env(parent=emptyenv())
NRLMSISE00.env$gsurf <- NULL
NRLMSISE00.env$re_ <- NULL
NRLMSISE00.env$dd <- NULL
NRLMSISE00.env$dm04 <- NULL
NRLMSISE00.env$dm16 <- NULL
NRLMSISE00.env$dm28 <- NULL
NRLMSISE00.env$dm32 <- NULL
NRLMSISE00.env$dm40 <- NULL
NRLMSISE00.env$dm01 <- NULL
NRLMSISE00.env$dm14 <- NULL
NRLMSISE00.env$meso_tn1 <- NULL
NRLMSISE00.env$meso_tn2 <- NULL
NRLMSISE00.env$meso_tn3 <- NULL
NRLMSISE00.env$meso_tgn1 <- NULL
NRLMSISE00.env$meso_tgn2 <- NULL
NRLMSISE00.env$meso_tgn3 <- NULL
NRLMSISE00.env$dfa <- NULL
NRLMSISE00.env$plg <- NULL
NRLMSISE00.env$ctloc <- NULL
NRLMSISE00.env$stloc <- NULL
NRLMSISE00.env$c2tloc <- NULL
NRLMSISE00.env$s2tloc <- NULL
NRLMSISE00.env$s3tloc <- NULL
NRLMSISE00.env$c3tloc <- NULL
NRLMSISE00.env$apdf <- NULL
NRLMSISE00.env$apt <- NULL
