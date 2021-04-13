# Constants and tables

Cnm <- matrix(0L, nrow=181, ncol=181)
Snm <- matrix(0L, nrow=181, ncol=181)

file_name <- "./R/hpop_files/GGM03S.txt"
lines_Cnm_matrix <- readLines(file_name)

line <- 1
for(n in 1:181) {
    for (m in 1:n) {
        line_content <- strsplit(trimws(lines_Cnm_matrix[line]), "\\s+")[[1]]
        Cnm[n, m] <- as.numeric(line_content[3])
        Snm[n, m] <- as.numeric(line_content[4])
        line <- line + 1
    }
} # verified Cnm[103, 30] = 7.629580795415e-10 , Snm[103, 30] = 1.68949165464e-09

earthPositions <- read.table("R/hpop_files/EOP-All.txt")
colnames(earthPositions) <- c("Date(year)-0h UTC",
                              "Date(month)-0h UTC",
                              "Date(day)-0h UTC",
                              "Modified Julian Date",
                              "x", "y",
                              "UT1 - UTC",
                              "LOD", "dPsi", "dEpsilon",
                              "dX", "dY", "DAT (TAI-UTC)")

spaceWeather <- read.table("R/hpop_files/SW-All.txt")
colnames(spaceWeather) <- c("Date(year)",
                            "Date(month)",
                            "Date(day)",
                            "BSRN",
                            "ND", paste(rep("Kp", 8), 1:8, sep=""),
                            "Sum", paste(rep("Ap", 8), 1:8, sep=""),
                            "Avg", "Cp", "C9", "ISN", "AdjF10.7", "Q",
                            "AdjCtr81", "AdjLst81", "ObsF10.7", "ObsCtr81",
                            "ObsLst81")

solarStorms <- read.table("R/hpop_files/SOLFSMY.TXT")
colnames(solarStorms) <- c("Year", "Day", "JulianDay", "F10", "F54b",
                           "S10", "S54b", "M10", "M54b", "Y10", "Y54b", "Ssrc")

geomagneticStormsDTC <- read.table("R/hpop_files/DTCFILE.TXT")
colnames(geomagneticStormsDTC) <- c("DTC", "Year", "Day", 
                                    paste(rep("DTC", 24), 1:24, sep=""))

geomagneticStormsAP <- read.table("R/hpop_files/SOLRESAP.TXT")  # different columns from Matlab implementation
colnames(geomagneticStormsAP) <- c("Day", "F10", "F10B", 
                                   paste(rep("Ap", 8), 1:8, sep=""),
                                   "Year", "F")

# Initial values of auxiliary parameters

Mjd_UTC <- 0
area_solar <- 0
area_drag <- 0
mass <- 0
Cr <- 0
Cd <- 0
n <- 0
m <- 0
sun <- 0
moon <- 0
sRad <- 0
drag <- 0
planets <- 0
SolidEarthTides <- 0
OceanTides <- 0
Relativity <- 0

# Mathematical constants
const_pi2 <- 2*pi                # 2pi
const_Rad  <- pi/180              # Radians per degree
const_Deg <- 180/pi              # Degrees per radian
const_Arcs <- 3600*180/pi         # Arcseconds per radian
DAS2R <- 4.848136811095359935899141e-6 # Arcseconds to radians 
TURNAS <- 1296000.0 # Arcseconds in a full circle 

# Date-related constants
DJ00 <- 2451545.0
DJC <- 36525.0
MJD_J2000 <- 51544.5 # MJD for J2000.0

# Acceleration equations constants
DE436coeffs <- read.csv("R/hpop_files/DE436Coeff.csv", header=FALSE, sep = ",", dec = ".")
# Table 6.5a IERS 2010
solidEarthTides_dC21dS21 <- read.csv("R/hpop_files/solidEarthTides_dC21dS21.csv", header=FALSE) 
# Table 6.5b IERS 2010
solidEarthTides_dC22dS22 <- read.csv("R/hpop_files/solidEarthTides_dC22dS22.csv", header=FALSE)
# Table 6.5c IERS 2010
solidEarthTides_dC20 <- read.csv("R/hpop_files/solidEarthTides_dC20.csv", header=FALSE)
gm <- 398600.4415e9 # GGM03S
r_ref <- 6378.1363e3 # radius of Earth in m
R_sun <- 696000e3 # radius of Sun in m
R_moon <- 1738e3 # radius of Moon in m
AU <- 149597870699.999988 # meters in an AU
c_light <- 299792457.999999984 # speed of light in m/s
solarPressureConst <- 1367/c_light # solar radiation pressure at 1 AU in N/m^2 = 1367 W/m^2)
omegaEarth <- 15.04106717866910/3600*(pi/180) # Earth rotation (derivative of GSMT at J2000) in rad/s

# Gravitational coefficients in m3/s2, all in DE436 except Earth in WGS84 system
GM_Earth <- 398600.4418e9
GM_Sun <- 132712440041.9394e9
GM_Moon <- GM_Earth / 81.3005682168675747
GM_Mercury <- 22031.78000000002e9
GM_Venus <- 324858.5920000001e9
GM_Mars <- 42828.37521400003e9
GM_Jupiter <- 126712764.1334462e9
GM_Saturn <- 37940585.20000001e9
GM_Uranus <- 5794556.465751793e9
GM_Neptune <- 6836527.100580024e9
GM_Pluto <- 975.5011758767654e9	  			  

# Ephemeris model constants
EMRAT = 81.3005682168675747
EMRAT1 = 1/(1+EMRAT)

# CHANGE THIS TO USER DEFINABLE VALUES
n <- 70
m <- 70

