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

# Legendre polynomials
n <- 70
m <- 70
