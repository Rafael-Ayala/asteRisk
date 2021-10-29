## SGP4 constants

# XKMPER <- 6378.135 # Earth radius in kilometers
earthRadius_SGP4 <- 6378.135 # Earth radius in kilometers, value for SGP/SDP. From WGS-72
ae <- 1 # distance units per earth radii
# earthEccentricity <- 0.0818191908 # Earth eccentricity

J2 <- 1.082616e-3 # second gravitational zonal harmonic of Earth
J3 <- -2.53881e-6 # third gravitational zonal harmonic of Earth
J4 <- -1.65597e-6 # fourth gravitational zonal harmonic of Earth
k2 <- 0.5 * J2 * ae ^ 2
k4 <- (-3 / 8) * J4 * ae ^ 4
# ke <- 7.43669161e-2 OLD IMPLEMENTATION OF WGS72
GM_Earth_WGS72 <- 398600.8
ke <- 60/(sqrt(earthRadius_SGP4^3/GM_Earth_WGS72))
A30 <- -J3 * ae ^ 3

q0 <- 120 / 6378.135 + 1 # q0 parameter for SGP4/SGP8 density function
s <- 78 / 6378.135 + 1 # s parameter for SGP4/SGP8 density function

## Constants required for SDP4

STEP <- 720.0
# MAX_INTEGRATE <- STEP * 10000
ZNS <- 1.19459E-5
C1SS <- 2.9864797e-6
ZES <- 0.01675
ZNL <- 1.5835218e-4
ZEL <- 0.05490
C1L <- 4.7968065e-7
ZSINIS <- 0.39785416
ZCOSIS <- 0.91744867
ZCOSGS <- 0.1945905
ZSINGS <- -0.98088458
THDT <- 4.37526908801129966e-3
ROOT22 <- 1.7891679e-6
ROOT32 <- 3.7393792e-7
ROOT44 <- 7.3636953e-9
ROOT52 <- 1.1428639e-7
ROOT54 <- 2.1765803e-9
G22 <- 5.7686396
G32 <- 0.95240898
G44 <- 1.8014998
G52 <- 1.0508330
G54 <- 4.4108898
Q22 <- 1.7891679e-6
Q31 <- 2.1460748e-6
Q33 <- 2.2123015e-7

## Constants required for calculation of acceleration and conversion of
## coordinates systems

# WGS84 constants
earthRadius_WGS84 <- 6.378137e6 # WGS84 Earth radius in meters
earthEccentricity_WGS84 <- 8.18191908426214947083e-2
earthFlatteningFactor_WGS84 <- 1/298.257223563
J2_WGS84 <- 0.00108262998905
# WGS84_E2D2 <-earthEccentricity_WGS84^2/2
# WGS84_E4D4 <- earthEccentricity_WGS84^4/4
# WGS84_INVA2 <- 1/earthRadius_WGS84^2
# WGS84_P1ME2 <- 1 - earthEccentricity_WGS84^2
# WGS84_P1ME2DA2 <- (1 - earthEccentricity_WGS84^2) / (earthRadius_WGS84^2)
# WGS84_HMIN <- 2.25010182030430273673e-14

# Mathematical constants
const_Arcs <- 3600*180/pi         # Arcseconds per radian
DAS2R <- 4.848136811095359935899141e-6 # Arcseconds to radians 
TURNAS <- 1296000.0 # Arcseconds in a full circle 
inv_cbr_2 <- 1/(2^(1/3)) # Inverse of cubic root of 2

# Date-related constants
JD_J2000_0 <- 2451545.0 # Julian Day of J2000.0 epoch.
DJC <- 36525.0
MJD_J2000 <- 51544.5 # MJD for J2000.0

# Gravitational coefficients in m3/s2, all in DE436 except Earth in WGS84 system
GM_Earth_TCB <- 398600.4418e9 # WGS84 system (TCB compatible)
GM_Earth_TT <- 398600.4415e9 # GGM03S - TT compatible
GM_Sun <- 132712440041.9394e9
GM_Moon <- GM_Earth_TCB / 81.3005682168675747
GM_Mercury <- 22031.78000000002e9
GM_Venus <- 324858.5920000001e9
GM_Mars <- 42828.37521400003e9
GM_Jupiter <- 126712764.1334462e9
GM_Saturn <- 37940585.20000001e9
GM_Uranus <- 5794556.465751793e9
GM_Neptune <- 6836527.100580024e9
GM_Pluto <- 975.5011758767654e9

# Other constants required for calculation of acceleration
earthRadius_EGM96 <- 6378.1363e3 # radius of Earth in m, EGM96 model
sunRadius <- 696000e3 # radius of Sun in m
moonRadius <- 1738e3 # radius of Moon in m
AU <- 149597870699.999988 # meters in one AU
c_light <- 299792457.999999984 # speed of light in m/s
solarPressureConst <- 1367/c_light # solar radiation pressure at 1 AU in N/m^2 = 1367 W/m^2)
omegaEarth <- 15.04106717866910/3600*(pi/180) # Earth rotation (derivative of GSMT at J2000) in rad/s

# Ephemeris model constants
EMRAT = 81.3005682168675747
EMRAT1 = 1/(1+EMRAT)
