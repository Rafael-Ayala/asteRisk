## SGP4 constants

# XKMPER <- 6378.135 # Earth radius in kilometers
earthRadius_SGP4 <- 6378.135 # Earth radius in kilometers, value for SGP/SDP. From WGS-72
ae_SGP4 <- 1 # distance units per earth radii
# earthEccentricity <- 0.0818191908 # Earth eccentricity

J2_SGP4 <- 1.082616e-3 # second gravitational zonal harmonic of Earth
J3_SGP4 <- -2.53881e-6 # third gravitational zonal harmonic of Earth
J4_SGP4 <- -1.65597e-6 # fourth gravitational zonal harmonic of Earth
k2_SGP4 <- 0.5 * J2_SGP4 * ae_SGP4 ^ 2
k4_SGP4 <- (-3 / 8) * J4_SGP4 * ae_SGP4 ^ 4
# ke <- 7.43669161e-2 OLD IMPLEMENTATION OF WGS72
GM_Earth_WGS72_SGP4 <- 398600.8
ke_SGP4 <- 60/(sqrt(earthRadius_SGP4^3/GM_Earth_WGS72_SGP4))
A30_SGP4 <- -J3_SGP4 * ae_SGP4 ^ 3

q0_SGP4 <- 120 / 6378.135 + 1 # q0 parameter for SGP4/SGP8 density function
s0_SGP4 <- 78 / 6378.135 + 1 # s (s0) parameter for SGP4/SGP8 density function
qzms2t_SGP4 <- ((120 - 78)/earthRadius_SGP4)^4

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

# Ephemeris model constants
EMRAT_DE440 = 81.3005682214972154
EMRAT1_DE440 = 1/(1+EMRAT_DE440)

# Gravitational coefficients in m3/s2, all from DE440 except Earth also in WGS84 system and GGM05C
GM_Earth_TDB <- 398600.4418e9 # WGS84 system (TDB compatible)
GM_Earth_TT <- 398600.4415e9 # GGM03S and GGM05C - TT compatible
GM_Earth_DE440 <- 398600.435507e9
GM_Sun_DE440 <- 132712440041.279419e9
GM_Moon_DE440 <- 4902.800118e9
GM_Mercury_DE440 <- 22031.78000000002e9
GM_Venus_DE440 <- 324858.5920000001e9
GM_Mars_DE440 <- 42828.375816e9
GM_Jupiter_DE440 <- 126712764.1e9
GM_Saturn_DE440 <- 37940584.841800e9
GM_Uranus_DE440 <- 5794556.4e9
GM_Neptune_DE440 <- 6836527.100580e9
GM_Pluto_DE440 <- 975.500000e9
GM_Moon_GRGM1200B <- 4902.8001224453001e9

# Other constants required for calculation of acceleration
earthRadius_EGM96 <- 6378.1363e3 # radius of Earth in m, EGM96 model
sunRadius <- 696000e3 # radius of Sun in m
moonRadius <- 1738e3 # radius of Moon in m
AU <- 149597870699.999988 # meters in one AU
c_light <- 299792457.999999984 # speed of light in m/s
solarPressureConst <- 1367/c_light # solar radiation pressure at 1 AU in N/m^2 = 1367 W/m^2)
omegaEarth <- 15.04106717866910/3600*(pi/180) # Earth rotation (derivative of GSMT at J2000) in rad/s
GMratioMoonEarth <- GM_Moon_DE440/GM_Earth_TT
GMratioSunEarth <- GM_Sun_DE440/GM_Earth_TT
GMratioEarthMoon <- GM_Earth_TT/GM_Moon_GRGM1200B
GMratioSunMoon <- GM_Sun_DE440/GM_Moon_GRGM1200B
# Elastic Earth solid tide Love numbers (n=2,3, m from 0 to n. n in rows, m in cols)
# note 0 in first row to make matrix square
knm0 <- matrix(c(0.29525, 0.29470, 0.29801, 0,
                 0.093, 0.093, 0.093, 0.094), 
               nrow=2, ncol=4, byrow=TRUE)
# k2m + values, used for corrections of C4m and S4m, m from 0 to 2
knmplus <- c(-0.00087, -0.00079, -0.00057) 
# Anelastic Earth solid tide Love numbers. Same comments as before apply
# In short, imaginary components are 0 in elastic Earth model, but not here when
# n=2 (but yes when n=3)
reKnm0An <- matrix(c(0.30190, 0.29830, 0.30102, 0,
                       0.093, 0.093, 0.093, 0.094), 
                     nrow=2, ncol=4, byrow=TRUE)
imKnm0An <- matrix(c(0, -0.00144, -0.00130, 0,
                     0, 0, 0, 0), 
                   nrow=2, ncol=4, byrow=TRUE)
knmplusAn <- c(-0.00089, -0.00080, -0.00057) 

# Radius of sphere of influence, from Table 4.1 Lunar Transfer Orbits by W. Seefelder
SOIradii <- c(0.112, 0.616, 0.929, 0.578, 48.2, 54.5, 51.9, 86.8, 34.1)*1e9
names(SOIradii) <- c("Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto")
SOIMoon <- 0.0661e9 # in relation to the Earth

# Vector of numbers with names set to the possible central bodies. This is necessary
# because the solver does not seem to allow to return values other than numbers

centralBodiesNum <- 1:11
names(centralBodiesNum) <- c("SunSSBarycentric", "Mercury", "Venus", "Earth", "Moon",
                             "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto")

# Constants for time conversions

# Strange: NAIF document says amplitude (K) is 1.658e-3, but their kernels use 
# 1.657e-3. ESA's documentation uses 1.658e-3.
# But there is additional small differences in the values of M0 and M1
# after converting ESA's values to radians. Additionally, ESA differentiates
# TDB and TCB
TDB_to_TT_K <- 1.657e-3

# ESA's docs use 0.0167 for next one
TDB_to_TT_EB <- 1.671e-2

TDB_to_TT_M0 <- 6.239996
TDB_to_TT_M1 <- 1.99096871e-7

# For reverse transformation: I take values from ESA's navipedia
# since that expression is dependent on TDT MJD, so I guess it works
# in that direction. Maybe that explains small difference between values?

TT_to_TDB_K <- 1.658e-3
TT_to_TDB_EB <- 1.67e-2
TT_to_TDB_M0 <- 6.2400407681
TT_to_TDB_M1 <- 1.9909687368e-07

# Reference to convert to UTC seconds since J2000
J2000_POSIXct <- as.POSIXct("2000-01-01 12:00:00", tz="UTC")

# Constants for anelastic solid Moon tide corrections according to 2015 paper,
# using Debye model for k2 and k2/Q

tau2 <- 15.7 # in days
Pref <- 27.212
vref <- (2*pi)/Pref
k2ByQref <- 6.4e-4
k2ref <- 0.02416

# constants for elastic solid Moon tide corrections
k20moon <- 0.02408
k21moon <- 0.02414
k22moon <- 0.02394
solidMoonTidesSimple <- matrix(c(1, 0, 0, 0, 0, 27.555, -1.89e-10, 0, 0, 3.28e-10, -4.41e-10,
                                 0, 0, 1, 0, 0, 27.212, 0, 4.72e-10, -0.02e-10, 0, 0,
                                 -1, 0, 0, 2, 0, 31.812, -0.36e-10, 0, 0, 0.58e-10, -0.88e-10,
                                 0, 0, 0, 2, 0, 14.765, -0.32e-10, 0, 0, 0.75e-10, -0.48e-10,
                                 2, 0, 0, 0, 0, 13.777, -0.16e-10, 0, 0, 0.51e-10, -0.51e-10,
                                 1, 0, 1, 0, 0, 13.691, 0, 0.64e-10, 0.26e-10, 0, 0), 
                               byrow=TRUE, nrow=6)
rownames(solidMoonTidesSimple) <- c("l", "F", "2D-l", "2D", "2l", "F+l")
colnames(solidMoonTidesSimple) <- c("l", "lprime", "F", "D", "Omega", "Period(days)", 
                                    "C20qNorm", "C21qNorm", "S21qnorm", "C22qNorm",
                                    "S22qNorm")

## File downloads URLs etc

JPLbinPlanetEphemeridesURL <- "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/"
asteRiskDataFolder <- R_user_dir("asteRisk")
