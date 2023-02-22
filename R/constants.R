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
qzms2t <- ((120 - 78)/earthRadius_SGP4)^4

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
