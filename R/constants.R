XKMPER <- 6378.135 # Earth radius in kilometers
ae <- 1 # distance units per earth radii
earthEccentricity <- 0.0818191908 # Earth eccentricity

J2 <- 1.082616e-3 # second gravitational zonal harmonic of Earth
J3 <- -2.53881e-6 # third gravitational zonal harmonic of Earth
J4 <- -1.65597e-6 # fourth gravitational zonal harmonic of Earth
k2 <- 0.5 * J2 * ae ^ 2
k4 <- (-3 / 8) * J4 * ae ^ 4
ke <- 7.43669161e-2
A30 <- -J3 * ae ^ 3

q0 <- 120 / 6378.135 + 1 # parameter for SGP4/SGP8 density function
s <- 78 / 6378.135 + 1 # another parameter for SGP4/SGP8 density function

# Julian Day of J2000.0 epoch.
JD_J2000_0 <- 2451545.0

## Constants required for SDP4

STEP <- 720.0
MAX_INTEGRATE <- STEP * 10000
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


## Initial values for SDP4

atime <- 0.0
xli <- 0.0
xni <- 0.0
xfact <- 0.0
ssl <- 0.0
ssg <- 0.0
ssh <- 0.0
sse <- 0.0
ssi <- 0.0
xlamo <- 0.0
gmst <- 0.0
del1 <- 0.0
del2 <- 0.0
del3 <- 0.0
fasx2 <- 0.0
fasx4 <- 0.0
fasx6 <- 0.0
d2201 <- 0.0
d2211 <- 0.0
d3210 <- 0.0
d3222 <- 0.0
d4410 <- 0.0
d4422 <- 0.0
d5220 <- 0.0
d5232 <- 0.0
d5421 <- 0.0
d5433 <- 0.0
xnddt <- 0.0
xndot <- 0.0
xldot <- 0.0
zmos <- 0.0
se2 <- 0.0
se3 <- 0.0
si2 <- 0.0
si3 <- 0.0
sl2 <- 0.0
sl3 <- 0.0
sl4 <- 0.0
sgh2 <- 0.0
sgh3 <- 0.0
sgh4 <- 0.0
sh2 <- 0.0
sh3 <- 0.0
zmol <- 0.0
ee2 <- 0.0
e3 <- 0.0
xi2 <- 0.0
xi3 <- 0.0
xl2 <- 0.0
xl3 <- 0.0
xl4 <- 0.0
xgh2 <- 0.0
xgh3 <- 0.0
xgh4 <- 0.0
xh2 <- 0.0
xh3 <- 0.0
pe <- 0.0
pinc <- 0.0
pgh <- 0.0
ph <- 0.0
pl <- 0.0
pgh0 <- 0.0
ph0 <- 0.0
pe0 <- 0.0
pinc0 <- 0.0
pl0 <- 0.0
se <- 0
si <- 0
sl <- 0
sgh <- 0
shdq <- 0

# WGS84 constants

WGS84_A <- 6.378137e6
WGS84_E <- 8.18191908426214947083e-2
WGS84_E2D2 <-WGS84_E^2/2
WGS84_E4D4 <- WGS84_E^4/4
WGS84_INVA2 <- 1/WGS84_A^2
WGS84_P1ME2 <- 1 - WGS84_E^2
WGS84_P1ME2DA2 <- (1 - WGS84_E^2) / (WGS84_A^2)
WGS84_HMIN <- 2.25010182030430273673e-14

# Other constants

inv_cbr_2 <- 1/(2^(1/3))
earth_mu <- 3.986004418e14
