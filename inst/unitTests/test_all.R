library(RUnit)

## Test parseTLElines

testTLElines <- c("ITALSAT 2", "1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600", "2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119")

checkEquals(parseTLElines(testTLElines)$epochRevolutionNumber, 3611)

## Test readTLE

testTLEs <- readTLE(paste0(path.package("asteRisk"), "/testTLE.txt"))
checkTrue(length(testTLEs) == 29)

## Test sgp4, sdp4 and sgpdp4

testSGDP4_1 <- sgdp4(
    n0 = testTLEs[[29]]$meanMotion * ((2 * pi) / (1440)),
    e0 = testTLEs[[29]]$eccentricity,
    i0 = asteRisk:::deg2rad(testTLEs[[29]]$inclination),
    M0 = asteRisk:::deg2rad(testTLEs[[29]]$meanAnomaly),
    omega0 = asteRisk:::deg2rad(testTLEs[[29]]$perigeeArgument),
    OMEGA0 = asteRisk:::deg2rad(testTLEs[[29]]$ascension),
    Bstar = testTLEs[[29]]$Bstar,
    initialDateTime = testTLEs[[29]]$dateTime,
    targetTime = 80
)

testSGDP4_2 <- sgdp4(
    n0 = testTLEs[[17]]$meanMotion * ((2 * pi) / (1440)),
    e0 = testTLEs[[17]]$eccentricity,
    i0 = asteRisk:::deg2rad(testTLEs[[17]]$inclination),
    M0 = asteRisk:::deg2rad(testTLEs[[17]]$meanAnomaly),
    omega0 = asteRisk:::deg2rad(testTLEs[[17]]$perigeeArgument),
    OMEGA0 = asteRisk:::deg2rad(testTLEs[[17]]$ascension),
    Bstar = testTLEs[[17]]$Bstar,
    initialDateTime = testTLEs[[17]]$dateTime,
    targetTime = 1440
)

checkEquals(testSGDP4_1$algorithm, "sgp4")
checkEquals(testSGDP4_2$algorithm, "sdp4")
checkEqualsNumeric(testSGDP4_1$position[1], 231.5594, tolerance=0.0001)
checkEqualsNumeric(testSGDP4_2$position[1], 5501.081, tolerance=0.0001)

# Test TEMEtoECEF

testECEF <- TEMEtoECEF(testSGDP4_2$position,
                       testSGDP4_2$velocity,
                       "2006-06-27 00:58:29")

checkEqualsNumeric(testECEF$position[1], 41021.5990, tolerance=0.0001)
