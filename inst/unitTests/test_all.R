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
checkEqualsNumeric(testSGDP4_1$position[1], 231.5594, tolerance=5e-7)
checkEqualsNumeric(testSGDP4_2$position[1], 5501.081, tolerance=2e-7)

# Test TEMEtoECEF

testECEF <- TEMEtoECEF(testSGDP4_2$position*1000,
                       testSGDP4_2$velocity*1000,
                       "2006-06-27 00:58:29")

checkEqualsNumeric(testECEF$position[1], -37325973.4, tolerance=3e-9)

# Test ECEFtoLATLON

testLATLON1 <- ECEFtoLATLON(testECEF$position)

checkEqualsNumeric(testLATLON1[1], 0.1891839, tolerance=6e-6)

# Test TEMEtoLATLON

testLATLON2 <- TEMEtoLATLON(testSGDP4_2$position*1000,
                            "2006-06-27 00:58:29")

checkEqualsNumeric(testLATLON2[1], 0.1891839, tolerance=6e-6)

# Test TEMEtoGCRF

testGCRF <- TEMEtoGCRF(testSGDP4_2$position*1000,
                       testSGDP4_2$velocity*1000,
                       "2006-06-27 00:58:29")

checkEqualsNumeric(testGCRF$position[1], 5560876.4, tolerance=2e-8)

# Test ECEFtoGCRF

testGCRF2 <- ECEFtoGCRF(testECEF$position,
                        testECEF$velocity,
                        "2006-06-27 00:58:29")

checkEqualsNumeric(testGCRF2$position[1], 5560876.4, tolerance=2e-8)

# Test GCRFtoECEF

testECEF2 <- GCRFtoECEF(testGCRF$position,
                        testGCRF$velocity,
                        "2006-06-27 00:58:29")

checkEqualsNumeric(testECEF2$position[1], -37325973.4, tolerance=3e-9)

# Test GCRFtoLATLON

testLATLON3 <- GCRFtoLATLON(testGCRF$position,
                            "2006-06-27 00:58:29")

checkEqualsNumeric(testLATLON3[1], 0.1891839, tolerance=6e-6)

# Test ECItoKOE

testKOE <- ECItoKOE(testGCRF$position,
                    testGCRF$velocity)

checkEqualsNumeric(testKOE$argumentPerigee[1], 5.4, tolerance=0.02)

# Test KOEtoECI

testECI <- KOEtoECI(a=testKOE$semiMajorAxis,
                    e=testKOE$eccentricity,
                    i=testKOE$inclination,
                    M=testKOE$meanAnomaly,
                    omega=testKOE$argumentPerigee,
                    OMEGA=testKOE$longitudeAscendingNode)

checkEqualsNumeric(testECI$position[1], 5560876.4, tolerance=2e-8)

# Test readGLONASSNavigationRINEX

testGLONASSnav <- readGLONASSNavigationRINEX(paste0(path.package("asteRisk"), 
                                                    "/testGLONASSRINEX.txt"))

checkTrue(length(testGLONASSnav$messages) == 5)

# Test readGPSNavigationRINEX

testGPSnav <- readGPSNavigationRINEX(paste0(path.package("asteRisk"), 
                                            "/testGPSRINEX.txt"))

checkTrue(length(testGPSnav$messages) == 3)

# Test hpop - only if asteRiskData available

if (requireNamespace("asteRiskData", quietly = TRUE)) {
    initialPosition <-c(-14568679.5026116, -4366250.78287623, 9417.9289105405)
    initialVelocity <- c(-3321.17428902497, -3205.49400830455, 4009.26862308806) 
    initialTime <- "2006-06-25 00:33:43"
    molniyaMass <- 1600
    molniyaCrossSection <- 15
    molniyaCr <- 1.2
    molniyaCd <- 2.2
    targetTimes <- c(0, 1)
    
    testHpop <- hpop(initialPosition, initialVelocity, initialTime, targetTimes, 
                     molniyaMass, molniyaCrossSection, molniyaCrossSection,
                     molniyaCr, molniyaCd)
    
    checkEqualsNumeric(testHpop[2, "X"], -14572000, tolerance=1e-6)
} 
