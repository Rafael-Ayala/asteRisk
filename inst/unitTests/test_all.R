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
checkEqualsNumeric(testSGDP4_1$position[1], 231.5655, tolerance=5e-7)
checkEqualsNumeric(testSGDP4_2$position[1], 5501.081, tolerance=2e-7)

# Test ITRFtoLATLON

# Provided values are those obtained by test of TEMEtoITRF
testLATLON1 <- ITRFtoLATLON(c(-37325973.4, 19151626.6, 138376.3))

checkEqualsNumeric(testLATLON1[1], 0.1891839, tolerance=6e-5)

# Test LATLONtoITRF

testECEF3 <- LATLONtoITRF(testLATLON1)

checkEqualsNumeric(testECEF3[1], -37325973.4, tolerance=3e-9)

# Test readGLONASSNavigationRINEX

testGLONASSnavV2 <- readGLONASSNavigationRINEX(paste0(path.package("asteRisk"), 
                                                    "/testGLONASSRINEXv2.txt"))

checkTrue(length(testGLONASSnavV2$messages) == 5)

testGLONASSnavV3 <- readGLONASSNavigationRINEX(paste0(path.package("asteRisk"), 
                                                      "/testGLONASSRINEXv3.txt"))

checkEquals(testGLONASSnavV3$messages[[2]]$tauCSource, "Ground")

# Test readGPSNavigationRINEX

testGPSnavV2 <- readGPSNavigationRINEX(paste0(path.package("asteRisk"), 
                                            "/testGPSRINEXv2.txt"))

checkTrue(length(testGPSnavV2$messages) == 3)

testGPSnavV3 <- readGPSNavigationRINEX(paste0(path.package("asteRisk"), 
                                                      "/testGPSRINEXv3.txt"))

checkEquals(testGPSnavV3$messages[[2]]$IODC, 389)

if (requireNamespace("asteRiskData", quietly = TRUE)) {
    # Test TEMEtoLATLON
    
    testLATLON2 <- TEMEtoLATLON(testSGDP4_2$position*1000,
                                "2006-06-27 00:58:29")
    
    checkEqualsNumeric(testLATLON2[1], 0.1891839, tolerance=6e-5)
    
    # Test TEMEtoITRF
    
    testECEF <- TEMEtoITRF(testSGDP4_2$position*1000,
                           testSGDP4_2$velocity*1000,
                           "2006-06-27 00:58:29")
    
    checkEqualsNumeric(testECEF$position[1], -37325973.4, tolerance=3e-9)
    
    # Test TEMEtoGCRF
    
    testGCRF <- TEMEtoGCRF(testSGDP4_2$position*1000,
                           testSGDP4_2$velocity*1000,
                           "2006-06-27 00:58:29")
    
    checkEqualsNumeric(testGCRF$position[1], 5560876.4, tolerance=2e-8)
    
    # Test ITRFtoGCRF
    
    testGCRF2 <- ITRFtoGCRF(testECEF$position,
                            testECEF$velocity,
                            "2006-06-27 00:58:29")
    
    checkEqualsNumeric(testGCRF2$position[1], 5560876.4, tolerance=2e-8)
    
    # Test GCRFtoITRF
    
    testECEF2 <- GCRFtoITRF(testGCRF$position,
                            testGCRF$velocity,
                            "2006-06-27 00:58:29")
    
    checkEqualsNumeric(testECEF2$position[1], -37325973.4, tolerance=3e-9)
    
    # Test GCRFtoLATLON
    
    testLATLON3 <- GCRFtoLATLON(testGCRF$position,
                                "2006-06-27 00:58:29")
    
    checkEqualsNumeric(testLATLON3[1], 0.1891839, tolerance=6e-5)
    
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
    
    # Test hpop 
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
    checkEqualsNumeric(testHpop[2, "positionX"], -14572000, tolerance=1e-7)
} 

## Tests of functions to read DAF, SPK

testDAF <- readBinDAF(paste0(path.package("asteRisk"), "/vgr2_jup230.bsp"))
checkEquals(testDAF$arrays[[1]]$arrayElements[1], -649364399)

testSPK1 <- readBinSPK(paste0(path.package("asteRisk"), "/vgr2_jup230.bsp"))
checkEquals(testSPK1$segments[[1]]$segmentSummary$SPKType, "Modified Difference Arrays")

testSPK2 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw02_ex1.bsp"))
checkEquals(testSPK2$segments[[1]]$segmentSummary$SPKType, "Chebyshev position (equal time steps)")

testSPK3 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw03_ex1.bsp"))
checkEquals(testSPK3$segments[[1]]$segmentSummary$SPKType, "Chebyshev position and velocity (equal time steps)")
checkEquals(dim(testSPK3$segments[[1]]$segmentData$chebyshevCoefficients), c(4, 21))

testSPK5 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw05_ex1.bsp"))
checkEquals(testSPK5[["segments"]][[1]][["segmentData"]][["stateVectors"]][[9,2]], 109)

testSPK8 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw08_ex1.bsp"))
checkEquals(testSPK8[["segments"]][[1]][["segmentData"]][["polynomialDegree"]], 3)

testSPK9 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw09_ex1.bsp"))
checkEqualsNumeric(testSPK9[["segments"]][[1]][["segmentData"]][["stateVectors"]][[3,4]], -5574.981, tolerance = 1e-6)

testSPK10 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw10_ex1.bsp"))
checkEqualsNumeric(testSPK10[["segments"]][[1]][["segmentData"]][["constants"]][["earthRadius"]], 6378.135, tolerance = 1e-6)

testSPK12 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw12_ex1.bsp"))
checkEquals(testSPK12[["segments"]][[1]][["segmentData"]][["windowSize"]], 2)

testSPK13 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw13_ex1.bsp"))
checkEqualsNumeric(diff(testSPK13[["segments"]][[1]][["segmentData"]][["stateVectors"]][,1])[2], 67.07107, tolerance = 1e-6)

testSPK14 <- readBinSPK(paste0(path.package("asteRisk"), "/spk14b_ex1.bsp"))
checkEqualsNumeric(testSPK14[["segments"]][[1]][["segmentData"]][["chebyshevCoefficients"]][[4,4]], 4.0101, tolerance = 1e-6)

testSPK15 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw15_ex1.bsp"))
checkEqualsNumeric(sqrt(sum(unlist(testSPK15[["segments"]][[1]][["segmentData"]][5:7])^2)), 1, tolerance = 1e-6)

testSPK17 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw17_ex1.bsp"))
checkEqualsNumeric(testSPK17[["segments"]][[1]][["segmentData"]][["meanLongitude"]], 2.699081, tolerance = 1e-6)

testSPK18 <- readBinSPK(paste0(path.package("asteRisk"), "/CORL_DL_001_02____A__00002.BSP"))
checkEqualsNumeric(testSPK18[["segments"]][[1]][["segmentData"]][["stateVectors"]][[13,4]], -82993704, tolerance = 1e-6)

testSPK19 <- readBinSPK(paste0(path.package("asteRisk"), "/CORL_DL_001_02_B_T19_00002.BSP"))
checkEquals(testSPK19[["segments"]][[1]][["segmentData"]][["minisegments"]][[1]][["subTypeCode"]], 1)

testSPK20 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw20_ex.bsp"))
checkEqualsNumeric(testSPK20[["segments"]][[1]][["segmentData"]][["chebyshevCoefficients"]][[3,15]], 4.21923298768136976e-15, tolerance = 1e-18)

testSPK21 <- readBinSPK(paste0(path.package("asteRisk"), "/spkw21_ex.bsp"))
checkEquals(nrow(testSPK21[["segments"]][[1]][["segmentData"]][[2]][["MDA"]]), 15)

## Test of kernel evaluators

evalSPK1ref <- matrix(c(3.63135281336917542e+06, -2.62205981397779584e+07, -8.51800527800576761e+06, -5.63779389229125738e-01, 7.81323494580312605e+00, 2.50388483296282116e+00,
                        -9.97893515585039556e+06, -4.00691830769035034e+05, 5.42991933363757678e+05, -8.92760256266241115e+00, -1.74787931490543658e+00, 4.56931508266145217e-02,
                        -1.00232005861960761e+07, -4.09360642742540163e+05, 5.43217746814036160e+05, -8.92137054840694965e+00, -1.74760903143473967e+00, 4.53613939011274364e-02,
                        -1.00236466531599630e+07, -4.09448023125729465e+05, 5.43220014800680219e+05, -8.92130800794156542e+00, -1.74760629610235507e+00, 4.53580719352145445e-02,
                        -1.00238250790699739e+07, -4.09482975240709609e+05, 5.43220921948832111e+05, -8.92128299332553354e+00, -1.74760520191039159e+00, 4.53567432730572623e-02,
                        -1.00240927169969641e+07, -4.09535403372146888e+05, 5.43222282621237100e+05, -8.92124547308360150e+00, -1.74760356055923860e+00, 4.53547504128853848e-02),
                      nrow=6, byrow=TRUE)

evalSPK1times <- c(-649364397, -645369860, -645364900, -645364850, -645364830, -645364800)

evalSPK1test <- asteRisk:::evaluateType1SPKSegment(testSPK1$segments[[1]], evalSPK1times)

evalSPK2ref <- matrix(c(1.01019999999999999e+00, 1.02020000000000000e+00, 1.03020000000000000e+00, -6.06199999999999933e-02, -6.12199999999999966e-02, -6.18200000000000069e-02,
                        -1.21415999999999968e-01, -1.22616000000000058e-01, -1.23815999999999926e-01, 4.03919999999999800e-03, 4.07919999999999811e-03, 4.11919999999999821e-03,
                        1.89118823999999974e+00, 1.90059624000000005e+00, 1.91000423999999991e+00, -1.17403519999999997e-01, -1.17987520000000026e-01, -1.18571519999999986e-01,
                        -1.99999999999977973e-04, -2.00000000000422062e-04, -1.99999999999977973e-04, 6.02040000000000075e-02, 6.04039999999999994e-02, 6.06039999999999981e-02,
                        4.01020000000000021e+00, 4.02020000000000000e+00, 4.03020000000000067e+00, -2.40619999999999973e-01, -2.41219999999999990e-01, -2.41820000000000035e-01,
                        1.20305999999999997e+01, 1.20606000000000009e+01, 1.20906000000000020e+01, 4.01027999999999996e-01, 4.02027999999999941e-01, 4.03028000000000053e-01),
                      nrow=6, byrow=TRUE)

evalSPK2times <- c(100, 140, 201, 350, 400, 500)

evalSPK2test <- asteRisk:::evaluateType2SPKSegment(testSPK2$segments[[1]], evalSPK2times)

evalSPK3ref <- matrix(c(1.01019999999999999e+00, 1.02020000000000000e+00, 1.03020000000000000e+00, 1.04020000000000001e+00, 1.05020000000000002e+00, 1.06020000000000003e+00,
                        -1.21415999999999968e-01, -1.22616000000000058e-01, -1.23815999999999926e-01, -1.25016000000000016e-01, -1.26215999999999884e-01, -1.27415999999999974e-01,
                        1.89118823999999974e+00, 1.90059624000000005e+00, 1.91000423999999991e+00, 1.91941223999999933e+00, 1.92882024000000052e+00, 1.93822823999999994e+00,
                        -1.99999999999977973e-04, -2.00000000000422062e-04, -1.99999999999977973e-04, -1.99999999999977973e-04, -1.99999999999977973e-04, -1.99999999999977973e-04,
                        4.01020000000000021e+00, 4.02020000000000000e+00, 4.03020000000000067e+00, 4.04020000000000046e+00, 4.05019999999999936e+00, 4.06020000000000003e+00,
                        1.20305999999999997e+01, 1.20606000000000009e+01, 1.20906000000000020e+01, 1.21205999999999996e+01, 1.21506000000000007e+01, 1.21806000000000001e+01),
                      nrow=6, byrow=TRUE)

evalSPK3times <- c(100, 140, 201, 350, 400, 500)

evalSPK3test <- asteRisk:::evaluateType2SPKSegment(testSPK3$segments[[1]], evalSPK3times)

evalSPK5ref <- matrix(c(1.01000000000000000e+02, 2.01000000000000000e+02, 3.01000000000000000e+02, 4.01000000000000000e+02, 5.01000000000000000e+02, 6.01000000000000000e+02,
                        1.01622449323048158e+02, 2.01285129088869525e+02, 3.00947808854690948e+02, -1.20175652301085080e+02, -5.33568772499322790e+02, -9.46961892697560302e+02,
                        7.73918870417931686e+01, 1.49374742098843939e+02, 2.21357597155894695e+02, -4.18905319939365108e+03, -8.47294482996910847e+03, -1.27568364605445659e+04,
                        1.40985329122283396e+01, 2.75118748824818127e+01, 4.09252168527352893e+01, 4.51422459370427191e+02, 2.57514294163681370e+03, 4.69886342390319714e+03,
                        1.04000000000000000e+02, 2.04000000000000000e+02, 3.04000000000000000e+02, 4.04000000000000000e+02, 5.04000000000000000e+02, 6.04000000000000000e+02,
                        1.09000000000000000e+02, 2.09000000000000000e+02, 3.09000000000000000e+02, 4.09000000000000000e+02, 5.09000000000000000e+02, 6.09000000000000000e+02),
                      nrow=6, byrow=TRUE)

evalSPK5times <- c(100, 140, 201, 350, 400, 900)

evalSPK5test <- asteRisk:::evaluateType5SPKSegment(testSPK5$segments[[1]], evalSPK5times)

evalSPK8ref <- matrix(c(1.01000000000000000e+02, 2.01000000000000000e+02, 3.01000000000000000e+02, 4.01000000000000000e+02, 5.01000000000000000e+02, 6.01000000000000000e+02,
                        1.01400000000000020e+02, 2.01400000000000006e+02, 3.01399999999999977e+02, 4.01399999999999920e+02, 5.01399999999999920e+02, 6.01400000000000091e+02,
                        1.02009999999999991e+02, 2.02009999999999991e+02, 3.02009999999999991e+02, 4.02009999999999991e+02, 5.02010000000000048e+02, 6.02009999999999877e+02,
                        1.03500000000000000e+02, 2.03500000000000000e+02, 3.03500000000000000e+02, 4.03500000000000000e+02, 5.03500000000000000e+02, 6.03500000000000000e+02,
                        1.04000000000000000e+02, 2.04000000000000000e+02, 3.04000000000000000e+02, 4.04000000000000000e+02, 5.04000000000000000e+02, 6.04000000000000000e+02,
                        1.09000000000000000e+02, 2.09000000000000000e+02, 3.09000000000000000e+02, 4.09000000000000000e+02, 5.09000000000000000e+02, 6.09000000000000000e+02),
                      nrow=6, byrow=TRUE)

evalSPK8times <- c(100, 140, 201, 350, 400, 900)

evalSPK8test <- asteRisk:::evaluateType8SPKSegment(testSPK8$segments[[1]], evalSPK8times)

evalSPK9ref <- matrix(c(4.81639065212962578e+03, -5.99882801237895819e+03, -5.55637829926284121e+03, 1.54905176943043421e+00, 1.42579867905110658e+00, -1.64122966776454049e-01,
                        4.90889631914680376e+03, -5.91274126481189069e+03, -5.56572394729042844e+03, 1.53442319076060718e+00, 1.44371749987201481e+00, -1.47392928575104443e-01,
                        4.91014578423952025e+03, -5.91156548483412280e+03, -5.56584388290529387e+03, 1.53422269777183429e+00, 1.44395896071349350e+00, -1.47165629974222334e-01,
                        4.93313137605107113e+03, -5.88987278160746973e+03, -5.56801995705583613e+03, 1.53052036713366868e+00, 1.44839822903281479e+00, -1.42977774037516298e-01,
                        -7.07793814943989491e+03, -6.20213148198037288e+03, 9.01846502150108904e+02, 1.02241342598092766e+00, -1.38721406892399957e+00, -1.23461929311409224e+00,
                        -7.06743266379557281e+03, -6.21634319567359307e+03, 8.89180957742123041e+02, 1.02609170847248743e+00, -1.38398472896571700e+00, -1.23508497862034705e+00),
                      nrow=6, byrow=TRUE)

evalSPK9times <- c(576020169.18565690517, 576020229.18565690517, 576020230, 576020245, 576068106, 576068116.25672471523)

evalSPK9test <- asteRisk:::evaluateType8SPKSegment(testSPK9$segments[[1]], evalSPK9times)

evalSPK10ref <- matrix(c(-6.11437427560608012e+03, 2.82963888361039153e+03, -2.65969245947526952e+03, -1.92233356409824707e+00, 2.23225970729102174e+00, 6.80137557250271740e+00,
                         -6.40309073901097236e+03, 3.38131047705330411e+03, -7.69857031050769081e+00, 5.31178177627154247e-01, 1.00166556861049560e+00, 7.32839972781822535e+00,
                         -6.96899792185227579e+03, 1.95494910972087905e+03, -1.85408855087806614e+01, 3.03362758253213516e-01, 1.09502585074905490e+00, 7.33115995576337998e+00,
                         2.98389301080913924e+02, 1.06473144273454295e+03, 7.12699041784661858e+03, 7.16674614552243483e+00, -1.98481262083497279e+00, -1.14465893103163074e-02,
                         -7.21693255591431625e+03, 4.72343829889546669e+02,-8.48129914191757628e+00, 8.47241913624728060e-02, 1.13561399183107725e+00, 7.33663639458283168e+00,
                         -6.67476620735959114e+03, 7.92743534012259261e+02, 2.64860514091211326e+03, 2.83161710848834947e+00, 8.97627546851406444e-01, 6.81160237805273994e+00),
                      nrow=6, byrow=TRUE)

evalSPK10times <- c(-382353604.79526937008, -382310404.79526937008, -381245795.91127812862, -381225911.01041328907, -380199544.05941164494, -380156344.05941164494)

evalSPK10test <- asteRisk:::evaluateType10SPKSegment(testSPK10$segments[[1]], evalSPK10times)

evalSPK12ref <- matrix(c(1.01000000000000000e+02, 2.01000000000000000e+02, 3.01000000000000000e+02, 4.01000000000000000e+02, 5.01000000000000000e+02, 6.01000000000000000e+02,
                         2.01655200000000104e+03, 2.59655200000000332e+03, 3.17655200000000332e+03, -1.76745599999999939e+02, -2.20745599999999939e+02, -2.64745599999999911e+02,
                         4.92010798000000364e+02, 6.89030798000000345e+02, 8.86050798000000441e+02, 3.78102094000000022e+02, 4.72162094000000025e+02, 5.66222093999999970e+02,
                         9.10000000000000000e+01, 1.91000000000000000e+02, 2.91000000000000000e+02, -2.01735000000000014e+02, -2.51735000000000014e+02, -3.01735000000000014e+02,
                         1.04000000000000000e+02, 2.04000000000000000e+02, 3.04000000000000000e+02, 4.04000000000000000e+02, 5.04000000000000000e+02, 6.04000000000000000e+02,
                         1.09000000000000000e+02, 2.09000000000000000e+02, 3.09000000000000000e+02, 4.09000000000000000e+02, 5.09000000000000000e+02, 6.09000000000000000e+02),
                       nrow=6, byrow=TRUE)

evalSPK12times <- c(100, 140, 201, 350, 400, 900)

evalSPK12test <- asteRisk:::evaluateType8SPKSegment(testSPK12$segments[[1]], evalSPK12times)

evalSPK13ref <- matrix(c(4.81639065212962578e+03, -5.99882801237895819e+03, -5.55637829926284212e+03, 1.54905176943043421e+00, 1.42579867905110658e+00, -1.64122966776454049e-01,
                         -6.66004611789999308e+03, -6.71064442757857614e+03, 4.25097524756124585e+02, 1.15590893997463762e+00, -1.26164265108543749e+00, -1.24743754267069207e+00,
                         -6.61938355904242144e+03, -6.75459375412028567e+03, 3.81424490725728845e+02, 1.16765326132798131e+00, -1.24973508475215156e+00, -1.24815148689115030e+00,
                         -6.57713178216069537e+03, -6.79936260120615862e+03, 3.36479027610012167e+02, 1.17965437130948647e+00, -1.23741003691396800e+00, -1.24880503983723679e+00,
                         -7.07621310482570516e+03, -6.20447087865108733e+03, 8.99763974013096686e+02, 1.02301870211446921e+00, -1.38668350623686232e+00, -1.23469632893220638e+00,
                         -7.06743266379557281e+03, -6.21634319567359307e+03, 8.89180957742123041e+02, 1.02609170847248743e+00, -1.38398472896571700e+00, -1.23508497862034705e+00),
                       nrow=6, byrow=TRUE)

evalSPK13times <- c(576020169.18565690517, 576040936.25672471523, 576040971.25672483444, 576041007.25672471523, 576068107.68672466278, 576068116.25672471523)

evalSPK13test <- asteRisk:::evaluateType8SPKSegment(testSPK13$segments[[1]], evalSPK13times)

evalSPK14ref <- matrix(c(1.01019999999999999e+00, 1.02020000000000000e+00, 1.03020000000000000e+00, 1.04020000000000001e+00, 1.05020000000000002e+00, 1.06020000000000003e+00,
                         -1.21415999999999968e-01,-1.22616000000000058e-01,-1.23815999999999926e-01,-1.25016000000000016e-01,-1.26215999999999884e-01,-1.27415999999999974e-01,
                         1.89118823999999974e+00, 1.90059624000000005e+00, 1.91000423999999991e+00, 1.91941223999999933e+00, 1.92882024000000052e+00, 1.93822823999999994e+00,
                         -1.99999999999977973e-04,-2.00000000000422062e-04,-1.99999999999977973e-04,-1.99999999999977973e-04,-1.99999999999977973e-04,-1.99999999999977973e-04,
                         4.01020000000000021e+00, 4.02020000000000000e+00, 4.03020000000000067e+00, 4.04020000000000046e+00, 4.05019999999999936e+00, 4.06020000000000003e+00,
                         1.20305999999999997e+01, 1.20606000000000009e+01, 1.20906000000000020e+01, 1.21205999999999996e+01, 1.21506000000000007e+01, 1.21806000000000001e+01),
                      nrow=6, byrow=TRUE)

evalSPK14times <- c(100, 140, 201, 350, 400, 500)

evalSPK14test <- asteRisk:::evaluateType2SPKSegment(testSPK14$segments[[1]], evalSPK14times)

evalSPK15ref <- matrix(c(-8.29865424318541045e+02,-3.23532908176456885e+03,-8.65822321183718486e+03,-4.98019116398637385e-01,-1.94289670737570486e+00, 8.00452146864809810e-01,
                         -8.44927195955275238e+02,-3.29411895149287238e+03,-8.63375739071445241e+03,-4.96657363315205047e-01,-1.93759568375504099e+00, 8.14492196053350548e-01,
                         -2.14509430965898582e+03,-8.37254013424238656e+03,-3.24577207989533690e+03,-1.89381656536035670e-01,-7.39226175056724499e-01, 2.03201031076371530e+00,
                         -2.14509430966461878e+03,-8.37254013426441816e+03,-3.24577207983478002e+03,-1.89381656532556453e-01,-7.39226175043144806e-01, 2.03201031076897998e+00,
                         -2.28052753988870063e+03,-8.90398495023613395e+03, 1.01914050097237237e+03, 5.54809480768173602e-02, 2.16811489618588205e-01, 2.15542471093481236e+00,
                         -2.03761043926578077e+03,-7.95689731224076240e+03, 4.31935052247881413e+03, 2.43757458090912427e-01, 9.52307501275065582e-01, 1.92283307124945635e+00),
                       nrow=6, byrow=TRUE)

evalSPK15times <- c(testSPK15$segments[[1]]$segmentSummary$initialEpoch, 251810830.3837097168, 251814400.08370968699, 251814400.0837097168, 251816400.0837097168, 251818000.0837097168)

evalSPK15test <- asteRisk:::evaluateType15SPKSegment(testSPK15$segments[[1]], evalSPK15times)

evalSPK17ref <- matrix(c(1.35916046708141279e+05, 4.11844398744404043e+03,-1.19928172762379872e+04,-6.08731562614768906e-01, 1.66572690289505800e+01,-1.17468913479814652e+00,
                         -1.34190630113009189e+05, 2.30007818215015905e+04, 9.84675467196407408e+03,-2.70972586845232533e+00,-1.64256688471936734e+01, 1.44302848555487762e+00,
                         9.43121942873224179e+04,-9.86846737847309269e+04,-8.41781872077903699e+02, 1.19958131092612827e+01, 1.14798984847375927e+01,-1.87735327466313473e+00,
                         1.05577615405620920e+05,-8.64885653573769814e+04,-2.70907685959710579e+03, 1.04956163578858455e+01, 1.28696463009662896e+01,-1.85070363006930894e+00,
                         -1.79353528795612765e+04,-1.34832130152131896e+05, 1.14748231793474588e+04, 1.65036151752925235e+01,-2.30154565951702672e+00,-1.24985000124509771e+00,
                         -1.01892824687807606e+04,-1.35687224835073052e+05, 1.08716132284716168e+04, 1.66022274945964163e+01,-1.35303363319574066e+00,-1.32820219757062752e+00),
                       nrow=6, byrow=TRUE)

evalSPK17times <- c(504878400, 704878400.467, 1041768000, 1041769000.5, 1578657132.1700000763, 1578657600)

evalSPK17test <- asteRisk:::evaluateType17SPKSegment(testSPK17$segments[[1]], evalSPK17times)

evalSPK18ref <- matrix(c(-5.54509787214519382e+08,-2.38763632192731798e+08,-6.73662690184233636e+07,-3.52141492616413210e+00,-1.16176987328556720e+01,-5.79026706118675527e+00,
                         -5.59388593663391352e+08,-2.55878632545996606e+08,-7.59294966991903484e+07,-3.04603246025743735e+00,-1.14064598481852926e+01,-5.72907045877663812e+00,
                         -5.62590195268378854e+08,-2.68493754758671701e+08,-8.22835562820288539e+07,-2.70539885578100670e+00,-1.12471268942839462e+01,-5.68099132760671299e+00,
                         -6.33420076869108826e+07, 1.77678313560346484e+08, 9.16534203512977958e+07,-3.17279400780320238e+01, 1.50109074899400230e-01, 1.49506361704912827e+00,
                         -6.94580986262519062e+07, 1.77660016547494560e+08, 9.19180075845600963e+07,-3.15459249113666793e+01,-3.37054897140248000e-01, 1.24338941518980062e+00,
                         -7.06243323894672394e+07, 1.77645852778359026e+08, 9.19631236015359163e+07,-3.15098097940438926e+01,-4.28664073814749469e-01, 1.19597861875049993e+00),
                       nrow=6, byrow=TRUE)

evalSPK18times <- c(328795200, 330281878.96174395084, 331395614.13873040676, 694036894.09602458954, 694230209.46981406212, 694267200)

evalSPK18test <- asteRisk:::evaluateType18SPKSegment(testSPK18$segments[[1]], evalSPK18times)

evalSPK19test <- asteRisk:::evaluateType19SPKSegment(testSPK19$segments[[1]], evalSPK18times)

evalSPK20ref <- matrix(c(1.39928333774887193e+07, 3.87450839020676538e+07, 1.92531139071464241e+07,-5.60412639875211838e+01, 1.29893173118727674e+01, 1.27813672830321590e+01,
                         -1.09790547406707034e+07, 3.92720684596228302e+07, 2.21400849410853125e+07,-5.71478738666573207e+01,-1.03985753062768111e+01, 4.12240482255438856e-01,
                         -2.28892338036810197e+07, 3.58829417017747462e+07, 2.15737223009344712e+07,-5.26544083167211312e+01,-2.06989758642910431e+01,-5.55491372791246096e+00,
                         -3.35521484454394542e+07, 3.04723233413670138e+07, 1.97981650279176012e+07,-4.57529494639037466e+01,-2.90388133382793079e+01,-1.07267091317408365e+01,
                         -3.35522157022484019e+07, 3.04722806542674378e+07, 1.97981492596289143e+07,-4.57528964303632861e+01,-2.90388626152653266e+01,-1.07267409684175057e+01,
                         -5.79761340734486952e+07,-2.22395848348223139e+06, 4.89387241236208100e+06,-1.04640995152464011e+01,-4.15975264229944628e+01,-2.11119016321740460e+01),
                       nrow=6, byrow=TRUE)

evalSPK20times <- c(-6699844800, -6699412800, -6699196800, -6698980800, -6698980798.529999733, -6698116800)

evalSPK20test <- asteRisk:::evaluateType20SPKSegment(testSPK20$segments[[1]], evalSPK20times)

evalSPK21ref <- matrix(c(3.63135450470729545e+06,-2.62206215794820711e+07,-8.51801278966005333e+06,-5.63779360411658992e-01, 7.81323449501207890e+00, 2.50388469111515688e+00,
                         3.63132660879073991e+06,-2.62202349790410697e+07,-8.51788889693441987e+06,-5.63779835731557988e-01, 7.81324193019232016e+00, 2.50388703069973007e+00,
                         3.63135394092793018e+06,-2.62206137662474997e+07,-8.51801028577533923e+06,-5.63779370017462034e-01, 7.81323464527566980e+00, 2.50388473839768011e+00,
                         3.63135336023546522e+06,-2.62206057186197378e+07,-8.51800770677529648e+06,-5.63779379911454215e-01, 7.81323480004718629e+00, 2.50388478709868911e+00,
                         3.63135223420620011e+06,-2.62205901133468710e+07,-8.51800270579857938e+06,-5.63779399097096956e-01, 7.81323510016717471e+00, 2.50388488153561495e+00,
                         3.63132687950709229e+06,-2.62202387308112085e+07,-8.51789009925333969e+06,-5.63779831118584762e-01, 7.81324185803644422e+00, 2.50388700799484454e+00),
                       nrow=6, byrow=TRUE)

evalSPK21times <- c(-649364400, -649364350.51981902122, -649364399, -649364397.97, -649364395.97271299362, -649364351)

evalSPK21test <- asteRisk:::evaluateType1SPKSegment(testSPK21$segments[[1]], evalSPK21times)

checkEqualsNumeric(evalSPK1test[,2:7], evalSPK1ref)
checkEqualsNumeric(evalSPK2test[,2:7], evalSPK2ref)
checkEqualsNumeric(evalSPK3test[,2:7], evalSPK3ref)
checkEqualsNumeric(evalSPK5test[,2:7], evalSPK5ref, tolerance = 2e-7)
checkEqualsNumeric(evalSPK8test[,2:7], evalSPK8ref)
checkEqualsNumeric(evalSPK9test[,2:7], evalSPK9ref)
checkEqualsNumeric(evalSPK10test[,2:7], evalSPK10ref*1000)
checkEqualsNumeric(evalSPK12test[,2:7], evalSPK12ref)
checkEqualsNumeric(evalSPK13test[,2:7], evalSPK13ref)
checkEqualsNumeric(evalSPK14test[,2:7], evalSPK14ref)
checkEqualsNumeric(evalSPK15test[,2:7], evalSPK15ref)
checkEqualsNumeric(evalSPK17test[,2:7], evalSPK17ref)
checkEqualsNumeric(evalSPK18test[,2:7], evalSPK18ref)
checkEqualsNumeric(evalSPK19test[,2:7], evalSPK18ref)
checkEqualsNumeric(evalSPK20test[,2:7], evalSPK20ref)
checkEqualsNumeric(evalSPK21test[,2:7], evalSPK21ref)
