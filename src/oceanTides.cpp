#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
List serialOceanTidesCorrections(StringVector rowTideNames,
                           NumericMatrix tideCnmSnmCoefficients,
                           NumericVector doodsonVars,
                           NumericMatrix Cnm,
                           NumericMatrix Snm,
                           double m1, double m2) {
    NumericMatrix correctedCnm = clone(Cnm);
    NumericMatrix correctedSnm = clone(Snm);
    NumericVector nvalues = tideCnmSnmCoefficients(_, 0);
    int max_degree_Cnm = Cnm.nrow()-1;
    int nrows = tideCnmSnmCoefficients.nrow();
    String newTideName = rowTideNames[0];
    String previousTideName;
    double thetaF = 0;
    int doodsonCoeff;
    double doodsonVar;
    for(int i = 0; i < 6; i++) {
        doodsonCoeff = tideCnmSnmCoefficients(0, i+2);
        doodsonVar = doodsonVars[i];
        thetaF += doodsonCoeff * doodsonVar;
    }
    double cosThetaF = cos(thetaF);
    double sinThetaF = sin(thetaF);
    double dCnmTide;
    double dSnmTide;
    int n;
    int m;
    for(int row = 0; row < nrows; row++) {
        n = tideCnmSnmCoefficients(row, 0);
        m = tideCnmSnmCoefficients(row, 1);
        if(n > max_degree_Cnm) {
            continue;
        }
        previousTideName = newTideName;
        newTideName = rowTideNames[row];
        String test_s = newTideName;
        if(newTideName != previousTideName) {
            thetaF = 0;
            for(int i = 0; i < 6; i++) {
                doodsonCoeff = tideCnmSnmCoefficients(row, i+2);
                doodsonVar = doodsonVars[i];
                thetaF += (doodsonCoeff * doodsonVar);
            }
            cosThetaF = cos(thetaF);
            sinThetaF = sin(thetaF);
        }
        dCnmTide = 1e-11*((tideCnmSnmCoefficients(row, 8) + tideCnmSnmCoefficients(row, 10))*cosThetaF +
            (tideCnmSnmCoefficients(row, 9) + tideCnmSnmCoefficients(row, 11))*sinThetaF);
        if(tideCnmSnmCoefficients(row, 1) != 0) {
            dSnmTide = 1e-11*((tideCnmSnmCoefficients(row, 10) - tideCnmSnmCoefficients(row, 8))*sinThetaF +
                (tideCnmSnmCoefficients(row, 9) - tideCnmSnmCoefficients(row, 11))*cosThetaF);
        } else {
            dSnmTide = 0;
        }
        correctedCnm(n, m) += dCnmTide;
        correctedSnm(n, m) += dSnmTide;
    }
    // Ocean pole tide: currently only for C21 and S21, but should extend to degree 10
    correctedCnm(2, 1) += -2.1778e-10*(m1 + 0.01724*m2);
    correctedSnm(2, 1) += -1.7232e-10*(m2 - 0.03365*m2);
    List output = List::create(correctedCnm, correctedSnm);
    return output;
}

// [[Rcpp::depends(RcppParallel)]]
struct OceanTidesCorrections : public Worker
{
    // source data
    const RMatrix<double> tideCnmSnmCoefficients;
    const RVector<double> rowTideNames;
    const RVector<double> doodsonVars;
    const size_t nmax;
    const size_t mmax;
    // accumulated values
    std::vector<double> dCnm;
    std::vector<double> dSnm;
    // constructors
    // Initialization constructor
    OceanTidesCorrections (
        const NumericMatrix tideCnmSnmCoefficients_in,
        const NumericVector rowTideNames_in,
        const NumericVector doodsonVars_in,
        const size_t nmax_in,
        const size_t mmax_in) :
        tideCnmSnmCoefficients(tideCnmSnmCoefficients_in),
        rowTideNames(rowTideNames_in), doodsonVars(doodsonVars_in),
        nmax(nmax_in), mmax(mmax_in), dCnm(), dSnm()
    {
        dCnm.resize((nmax+1)*(mmax+1), 0.0);
        dSnm.resize((nmax+1)*(mmax+1), 0.0);
    }
    // SPLIT constructor
    OceanTidesCorrections (
        const OceanTidesCorrections &oceanTidesCorrections, Split) :
        tideCnmSnmCoefficients(oceanTidesCorrections.tideCnmSnmCoefficients),
        rowTideNames(oceanTidesCorrections.rowTideNames), doodsonVars(oceanTidesCorrections.doodsonVars),
        nmax(oceanTidesCorrections.nmax), mmax(oceanTidesCorrections.mmax),
        dCnm(), dSnm()
    {
        dCnm.resize((nmax+1)*(mmax+1), 0.0);
        dSnm.resize((nmax+1)*(mmax+1), 0.0);
    }
    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end) {
        double newTideName = rowTideNames[0];
        double previousTideName;
        double thetaF = 0;
        int doodsonCoeff;
        double doodsonVar;
        for(size_t i = 0; i < 6; i++) {
            doodsonCoeff = tideCnmSnmCoefficients(0, i+2);
            doodsonVar = doodsonVars[i];
            thetaF += doodsonCoeff * doodsonVar;
        }
        double cosThetaF = cos(thetaF);
        double sinThetaF = sin(thetaF);
        double dCnmTide;
        double dSnmTide;
        int n;
        int m;
        for (size_t row = begin; row < end; row++) {
            n = tideCnmSnmCoefficients(row, 0);
            m = tideCnmSnmCoefficients(row, 1);
            if(n > nmax) {
                continue;
            }
            previousTideName = newTideName;
            newTideName = rowTideNames[row];
            if(newTideName != previousTideName) {
                thetaF = 0;
                for(size_t i = 0; i < 6; i++) {
                    doodsonCoeff = tideCnmSnmCoefficients(row, i+2);
                    doodsonVar = doodsonVars[i];
                    thetaF += (doodsonCoeff * doodsonVar);
                }
                cosThetaF = cos(thetaF);
                sinThetaF = sin(thetaF);
            }
            dCnmTide = 1e-11*((tideCnmSnmCoefficients(row, 8) + tideCnmSnmCoefficients(row, 10))*cosThetaF +
                (tideCnmSnmCoefficients(row, 9) + tideCnmSnmCoefficients(row, 11))*sinThetaF);
            if(tideCnmSnmCoefficients(row, 1) != 0) {
                dSnmTide = 1e-11*((tideCnmSnmCoefficients(row, 10) - tideCnmSnmCoefficients(row, 8))*sinThetaF +
                    (tideCnmSnmCoefficients(row, 9) - tideCnmSnmCoefficients(row, 11))*cosThetaF);
            } else {
                dSnmTide = 0;
            }
            dCnm[n*(mmax+1) + m] += dCnmTide;
            dSnm[n*(mmax+1) + m] += dSnmTide;
        }
    }
    void join(const OceanTidesCorrections &rhs) {
        for (size_t i = 0; i < ((nmax+1)*(mmax+1)); i++) {
            dCnm[i] += rhs.dCnm[i];
            dSnm[i] += rhs.dSnm[i];
        }
    }
};

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
List parallelOceanTidesCorrections(NumericVector rowTideNames,
                                   NumericMatrix tideCnmSnmCoefficients,
                                   NumericVector doodsonVars,
                                   NumericMatrix Cnm,
                                   NumericMatrix Snm,
                                   double m1, double m2) {
    NumericMatrix correctedCnm = clone(Cnm);
    NumericMatrix correctedSnm = clone(Snm);
    NumericVector nvalues = tideCnmSnmCoefficients(_, 0);
    const int max_degree_Cnm = Cnm.nrow() - 1;
    const int max_degree_corrections = max(nvalues);
    const size_t nmax = std::min(max_degree_Cnm, max_degree_corrections);
    const size_t nrows = tideCnmSnmCoefficients.nrow();
    OceanTidesCorrections oceanTidesCorrections(tideCnmSnmCoefficients,
                                                rowTideNames,
                                                doodsonVars,
                                                nmax, nmax);
    parallelReduce(0, nrows, oceanTidesCorrections);
    for(int i = 0; i < oceanTidesCorrections.dCnm.size(); i++) {
        correctedCnm(i/(nmax+1), i%(nmax+1)) += oceanTidesCorrections.dCnm[i];
        correctedSnm(i/(nmax+1), i%(nmax+1)) += oceanTidesCorrections.dSnm[i];
    }
    correctedCnm(2, 1) += -2.1778e-10*(m1 + 0.01724*m2);
    correctedSnm(2, 1) += -1.7232e-10*(m2 - 0.03365*m2);
    List output = List::create(correctedCnm, correctedSnm);
    return output;
}
