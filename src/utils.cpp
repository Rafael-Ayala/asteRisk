#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
List legendre(int n, int m, double angle) {
    // using recursive relationships of Eqs 16 and on from https://mediatum.ub.tum.de/doc/1365862/1365862.pdf
    // This is equivalent to Lear's algorithm. 
    // The derivative values are actually values of P' * cos(angle)
    // The angle must be latitude, measured North and South from equator of sphere
    NumericMatrix pnm(n+1, m+1);
    NumericMatrix dpnm(n+1, m+1);
    double cosTheta = cos(angle);
    double sinTheta = sin(angle);
    pnm(0,0) = 1;
    pnm(1,1) = sqrt(3) * cosTheta;
    dpnm(1,1) = -sqrt(3) * sinTheta;
    // Diagonal coefficients
    for(double i=2; i <= n; i++) {
        pnm(i,i) = sqrt((2*i+1)/(2*i))*cosTheta*pnm(i-1,i-1);
        dpnm(i,i) = sqrt((2*i+1)/(2*i))*(cosTheta*dpnm(i-1,i-1) - sinTheta*pnm(i-1,i-1));
    }
    // First horizontal step of recurrence
    for(double i=1; i <= n; i++) {
        pnm(i,i-1) = sqrt(2*i+1)*sinTheta*pnm(i-1,i-1);
        dpnm(i,i-1) = sqrt(2*i+1)*(cosTheta*pnm(i-1,i-1) + sinTheta*dpnm(i-1,i-1));
    }
    // Second horizontal step of recurrence
    int j = 0;
    int k = 2;
    while((j <= m) & (k <= n)) {
        for(double i = k; i <= n; i++) {
            pnm(i, j) = sqrt((2*i+1)/((i-j)*(i+j)))*
                ((sqrt(2*i-1)*sinTheta*pnm(i-1,j)) -
                (sqrt(((i+j-1)*(i-j-1))/(2*i-3))*pnm(i-2,j)));
            dpnm(i, j) = sqrt((2*i+1)/((i-j)*(i+j)))*
                ((sqrt(2*i-1)*sinTheta*dpnm(i-1,j)) +
                (sqrt(2*i-1)*cosTheta*pnm(i-1,j)) -
                (sqrt(((i+j-1)*(i-j-1))/(2*i-3))*dpnm(i-2,j)));
        }
        j += 1;
        k += 1;
    }
    List output = List::create(pnm, dpnm);
    return output;
}

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector gravityGradientSphericalCoords(const NumericMatrix Pnm,
                                             const NumericMatrix dPnm,
                                             const NumericMatrix Cnm,
                                             const NumericMatrix Snm,
                                             const double lat,
                                             const double lon,
                                             const double d,
                                             const double R,
                                             const double GM,
                                             const int n,
                                             const int m) {
    // Using equations A1 to A10 from "Gravitational gradient changes following 
    // the 2004 December 26 Sumatra–Andaman Earthquake inferred from GRACE"
    double dUr = 0;
    double dUlat = 0;
    double dUlon = 0;
    double g1 = -GM/pow(d,2);
    double g2 = GM/d;
    double x1;
    double x2;
    for(double i = 0 ; i <= n ; i++) {
        // for(double i = 2 ; i < n ; i++) {
        x1 = g1 * (i + 1) * pow((R/d), i);
        x2 = g2 * pow((R/d), i);
        for(double j = 0; j <= i; j++) {
            dUr += x1*Pnm(i,j)*(Cnm(i,j)*cos(j*lon) + Snm(i,j)*sin(j*lon));
            // Note the derivative with respect to latitude should include
            // multiplication by cos(lat) due to chain rule, but this is 
            // already included in the values calculated by function legendre
            // so I don't include it here
            dUlat += x2*dPnm(i,j)*(Cnm(i,j)*cos(j*lon) + Snm(i,j)*sin(j*lon));
            // dUlat += x2*(Pnm(i,j+1)-j*tan(lat)*Pnm(i,j))*(Cnm(i,j)*cos(j*lon) + Snm(i,j)*sin(j*lon));
            dUlon += x2* j * Pnm(i,j) * (Snm(i,j)*cos(j*lon) - Cnm(i,j)*sin(j*lon));
        }
    }
    NumericVector output = {dUr, dUlat, dUlon};
    return output;
}
