#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
using namespace Rcpp;
using Rcpp::as;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::Dynamic;
using Eigen::AngleAxisd;
//typedef Eigen::Matrix<double, 3, Dynamic> Matrix3Colsd;
//typedef Map<Matrix3Colsd> MapMatrix3Colsd;
//typedef Eigen::Matrix< double, 6, 1 > 	Vector6d;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
double angularSeparation(NumericVector x, NumericVector y) {
    // using standard acos instead of SPICE's VSEP routine
    // tests give same results as written in SPICE example
    // using Eigen library
    Map<VectorXd> xEigen(as<Map<VectorXd> >(x));
    Map<VectorXd> yEigen(as<Map<VectorXd> >(y));
    double angle;
    double magX = xEigen.norm();
    double magY = yEigen.norm();
    if(magX == 0. || magY == 0.) {
        angle = 0;
        return angle;
    }
    VectorXd xEigenUnit = xEigen/magX;
    VectorXd yEigenUnit = yEigen/magY;
    double xyDot = xEigenUnit.dot(yEigenUnit);
    angle = acos(xyDot);
    return angle;
}

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericMatrix rotate3DVectorsAngleAxis(NumericMatrix x, NumericVector axis, double angle) {
    // Note: x should be vectors to rotate, each column is 1 vector
    Map<MatrixXd> xEigen(as<Map<MatrixXd> >(x));
    //Map<Matrix3Colsd> xEigen3Cols(xEigen, xEigen.size());
    Map<VectorXd> axisEigen(as<Map<VectorXd> >(axis));
    Matrix3d rotation;
    rotation = AngleAxisd(angle, axisEigen.normalized()).toRotationMatrix();
    Eigen::Matrix<double, 3, Dynamic> rotatedVectors;
    rotatedVectors = rotation * xEigen.transpose();
    return wrap(rotatedVectors.transpose());
}

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector rotate3DVectorAngleAxis(NumericVector x, NumericVector axis, double angle) {
    // Note: x should be 1 single vector to rotate
    Map<VectorXd> xEigen(as<Map<VectorXd> >(x));
    Map<VectorXd> axisEigen(as<Map<VectorXd> >(axis));
    Matrix3d rotation;
    rotation = AngleAxisd(angle, axisEigen.normalized()).toRotationMatrix();
    Vector3d rotatedVector;
    rotatedVector = rotation * xEigen;
    return wrap(rotatedVector);
}

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector dvhat(NumericVector position, NumericVector velocity) {
    Map<VectorXd> positionEigen(as<Map<VectorXd> >(position));
    Map<VectorXd> velocityEigen(as<Map<VectorXd> >(velocity));
    Vector3d unitVector;
    Vector3d unitVectorRate;
    Vector3d parallelComponentVelocity;
    Vector3d orthComponentVelocity;
    unitVector = positionEigen.normalized();
    if(unitVector.norm() == 0) {
        unitVectorRate = velocityEigen;
    } else {
        parallelComponentVelocity = velocityEigen.dot(unitVector) * unitVector;
        orthComponentVelocity = velocityEigen - parallelComponentVelocity;
        unitVectorRate = orthComponentVelocity / positionEigen.norm();
    }
    VectorXd result(6);
    result << unitVector, unitVectorRate;
    return wrap(result);
}

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector vectorCrossProduct3D_cpp(NumericVector u, NumericVector v) {
    if(v.length() != 3 || u.length() != 3) {
        stop("Input vectors must be length 3");
    }
    NumericVector result = {u[1]*v[2] - u[2]*v[1],
                            u[2]*v[0] - u[0]*v[2],
                            u[0]*v[1] - u[1]*v[0]};
    return result;
}
