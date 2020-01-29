#ifndef mam_h
#define mam_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <list>
#define MIN(a, b) (a<b?a:b)
#define MAX(a, b) (a>b?a:b) 
#define IDEX(a, b) (a>b?1:0) 
using namespace Rcpp;
using namespace Eigen;

//----------------------------------------------------------------**
//***----------------------UpTriangularInv------------------------**
MatrixXd UpTriangularInv(MatrixXd A);
//----------------------------------------------------------------**
//***---------------------- Q*R of qr decomposition --------------**
MatrixXd QbyR(MatrixXd Q, MatrixXd R, int isQR);
//----------------------------------------------------------------**
//***--------------------penalty----------------------------------**
double penalties(double z, double v, double lambda, double alpha, double gamma, int penalty);

#endif