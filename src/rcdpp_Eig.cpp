#include "rcdpp_Eig.h"

using namespace Rcpp;

void dpp_Eig::computeIndex() {

  // std::vector<int> coord (mDim, 0);              // vector of a possible permutations of {-k,...,k} initialized with -k
  // std::vector<double> temp(mDim, 0.);                // subvector of K coresponding to coord

  // std::cout<<"In ComputeIndex"<<std::endl;
  if (mIsProj) {
    this->setIndex(mIndextot);
  } else {

    std::vector< std::vector<int> > res;

    // std::vector<double> temp = mEig;

    double tp, tpbool;

    for (int i = 0; i < mIndextot.size(); ++i) {
      // tp = temp[i];
      // std::cout<<"In ComputeIndex"<<std::endl;
      tp = mEig[i];
      // std::cout<<"mEig[" <<i<< "]=" << tp <<std::endl;
      tpbool = rbinom(1, 1., tp)[0];  // Bernoulli draw
      // std::cout<<"tpbool=" << tpbool <<std::endl;
      if (tpbool == 1) {
        // coord = mIndextot[i];
        res.push_back(mIndextot[i]);
      }
    }

    this->setIndex(res);
  }
}





ComplexMatrix dpp_Eig::computeKernelR(const NumericMatrix &PP) {

  int np = PP.nrow();
  ComplexMatrix res (np,np);
  NumericVector tpX (mDim);
  // NumericVector tpY (mDim);
  Rcomplex diag, tpR;
  diag.r = mInt; diag.i = 0.;
  std::complex<double> tp;


  for (int i = 0; i < np; ++i) {
    res(i,i) = diag;
    tpX = PP(i,_);
    for (int j = 0; j < i; ++j) {
      // if (j != i) {
        // tpY = PP(j,_);
        tp = computeKernel(tpX,PP(j,_));
        tpR.r = tp.real(); tpR.i = tp.imag();
        res(i,j) = tpR;
        tpR.i *= -1;
        res(j,i) = tpR;
      // }
    }
  }
  return res;
}

std::complex<double> dpp_Eig::computeKernel(const NumericVector& X, const NumericVector& Y) {

  int n = mIndextot.size();
  std::vector<int> K;
  std::complex<double> res;
  NumericVector Z (mDim);
  int i;
  for (i = 0; i < mDim; ++i) Z[i] = X[i]-Y[i];

  for (i = 0; i < n; ++i) {       // Basis-functions evaluated at first point
    K = mIndextot[i];
    res += mEig[i]*computeFourierbasis(K, Z, mWscale);
  }
  // res.r = tp.real(); res.i = tp.imag();
  return res;

}

// arma::cx_mat dpp_Eig::computeKernel(const NumericMatrix &PP) {
//
//   int np = PP.nrow();
//   arma::cx_mat res (np,np);
//   // NumericVector tpX (mDim);
//   // NumericVector tpY (mDim);
//   std::complex<double> tp;
//   std::complex<double> diag(mInt,0.);
//   // diag.r = mInt; diag.i = 0.;
//
//
//   for (int i = 0; i < np; ++i) {
//     res(i,i) = diag;
//     NumericVector tpX = PP(i,_);
//     for (int j = 0; j < i; ++j) {
//       // if (j != i) {
//         // tpY = PP(j,_);
//         tp = computeKernel(tpX,PP(j,_));
//         res(i,j) = tp;
//         res(j,i) = conj(tp);
//       // }
//     }
//   }
//   return res;
// }


// arma::mat dpp_Eig::computePCF(const NumericMatrix& PP) {
//
//   int np = PP.nrow();
//   arma::mat res (np,np);
//   arma::cx_mat K = computeKernel(PP);
//   // NumericVector tpX (mDim);
//   // NumericVector tpY (mDim);
//   std::complex<double> tpc;
//   double tp;
//   // diag.r = mInt; diag.i = 0.;
//
//
//   for (int i = 0; i < np; ++i) {
//     res(i,i) = 0.;
//     NumericVector tpX = PP(i,_);
//     for (int j = 0; j < i; ++j) {
//       // if (j != i) {
//         // tpY = PP(j,_);
//         tpc = K(i,j);
//         tp = pow(abs(tpc), 2);
//         tp /= pow(mInt,2);
//         res(i,j) = tp; res(j,i) = tp;
//       // }
//     }
//   }
//   return res;
// }
//
// NumericMatrix dpp_Eig::computePCFR(const NumericMatrix& PP) {
//
//   int np = PP.nrow();
//   NumericMatrix res (np,np);
//   ComplexMatrix K = computeKernelR(PP);
//   // NumericVector tpX (mDim);
//   // NumericVector tpY (mDim);
//   Rcomplex tpc;
//   double tp;
//   // diag.r = mInt; diag.i = 0.;
//
//
//   for (int i = 0; i < np; ++i) {
//     res(i,i) = 0.;
//     NumericVector tpX = PP(i,_);
//     for (int j = 0; j < i; ++j) {
//       // if (j != i) {
//         // tpY = PP(j,_);
//         tpc = K(i,j);
//         tp = pow(tpc.r,2)+pow(tpc.i,2);
//         tp /= pow(mInt,2);
//         res(i,j) = tp; res(j,i) = tp;
//       // }
//     }
//   }
//   return res;
//
// }
//
//
// double dpp_Eig::computekInt(const NumericMatrix& PP) {
//
//   arma::cx_mat K = computeKernel(PP);
//   std::complex<double> tp = det(K);
//   double res = tp.real();
//   return res;
// }
