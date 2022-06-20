#include "MATH_Kernel1D.h"


/// Computation for Gaussian 1-d kernels
double MATH_KernelGauss1D::computeEigen1D(const int k, double wsc) {
  double res;
  // double wsc = mWscale[i];
  // if (k == 0) {
  //   if (mType.compare("sum") == 0) return 0.;
  // }
  res = mRho*sqrt(M_PI)*mAlpha*exp(-pow(k/wsc*mAlpha*M_PI,2));
  return res;
}

std::complex<double> MATH_KernelGauss1D::computeExactKernel1D(const double x, const double y) {

  // double norm = 0.;
  double tpres;

  tpres = mRho*exp(-pow((x-y)/mAlpha, 2));

  //
  // for (int i =0; i < mDim; ++i) norm += pow(X[i]-Y[i], 2);
  //
  // tpres = mRho*exp(-norm/pow(mAlpha, 2));
  std::complex<double> res (tpres, 0.);
  return res;
}

/// Computation for L1Exponential 1-d kernels
double MATH_KernelL1Exp1D::computeEigen1D(const int k, double wsc) {
  double res;
  // double wsc = mWscale[i];
  // if (k == 0) {
  //   if (mType.compare("sum") == 0) return 0.;
  // }
  res = mRho*2*mAlpha/(1.+pow(2*M_PI*mAlpha*k/wsc, 2));
  return res;
}

std::complex<double> MATH_KernelL1Exp1D::computeExactKernel1D(const double x, const double y) {

  // double norm = 0.;
  double tpres;
  // for (int i =0; i < mDim; ++i) norm += std::abs(X[i]-Y[i]);
  tpres = mRho*exp(-std::abs(x-y)/mAlpha);
  std::complex<double> res (tpres, 0.);
  return res;
}
