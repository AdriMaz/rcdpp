#include "rcdpp_Dir.h"

using namespace Rcpp;

void dpp_Dir::computeEigenDir() {


  // std::cout<<"In computeEigenVec"<<std::endl;
  std::vector<double> eig;
  // double tpInt = 0.;
  int ni;
  // int n1=0,n0=0;

  if (mIsCube) {

    std::vector<double> eig(mK, 1.);

    this->setEigenDir(eig, 0);

  } else {
    for (int i = 0; i < mDim; ++i) {
      ni = mN[i];
      eig.clear();
      for (int tpk = 0; tpk < mK; tpk++) {
        if (tpk < ni) {
          eig.push_back(1.);
          // n1++;
          // tpInt++;
        } else {
          eig.push_back(0.);
        }
      }

      this->setEigenDir(eig, i);

      }
    }
    // std::cout<<"computeEigenVec done."<<std::endl;
}


void dpp_Dir::computeIndex() {

    // std::cout<<"In computeIndex: call computeEigenVec"<<std::endl;
    // this->computeEigenDir();
    // std::cout<<"In computeIndex: computeEigenVec done"<<std::endl;

    std::vector<int> coord (mDim, 0);              // vector of a possible permutations of {0,...,k-1}^d initialized with 0^d
    std::vector<double> temp(mDim, 0.);                // subvector of K coresponding to coord
    std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d

    double tp = 1;
    int curk = 0;
    // select(coord, mEig, temp);
    this->select(coord, temp);
    if (mDim == 1) tp = temp[0];
    else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());

    if (mWithKernel & !mIsEigSet) {
      mEig[curk] = tp;
      // if (tp > 0) mInt+=tp;
    }

    if (tp == 1) {
      res.push_back(coord);
      // std::cout<<"Keep this coord: ";
      // print_vector(coord);
    }
    // }

    // int max = *std::max_element(mN.begin(), mN.end());
    int max = mK-1;
    // std::cout<<"totk = "<<totk<<std::endl;
    // if (!mIsOdd) --totk;

    while (next_variation(coord.begin(), coord.end(), 0, max)) {
      ++curk;
    // while (next_variation(coord, totk)) {
      this->select(coord, temp);

      if (mDim == 1) tp = temp[0];
      else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());

      if (mWithKernel & !mIsEigSet) {
        mEig[curk] = tp;
        // if (tp > 0) ++mInt;
      }

      if (tp == 1) {
        res.push_back(coord);
        // std::cout<<"Keep this coord: ";
        // print_vector(coord);
      }
    }

    // std::cout<<"computeIndex done: size of res ="<<res.size()<<std::endl;

    this->setIndex(res);
    if (mWithKernel & !mIsEigSet) mIsEigSet = true;
    // std::cout<<"computeIndex done."<<std::endl;

}



std::complex<double> dpp_Dir::computeKernel(const NumericVector& X, const NumericVector& Y) {

  // int n = mIndextot.size();
  // std::vector<int> K;
  // int n = 2*mK;
  int minK = 0;
  int maxK = mK-1;
  std::vector<int> coord (mDim, minK);               // vector of a possible permutations of {-k,...,k}
  std::complex<double> res;
  // std::complex<double> tpcx;

  int curK = 0;
  double tpEig;

  NumericVector Z (mDim);

  int i;
  // double tp;

  for (i = 0; i < mDim; ++i) Z[i] = X[i]-Y[i];
  tpEig = mEig[curK];
  res = tpEig*computeFourierbasis(coord, Z, mWscale);
  // res += tpcx;

  while (next_variation(coord.begin(), coord.end(), minK, maxK)) {
    ++curK;
    tpEig = mEig[curK];
    res += tpEig*computeFourierbasis(coord, Z, mWscale);
    // res += tpcx;

  }
  // res.r = tp.real(); res.i = tp.imag();
  return res;

}
