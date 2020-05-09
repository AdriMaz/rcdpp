#include "rcdpp_Prod.h"

using namespace Rcpp;


void dpp_Prod::computeEigenForKernel() {


      std::vector<int> coord (mDim, 0);               // vector of a possible permutations of {-k,...,k}
      int curk = 0;


      double tp;

      std::vector<double> temp (mDim, 0.);          // vector of eigenvalues coresponding to coord
      select(coord, temp);
      if(mDim == 1) tp = temp[0];
      else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());    // Eigenvalue associated to  coord

      // if (mWithKernel & !mIsEigSet) {
      mEig[curk] = tp;
      if (tp > 0) mInt = tp;
      // }
      // std::cout << "Corresponding eigenvalue.: " << tp << std::endl;

      // if (mWithKernel) mEig.push_back(tp);
      // Repeat for all the permutations
      int totk = 2*mK;
      while (next_variation(coord.begin(), coord.end(), 0, totk)) {
        ++curk;
        select(coord, temp);
        if(mDim == 1) tp = temp[0];
        else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());

        mEig[curk] = tp;
        if (tp > 0) mInt+=tp;
      }

      if (!mIsEigSet) mIsEigSet = true;

}

void dpp_Prod::computeEigenDir() {
  // std::cout<<"In computeEigenDir"<<std::endl;
  std::vector<double> eig;

  double tpeig;
  double wsc;
  int curk=0;
  // double tpInt=0;
  // double sumEig;
  if (mIsCube) {
    // wsc = mWscale[0];
    // std::cout<<"Wscale = " << wsc << std::endl;
    for (int tpk = -mK; tpk < mK+1; tpk++)  {
      tpeig = computeEigen(tpk, 0);
      eig.push_back(tpeig);
      // if (mWithKernel) {
      //   tpeig = pow(tpeig, mDim);
      //   mEig[curk] = tpeig;
      //   tpInt += tpeig;
      //   curk++;
      // }
    }
    setEigenDir(eig, 0);
    // sumEig = std::accumulate(eig.begin(), eig.end(), 0.);
    // std::cout<<"sum(eig) = " << sumEig << std::endl;
  } else {
    double wscm = 0.;
    // double wsc = 1.;
    for (int i = 0; i < mDim; ++i) {
      wsc = mWscale[i];
      if (wsc != wscm) {   // Compute only if current wsc is different from previous one
        eig.clear();
        curk = 0;
        for (int tpk = -mK; tpk < mK+1; tpk++)  {
          tpeig = computeEigen(tpk, i);
          eig.push_back(tpeig);
          // if (mWithKernel) {
          //   mEig[curk]*=tpeig;
          //   curk++;
          // }
        }
        wscm = wsc;
      }
      // tpInt = std::accumulate(mEig.begin(), mEig.end(), 0.);
      setEigenDir(eig, i);
    }
  }
  // setInt(tpInt);
}


void dpp_Prod::computeIndex() {

  // if (mEigDir.size() == 0) computeEigenDir();

    std::vector<int> coord (mDim, 0);               // vector of a possible permutations of {-k,...,k}
    std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d

    int curk = 0;
    // double tpInt =0. ;
      // int totk = 2*k+1;


    double tp;

    std::vector<double> temp (mDim, 0.);          // vector of eigenvalues coresponding to coord
    select(coord, temp);
    if(mDim == 1) tp = temp[0];
    else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());    // Eigenvalue associated to  coord

    if (mWithKernel & !mIsEigSet) {
      mEig[curk] = tp;
      if (tp > 0) mInt = tp;
    }
    // std::cout << "Corresponding eigenvalue.: " << tp << std::endl;

    // if (mWithKernel) mEig.push_back(tp);
    double tpbool;
    if (tp > 0) {
      tpbool = rbinom(1, 1, tp)[0];  // Bernoulli draw
      if (tpbool == 1) res.push_back(coord);     // Keep the elements with proba tp
    }

    // Repeat for all the permutations
    int totk = 2*mK;
    while (next_variation(coord.begin(), coord.end(), 0, totk)) {
      ++curk;
    // while (next_variation(coord, totk)) {
      // std::cout << "Next permutation" <<std::endl;
      // ++cpt;
      // select(coord, mEig, temp);
      select(coord, temp);
      // std::cout<<"Proposed point:  ("<<temp[0]<<" , "<<temp[1]<<")"<<std::endl;
      // std::cout << "With coord.: ("<<coord[0]<<" , "<<coord[1]<<")"<<std::endl;
      if(mDim == 1) tp = temp[0];
      else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());

      if (mWithKernel & !mIsEigSet) {
        mEig[curk] = tp;
        if (tp > 0) mInt+=tp;
      }
      // std::cout << "Corresponding eigenvalue.: " << tp << std::endl;
      // tp = rbinom(1, 1, tp)[0];
      // if (tp == 1.) res.push_back(coord);
      if (tp > 0) {
        // tpInt+=tp;
        tpbool = rbinom(1, 1, tp)[0];
        // std::cout<< "Result of the draw: " << tpbool <<std::endl;
        if (tpbool == 1) res.push_back(coord);     // Keep the elements with proba tp
      }

    }
    // std::cout << "In computeIndex :  Size of the res = " << res.size() << std::endl;
    // print_vector(res);

    // std::cout << "In computeIndex :  Size of the mIndex = " << mIndex.size() << std::endl;

  // } else {    // Until now: only for 'most' repulsive stationary DPP


      // }
    setIndex(res);
    if (mWithKernel & !mIsEigSet) mIsEigSet = true;
    // setInt(tpInt);

}


ComplexMatrix dpp_Prod::computeKernelR(const NumericMatrix& PP) {

    if (!mIsEigSet & (mK > 0)) computeEigenForKernel();

    int np = PP.nrow();
    ComplexMatrix res (np,np);
    NumericVector tpX (mDim);

    Rcomplex diag, tpR;
    diag.r = mInt; diag.i = 0.;
    std::complex<double> tp;

    for (int i = 0; i < np; ++i) {
      res(i,i) = diag;
      tpX = PP(i,_);
      for (int j = 0; j < i; ++j) {
        // if (j != i) {
          // tpY = PP(j,_);
          tp = computeKernel(tpX, PP(j,_));
          tpR.r = tp.real(); tpR.i = tp.imag();
          res(i,j) = tpR;
          tpR.i *= -1;
          res(j,i) = tpR;
        // }
      }
    }

  return res;

}


std::complex<double> dpp_Prod::computeKernel(const NumericVector& X, const NumericVector& Y) {

  std::complex<double> res = 0.;
  if (mK > 0){
    int minK = -mK;
    int maxK = mK;
    std::vector<int> coord (mDim, minK);               // vector of a possible permutations of {-k,...,k}

    std::complex<double> tpcx;

    int curK = 0;
    double tpEig;

    NumericVector Z (mDim);

    int i;
    // double tp;

    for (i = 0; i < mDim; ++i) Z[i] = X[i]-Y[i];

    tpEig = mEig[curK];
    tpcx = tpEig*computeFourierbasis(coord, Z, mWscale);
    res += tpcx;

    while (next_variation(coord.begin(), coord.end(), minK, maxK)) {
      ++curK;
      tpEig = mEig[curK];
      tpcx = tpEig*computeFourierbasis(coord, Z, mWscale);
      res += tpcx;


    }
  } else {
    res = computeExactKernel(X,Y);
  }
  // res.r = tp.real(); res.i = tp.imag();
  return res;

}
