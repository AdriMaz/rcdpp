#ifndef RCDPP_DIR_H
#define RCDPP_DIR_H

#include "rcdpp_Prod.h"

using namespace Rcpp;
// Dirichlet DPP

class dpp_Dir : public dpp_Prod {

  private:
    IntegerVector mN;

  public:

    dpp_Dir() : dpp_Prod() {
      // mAsprod = false;
      mIsProj = true;
    };

    dpp_Dir(List args) : dpp_Prod(args) {
      // mTau = (pow(mRho, 1./mDim)-1.)/2.;
      // std::cout<<"Tau = "<<mTau << std::endl;
      // mAsprod = false;
      // std::cout<<"mN = "<<mN<<", mIsOdd = "<<mIsOdd<<std::endl;

      // std::cout<<"Read 'N' arguments"<<std::endl;
      mN = args["N"];
      if (mWithKernel) mEig.resize(pow(mK,mDim));
      mInt = std::accumulate(mN.begin(),mN.end(), 1., std::multiplies<double>());
      double tp = std::accumulate(mWscale.begin(), mWscale.end(), 1., std::multiplies<double>());
      mInt /= tp;
      computeEigenDir();
    };


    ~dpp_Dir() { };


    void computeEigenForKernel();
    void computeIndex();

    void computeEigenDir();

    std::complex<double> computeKernel(const NumericVector& X, const NumericVector& Y) ;

    double computeEigen(const int k, int i) {
      double res;
      if (k < mN[i]) res = 1.;
      else res = 0.;
      return res;
    };


    std::complex<double> computeExactKernel(const NumericVector& X, const NumericVector& Y){

      return computeKernel(X,Y);

    };


};



#endif
