#ifndef MATH_KERNEL1D_H
#define MATH_KERNEL1D_H

// #include "rcdpp_All.h"
#define _USE_MATH_DEFINES
#include <RcppArmadillo.h>

using namespace Rcpp;

/* Class of 1-d kernels*/

class MATH_Kernel1D {
  // protected:
  //   // std::vector<double> mParams;
  //   std::string mType;

  public:

  MATH_Kernel1D () {};

  // MATH_Kernel1D (std::string type) : mType(type) {};


  virtual double computeEigen1D(const int k, double wsc) = 0;

  virtual std::complex<double> computeExactKernel1D(const double x, const double y) = 0;

};

class MATH_KernelGauss1D : public MATH_Kernel1D {

  private:

    double mRho;    // Intensity
    double mAlpha;  // alpha parameter

  public:

    MATH_KernelGauss1D() : MATH_Kernel1D() {
       // mAsprod = true;
       // mIsProj = false;
     };


    MATH_KernelGauss1D(double rho, double alpha
      // , std::string type
                      ) : MATH_Kernel1D(),
                          mRho(rho),
                          mAlpha(alpha) {
                          // mAsprod = true;
                          // mIsProj = false;
                          // if (mWithKernel) mEig.resize(pow(2*mK+1, mDim));
                          // if (mK > 0) computeEigenDir();
                          // else mInt = mRho;
                         };


    // ~dpp_Gauss() { };
    double computeEigen1D(const int k, double wsc);
    std::complex<double> computeExactKernel1D(const double x, const double y);
};

class MATH_KernelL1Exp1D : public MATH_Kernel1D {

  private:

    double mRho;    // Intensity
    double mAlpha;  // alpha parameter

    // std::vector<std::vector<int> > mIndex;

    // std::vector<double> mEig;


  public:

    MATH_KernelL1Exp1D() : MATH_Kernel1D() {
      // mAsprod = true;
      // mIsProj = false;

    };

    MATH_KernelL1Exp1D(double rho, double alpha
      // , std::string type
                      ) : MATH_Kernel1D(),
                         mRho(rho),
                         mAlpha(alpha) {
      // mIsProj = false;
      // Rcout<<"computeEigenDir done\n";
    };


    // void computeEigen(const int k);
    double computeEigen1D(const int k, double wsc) ;

    std::complex<double> computeExactKernel1D(const double x, const double y) ;


};



#endif
