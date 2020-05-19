#ifndef RCDPP_PROD_H
#define RCDPP_PROD_H


#include "rcdpp_All.h"

using namespace Rcpp;

// Class of d-dimensional DPPs whose each eigenvalue (indexed by Z^d) is product of values indexed by Z
// contains: Gaussian, L1Exponential and MRProd DPPs
class dpp_Prod : public dpp_All {

  protected:
    int mK;           // Lattice

    std::vector<std::vector<double> > mEigDir;

    void setEigenDir(const std::vector<double>& eigen, const int i) {
      (mEigDir[i]).resize(eigen.size());
      mEigDir[i] = eigen;
    };


    void select(const std::vector<int>& coord, std::vector<double>& res) {

    // int n = coord.size();
      if (mIsCube) {
        for (int i = 0; i < mDim; ++i)  res[i] = (mEigDir[0])[coord[i]];
      } else {
        for (int i = 0; i < mDim; ++i)  res[i] = (mEigDir[i])[coord[i]];
      }

    } ;

  public:

    dpp_Prod() : dpp_All() {
      // mIsProd = true;
    };

    dpp_Prod(List args) : dpp_All(args),
                          mK(as<int>(args["k"])) {
                            if (!mIsCube) mEigDir.resize(mDim);
                            else mEigDir.resize(1);
                            // std::cout<<"In constructorof dpp_Prod: computation of eigenvalues"<<std::endl;
                            // computeEigenDir();
    };

    // dpp_Prod(double rho, int dim) : dpp_All(rho, dim) {
      // mIsProd = true;
      // if (!mIsCube) mEig.resize(mDim);
      // else mEig.resize(1);
    // };

    // ~dpp_Prod() { };
//

    // NumericMatrix getEigen() {
    //   // computeEigenDir();
    //   int n = mEigDir[0].size();
    //   NumericVector tp (n);
    //   NumericMatrix res (n, mDim);
    //
    //   if (mIsCube) {
    //     tp = NumericVector(mEigDir[0].begin(), mEigDir[0].end());
    //     for (int i=0; i < mDim; ++i)  res(_,i) = tp;
    //   } else {
    //     for (int i=0; i < mDim; ++i) {
    //       tp = NumericVector(mEigDir[i].begin(), mEigDir[i].end());
    //       res(_,i) = tp;
    //     }
    //   }
    //   return res;
    // };

    NumericMatrix getEigen() {
      // computeEigenDir();
      int n = mEigDir[0].size();
      NumericVector tp (n);
      NumericMatrix res (n, mDim);

      if (mIsCube) {
        tp = NumericVector(mEigDir[0].begin(), mEigDir[0].end());
        for (int i=0; i < mDim; ++i)  res(_,i) = tp;
      } else {
        for (int i=0; i < mDim; ++i) {
          tp = NumericVector(mEigDir[i].begin(), mEigDir[i].end());
          res(_,i) = tp;
        }
      }
      return res;
    };

    virtual void computeEigenForKernel();
    virtual void computeEigenDir();

    virtual void computeIndex();

    virtual double computeEigen(const int k, int i) = 0;

    virtual std::complex<double> computeExactKernel(const NumericVector& X, const NumericVector& Y) = 0;
    virtual std::complex<double> computeKernel(const NumericVector& X, const NumericVector& Y);
    ComplexMatrix computeKernelR(const NumericMatrix& PP);


};

class dpp_Gauss : public dpp_Prod {

  private:

    double mRho;    // Intensity
    double mAlpha;  // alpha parameter


    // std::vector<std::vector<int> > mIndex;
    //
    // std::vector<double> mEig;


  public:

    dpp_Gauss() : dpp_Prod() {
       // mAsprod = true;
       mIsProj = false;
     };


    dpp_Gauss(List args) : dpp_Prod(args),
                          mRho(as<double>(args["rho"])),
                          mAlpha(as<double>(args["alpha"])) {
                          // mAsprod = true;
                          mIsProj = false;
                          if (mWithKernel) mEig.resize(pow(2*mK+1, mDim));
                            // std::cout << "Construction of dpp_Gauss" << std::endl;
                            // std::cout << "Rho  = " << mRho << std::endl;
                            // std::cout << "Alpha = " << mAlpha << std::endl;
                            // std::cout << "Dim = " << mDim << std::endl;
                            if (mK > 0) computeEigenDir();
                            else mInt = mRho;
                         };


    // ~dpp_Gauss() { };

    double computeEigen(const int k, int i) {
      double res;
      double wsc = mWscale[i];
      res = pow(mRho, 1./mDim)*sqrt(M_PI)*mAlpha*exp(-pow(k/wsc*mAlpha*M_PI,2));
      return res;
    };

    std::complex<double> computeExactKernel(const NumericVector& X, const NumericVector& Y) {

      double norm = 0.;
      double tpres;
      for (int i =0; i < mDim; ++i) norm += pow(X[i]-Y[i], 2);

      tpres = mRho*exp(-norm/pow(mAlpha, 2));
      std::complex<double> res (tpres, 0.);
      return res;
    };

};


class dpp_L1Exp : public dpp_Prod {

  private:

    double mRho;    // Intensity
    double mAlpha;  // alpha parameter

    // std::vector<std::vector<int> > mIndex;

    // std::vector<double> mEig;


  public:

    dpp_L1Exp() : dpp_Prod() {
      // mAsprod = true;
      mIsProj = false;

    };

    dpp_L1Exp(List args) : dpp_Prod(args),
                         mRho(as<double>(args["rho"])),
                         mAlpha(as<double>(args["alpha"])) {
      mIsProj = false;
      if (mWithKernel) mEig.resize(pow(2*mK+1, mDim));
      // Rcout<<"Call computeEigenDir\n";
      if (mK > 0) computeEigenDir();
      else mInt = mRho;
      // Rcout<<"computeEigenDir done\n";
    };


    // dpp_L1Exp(double rho, double alpha, int dim) : dpp_Prod(rho, dim),
    //                                              mAlpha(alpha) {
    //   // mAsprod = true;
    //   mIsProj = false;
    //
    //  };

    // ~dpp_L1Exp() { };


    // void computeEigen(const int k);
    double computeEigen(const int k, int i) {
      double res;
      double wsc = mWscale[i];
      res = pow(mRho, 1./mDim)*2*mAlpha/(1.+pow(2*M_PI*mAlpha*k/wsc, 2));
      return res;
    };

    std::complex<double> computeExactKernel(const NumericVector& X, const NumericVector& Y) {

      double norm = 0.;
      double tpres;
      for (int i =0; i < mDim; ++i) norm += std::abs(X[i]-Y[i]);
      tpres = mRho*exp(-norm/mAlpha);
      std::complex<double> res (tpres, 0.);
      return res;
    };


};


#endif
