#ifndef RCDPP_COMB_H
#define RCDPP_COMB_H

#include "rcdpp_All.h"
using namespace Rcpp;


/* Class of d-dimensional DPPs whose kernels is a product of 1-d kernels
  else a sum of 1-d kernels*/

class dpp_Comb : public dpp_All {

  protected:
    MATH_Kernel1D* mKernel1D=NULL;  //
    std::string mType;   // Product of sum ?

    int mK;  // Lattice

    std::vector<std::vector<double> > mEigDir;   // vectors of eigen values (one for each direction)
                                              // if domain is a cube, eigen values are identical for each direction
                                              // => only one vector

    void setEigenDir(const std::vector<double>& eigen, const int i) {
      (mEigDir[i]).resize(eigen.size());
      mEigDir[i] = eigen;
    };


    // select the eigenvalues correspondign to the i° direction
    // and store it in res
    void select(const std::vector<int>& coord, std::vector<double>& res) {

    // int n = coord.size();
      if (mIsCube) {
        for (int i = 0; i < mDim; ++i)  res[i] = (mEigDir[0])[coord[i]];
      } else {
        for (int i = 0; i < mDim; ++i)  res[i] = (mEigDir[i])[coord[i]];
      }

    } ;

  public:

    dpp_Comb() : dpp_All() {
      if (!mKernel1D) delete mKernel1D;
      mKernel1D = NULL;
      // mIsProd = true;
    };

    ~dpp_Comb() {
      if (!mKernel1D) delete mKernel1D;
    }

    dpp_Comb(List args) : dpp_All(args),
                          mK(as<int>(args["k"]))
                           {
                            if (!mKernel1D) delete mKernel1D;
                            std::string s = args["type"];
                            mType = s;
                            if (!mIsCube) mEigDir.resize(mDim);
                            else mEigDir.resize(1);
                            std::string model = args["model"];


                            if (mWithKernel) {
                              // std::cout<<"Redim mEig"<<std::endl;
                              if (mK >= 0) mEig.resize(pow(2*mK+1, mDim));
                              // std::cout<<"Fait"<<std::endl;
                            }

                            MATH_Kernel1D* K1D = NULL;
                            if (model.compare("G") == 0) {
                              // std::cout<<"Gaussian DPP"<<std::endl;
                              mIsProj = false;
                              List par;
                              par = args["param"];
                              double rho = as<double>(par["rho"]);
                              if (mType.compare("prod") == 0) par["rho"] = pow(rho, 1./mDim);
                              else if (mType.compare("sum") == 0) par["rho"] = rho/mDim;

                              // std::cout<<"rho1d ="<<as<double>(par["rho"])<<std::endl;
                              K1D = new MATH_KernelGauss1D(par);
                              // mK < 0 uniquement si calcul exact du noyau.
                              if (mK < 0) mInt = rho;

                            } else if (model.compare("L1E") == 0) {
                              mIsProj = false;
                              List par ;
                              par = args["param"];
                              double rho = as<double>(par["rho"]);
                              if (mType.compare("prod") == 0) par["rho"] = pow(rho, 1./mDim);
                              else if (mType.compare("sum") == 0) par["rho"] = rho/mDim;

                              K1D = new MATH_KernelL1Exp1D(par);
                              // mK < 0 uniquement si calcul exact du noyau.
                              if (mK < 0) mInt = rho;

                            } else K1D = NULL;

                            mKernel1D = K1D;
                            // std::cout<<"model.compare('D') ="<<model.compare("D")<<std::endl;
                            if (model.compare("D") != 0) {
                              if (mK >= 0) computeEigenDir();
                            }


                            // std::cout<<"In constructorof dpp_Prod: computation of eigenvalues"<<std::endl;
                            // computeEigenDir();

                            // std::cout<<"Construction ok."<<std::endl;
    };

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


    IntegerMatrix getIndex() {

      int n = mIndex[0].size();
      IntegerVector tp (n);
      IntegerMatrix res (n, mDim);

      for (int i=0; i < mDim; ++i) {
        tp = IntegerVector(mIndex[i].begin(), mIndex[i].end());
        res(_,i) = tp;
      }
      return res;
    };

    double methodComputationEigen (const std::vector<int>& coord, std::vector<double>& vec) ;

    virtual void computeEigenForKernel();
    virtual void computeEigenDir();

    virtual void computeIndex();

    virtual std::complex<double> computeKernel(const NumericVector& X, const NumericVector& Y);

    std::complex<double> computeExactKernel(const NumericVector& X, const NumericVector& Y) {
      std::complex<double> res (0., 0.);
      res = mKernel1D->computeExactKernel1D(X[0], Y[0]);
      // std::cout<<"current value (i=0) of res ="<<res<<std::endl;

      if (mDim > 1) {
        if (mType.compare("prod") == 0) {
          // std::cout<<"Product case"<<std::endl;
          for (int i = 1; i < mDim; ++i) {
            res *= mKernel1D->computeExactKernel1D(X[i], Y[i]);
            // std::cout<<"current value (i="<<i<<") of res ="<<res<<std::endl;
          }
        } else if (mType.compare("sum") == 0) {
          // std::cout<<"Sum case"<<std::endl;
          for (int i = 1; i < mDim; ++i) {
           res += mKernel1D->computeExactKernel1D(X[i], Y[i]);
           // std::cout<<"current value (i="<<i<<") of res ="<<res<<std::endl;
         }
        }
      }

      return res;
    };
    ComplexMatrix computeKernelR(const NumericMatrix& PP);
    // ComplexMatrix computeExactKernelR(const NumericMatrix& PP);
    // List computeOnlyExactKernelR(const List& PP) ;
};

#endif
