#ifndef RCDPP_ALL_H
#define RCDPP_ALL_H

// #define _USE_MATH_DEFINES     // get M_PI
#include "MATH_Kernel1D.h"
// #include <Rcpp.h>
// #include <RcppArmadillo.h>
// #include <math.h>
// #include <iostream>
// #include <vector>
// #include <random>
// #include <complex>


using namespace Rcpp;

// Mother class of DPPs
class dpp_All {

  protected:

    // double mRho;   // Intensity
    int mDim;         // Dimension
    double mInt;

    // NumericVector mBinfs;
    // NumericVector mBsups;
    NumericVector mWscale;
    NumericVector mWcenter;

    bool mIsEigSet;
    std::vector<double> mEig;

    // double mPws;        // prod(mWscale)

    int mProg;      // Progress of the i° sampling
    int mProgSim;   // Progress of the sampling

    bool mIsProj;      // Is it a Projection DPP or not?
    // bool mIsProd;      // Is it a d-dimensional DPP with a kernel which is a product of d identical one-deimensional kernels?

    bool mIsCube;    // Is the domain a cube (i.e. product of d identical intervals)?

    std::vector<std::vector<int> > mIndex;   // Indexes in the Mercer's decomposition kept after Bernoulli drawns

    bool mWithKernel;
    // std::vector<double> mEig;             // Eigenvalues (i.e. Fourier Transform computed on [-k ; k])


    void setIndex(const std::vector<std::vector<int> >& index) {
      mIndex.resize(index.size());
      mIndex = index;
    };

    void resetIndex() {
      mIndex.clear();
    };

    void setInt(const double rho) {
      mInt = rho;
    };

    // void setEigen(const std::vector<double>& eigen) {
    //   (mEig[0]).resize(eigen.size());
    //   mEig[0] = eigen;
    // };

    // template <typename T1, typename T2> void select(const std::vector<T1>& coord, const std::vector<T2>& v, std::vector<T2>& res) ;
    // void select(const std::vector<int>& coord, std::vector<double>& res) {
    //
    // // int n = coord.size();
    //   if (mIsCube) {
    //     for (int i = 0; i < mDim; ++i)  res[i] = (mEig[0])[coord[i]];
    //   } else {
    //     for (int i = 0; i < mDim; ++i)  res[i] = (mEig[i])[coord[i]];
    //
    //   }
    //
    // } ;


  public:

        dpp_All() { };

        dpp_All(List args) :
        // mRho(as<double>(args["rho"])),
                            mDim(as<int>(args["dim"])),
                            mProg(as<int>(args["progress"])),
                            mProgSim(as<int>(args["simprogress"])),
                            mIsCube(as<bool>(args["ic"])),
                            mWithKernel(as<bool>(args["wk"]))
                            {
                            mWscale = args["Wscale"];
                            mWcenter = args["Wcenter"];

                            // std::cout<<"Size of args :"<<args.size()<<std::endl;

                            mIsEigSet = false;
                          };

        // dpp_All(double rho, int dim) :
        // // mRho(rho),
        //                               mDim(dim) { };

        virtual ~dpp_All() { };

    // NumericVector getEigen(const int j) {
    //   // std::cout << " Call computeEigenVec"
    //   computeEigenVec();
    //   return(NumericVector(mEig[j].begin(), mEig[j].end()));
    // };

    // NumericMatrix getIndex() {
    //   // std::cout << " Call computeEigenVec"
    //   computeIndex();
    //   int s = mIndex.size();
    //   NumericMatrix res (s,mDim);
    //   NumericVector tp (s);
    //   for (int i=0; i < s; ++i) {
    //     tp = NumericVector(mIndex[i].begin(), mIndex[i].end());
    //     res(i,_) = tp;
    //   }
    //   return(res);
    // };


    // NumericVector test_binom(const int k, const int n, const double p) {
    //   // Set seed
    //   RNGScope scope;
    //   return rbinom(k, n, p);
    // };
    // virtual void computeEigen(const std::vector<int> K) = 0;
    // virtual void computeEigenForKernel() = 0;
    virtual void computeIndex() = 0 ;

    // virtual NumericVector computeKernel() = 0;

    virtual ComplexMatrix computeKernelR(const NumericMatrix& PP) = 0;
    // virtual ComplexMatrix computeExactKernelR(const NumericMatrix& PP) = 0;

    List computeOnlyKernelR(const List& PP);
    NumericMatrix computePCF(const NumericMatrix& PP);
    List computePCFR(const List& PP);


    NumericMatrix computeSample();
    List computeListSamples(const int nsim);



};


template <typename T> void print_vector(std::vector<T> v);
void print_vector(NumericVector v);


bool next_variation(std::vector<int>::iterator first, std::vector<int>::iterator last, const int init, const int max) ;
std::complex<double> computeFourierbasis(const std::vector<int>& k, const NumericVector& x, const NumericVector& boxlengths) ;
std::complex<double> computeFourierbasis(const std::vector<int>& k, const NumericVector& x) ;


#endif
