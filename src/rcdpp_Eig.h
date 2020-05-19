#ifndef RCDPP_EIG_H
#define RCDPP_EIG_H

#include "rcdpp_All.h"

class dpp_Eig : public dpp_All {

  private:
    // int mDim;
    // double mInt;

    // NumericVector mWscale;
    // NumericVector mWcenter;

    // double mPws;        // prod(mWscale)

    // int mProg;      // Progress of the i° sampling
    // int mProgSim;   // Progress of the sampling

    // bool mIsProj;      // Is it a Projection DPP or not?
    // bool mIsProd;      // Is it a d-dimensional DPP with a kernel which is a product of d identical one-deimensional kernels?

    // bool mIsCube;    // Is the domain a cube (i.e. product of d identical intervals)?
    // bool mWithKernel;
    // std::vector<std::vector<int> > mIndex;   // Indexes in the Mercer's decomposition kept after Bernoulli drawns

    // std::vector<double> mEig;
    // std::vector<double> mEig;

    std::vector<std::vector<int> > mIndextot;

    // void setIndex(const std::vector<std::vector<int> >& index) {
    //   mIndex.resize(index.size());
    //   mIndex = index;
    // };
    //
    //
    // void resetIndex() {
    //   mIndex.clear();
    // };

  public:

    dpp_Eig() : dpp_All() { };
    // dpp_Eig() { };
    // dpp_Eig(List args) : dpp_All(args)
    dpp_Eig(List args) : dpp_All(args)
                        {
                          mIsEigSet = true;
                          // std::cout<<"Read 'Wscale'"<<std::endl;
                          // mWscale = args["Wscale"];
                          // std::cout<<"'Wscal' OK"<<std::endl;
                          // std::cout<<"Read 'Wcenter'"<<std::endl;
                          // mWcenter = args["Wcenter"];
                          // std::cout<<"'Wcenter' OK"<<std::endl;
                          // mIsProd = false;
                          // std::cout<<"Read 'ip'"<<std::endl;
                          // std::cout<<"'ip' OK"<<std::endl;
                          // std::cout<<"Read 'eigenval'"<<std::endl;
                          NumericVector tpEig = args["eigenval"];
                          // std::cout<<"'eigenval' OK"<<std::endl;
                          // std::cout<<"Read 'index'"<<std::endl;

                          int n = tpEig.size();
                          // mInt = std::accumulate(tpE.begin(), mEig.end(), 0.);
                          double tpi;
                          NumericMatrix temp = args["index"];
                          // std::cout<<"'index' OK"<<std::endl;
                          // std::cout<<"Construct 'mIndextot'"<<std::endl;
                          std::vector<int> tp (mDim);
                          for (int i = 0; i < n; ++i) {
                            tpi = tpEig[i];
                            mEig.push_back(tpi);
                            mInt += tpi;
                            for (int j = 0; j < mDim; ++j) tp[j] = (int)(temp(i,j));
                            mIndextot.push_back(tp);
                          }
                          // std::cout<<"'mIndextot' OK"<<std::endl;

    };

    // ~dpp_Eig() { };


    // void computeEigenVec(const int k) {  };


    void computeIndex();

    ComplexMatrix computeKernelR(const NumericMatrix& PP);
    // NumericMatrix computePCFR(const NumericMatrix& PP);

    std::complex<double> computeKernel(const NumericVector& X, const NumericVector& Y);
    // Rcomplex computeKernel(const NumericVector& X, const NumericVector& Y);
    // std::complex<double> computeKernel(const NumericVector& X, const NumericVector& Y);
    // arma::cx_mat computeKernel(const NumericMatrix& PP);

    // double computekInt(const NumericMatrix& PP);
    //
    // arma::mat computePCF(const NumericMatrix& PP);
    // List computeListSamples(const int nsim);
    // NumericMatrix computeSample();

};


#endif
