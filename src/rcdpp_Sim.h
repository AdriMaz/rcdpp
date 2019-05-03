#ifndef RCDPP_SIM_H
#define RCDPP_SIM_H


#define _USE_MATH_DEFINES     // get M_PI

#include <Rcpp.h>
// #include <math.h>
// #include <iostream>
// #include <vector>
// #include <random>
// #include <complex>


using namespace Rcpp;


// Mother class of DPPs
class dpp_All {

  protected:

    // double mRho;    // Intensity
    int mDim;         // Dimension

    // NumericVector mBinfs;
    // NumericVector mBsups;
    NumericVector mWscale;
    NumericVector mWcenter;

    // double mPws;        // prod(mWscale)

    int mProg;      // Progress of the i° sampling
    int mProgSim;   // Progress of the sampling

    bool mIsProj;      // Is it a Projection DPP or not?
    bool mIsProd;      // Is it a d-dimensional DPP with a kernel which is a product of d identical one-deimensional kernels?

    bool mIsCube;    // Is the domain a cube (i.e. product of d identical intervals)?

    std::vector<std::vector<int> > mIndex;   // Indexes in the Mercer's decomposition kept after Bernoulli drawns

    std::vector<std::vector<double> > mEig;                // Eigenvalues (i.e. Fourier Transform computed on [-k ; k])
    // std::vector<double> mEig;

    void setIndex(const std::vector<std::vector<int> >& index) {
      mIndex.resize(index.size());
      mIndex = index;
    };

    void resetIndex() {
      mIndex.clear();
    };

    virtual void setEigen(const std::vector<double>& eigen, const int i) {
      (mEig[i]).resize(eigen.size());
      mEig[i] = eigen;
    };

    // void setEigen(const std::vector<double>& eigen) {
    //   (mEig[0]).resize(eigen.size());
    //   mEig[0] = eigen;
    // };

    // template <typename T1, typename T2> void select(const std::vector<T1>& coord, const std::vector<T2>& v, std::vector<T2>& res) ;
      void select(const std::vector<int>& coord, std::vector<double>& res) {

      // int n = coord.size();
      if (mIsCube) {
        for (int i = 0; i < mDim; ++i)  res[i] = (mEig[0])[coord[i]];
      } else {
        for (int i = 0; i < mDim; ++i)  res[i] = (mEig[i])[coord[i]];

      }

    } ;


  public:

        dpp_All() { };

        dpp_All(List args) :
        // mRho(as<double>(args["rho"])),
                            mDim(as<int>(args["dim"])),
                            mProg(as<int>(args["progress"])),
                            mProgSim(as<int>(args["simprogress"])),
                            mIsCube(as<bool>(args["ic"]))
                            {
                            // mBinfs = args["binfs"];
                            // mBsups = args["bsups"];
                            mWscale = args["Wscale"];
                            mWcenter = args["Wcenter"];

                            if (!mIsCube) mEig.resize(mDim);
                            else mEig.resize(1);
                            // std::cout<<"Wscale = ("<<mWscale[0]<<","<<mWscale[1]<<")"<<std::endl;
                            // std::cout<<"Wcenter = ("<<mWcenter[0]<<","<<mWcenter[1]<<")"<<std::endl;
                            // mPws = std::accumulate(mWscale.begin(), mWscale.end(), 1., std::multiplies<double>());
                            // mEig.resize(mDim);
                            // mIsCube = (min(mBinfs) == max(mBinfs)) && (min(mBsups) == max(mBsups));

                            // if(mIsCube) std::cout<<"Domain is cubic"<<std::endl;

                          };

        dpp_All(double rho, int dim) :
        // mRho(rho),
                                      mDim(dim) { };

        ~dpp_All() { };

    NumericVector getEigen(const int k, const int j) {
      // std::cout << " Call computeEigenVec"
      computeEigenVec(k);
      return(NumericVector(mEig[j].begin(), mEig[j].end()));
    };


    // NumericVector test_binom(const int k, const int n, const double p) {
    //   // Set seed
    //   RNGScope scope;
    //   return rbinom(k, n, p);
    // };
    virtual void computeEigenVec(const int k) = 0 ;
    // virtual void computeEigen(const std::vector<int> K) = 0;
    virtual void computeIndex(const int k) = 0 ;

    NumericMatrix computeSample(const int k);

    // virtual NumericMatrix testSample();

    List computeListSamples(const int k, const int nsim);



};


// Class of d-dimensional DPPs whose each eigenvalue (indexed by Z^d) is product of values indexed by Z
// contains: Gaussian, L1Exponential and MRProd DPPs
class dpp_Prod : public dpp_All {


public:

  dpp_Prod() : dpp_All() {
    mIsProd = true;
  };

  dpp_Prod(List args) : dpp_All(args) {

    mIsProd = true;

    // If the definition set is a ``cube'', margins eigen values are identical on each direction
    // if (!mIsCube) mEig.resize(mDim);
    // else mEig.resize(1);
  };

  dpp_Prod(double rho, int dim) : dpp_All(rho, dim) {
    mIsProd = true;
    // if (!mIsCube) mEig.resize(mDim);
    // else mEig.resize(1);
  };

  ~dpp_Prod() { };


  void computeEigenVec(const int k);

  void computeIndex(const int k);

  virtual double computeEigen(const int k, double wsc) = 0;

};

class dpp_Gauss : public dpp_Prod {

  private:

    double mRho;    // Intensity
    double mAlpha;  // alpha parameter
    // int mDim;         // Dimension

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
                            // std::cout << "Construction of dpp_Gauss" << std::endl;
                            // std::cout << "Rho  = " << mRho << std::endl;
                            // std::cout << "Alpha = " << mAlpha << std::endl;
                            // std::cout << "Dim = " << mDim << std::endl;
                         };


    // dpp_Gauss(double rho, double alpha, int dim) : dpp_Prod(rho, dim),
    //                                               mAlpha(alpha) {
    //                                               // mAsprod = true;
    //                                               mIsProj = false;
    //                                             };

    ~dpp_Gauss() { };

    // void computeEigenVec(const int k);

    double computeEigen(const int k, double wsc) {
      double res;
      res = pow(mRho, 1./mDim)*sqrt(M_PI)*mAlpha*exp(-pow(k/wsc*mAlpha*M_PI,2));
      return res;
    };

};


class dpp_L1Exp : public dpp_Prod {

  private:

    double mRho;    // Intensity
    double mAlpha;  // alpha parameter
    // int mDim;         // Dimension

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

    };


    // dpp_L1Exp(double rho, double alpha, int dim) : dpp_Prod(rho, dim),
    //                                              mAlpha(alpha) {
    //   // mAsprod = true;
    //   mIsProj = false;
    //
    //  };

    ~dpp_L1Exp() { };


    // void computeEigen(const int k);
    double computeEigen(const int k, double wsc) {
      double res;
      res = pow(mRho, 1./mDim)*2*mAlpha/(1.+pow(2*M_PI*mAlpha*k/wsc, 2));
      return res;
    };


};




// Most repulsive stationary DPP

// class dpp_MR : public dpp_All {
//
//   private:
//     double mTau;
//
//   public:
//
//     dpp_MR() : dpp_All() {
//       // mAsprod = false;
//       mIsProj = true;
//       mIsProd = false;
//
//       if (!mIsCube) mEig.resize(mDim);
//       else mEig.resize(1);
//     };
//
//     dpp_MR(List args) : dpp_All(args),
//                         mTau(as<double>(args["tau"])) {
//
//       // std::cout<<"Tau = "<<mTau << std::endl;
//       // mAsprod = false;
//       mIsProj = true;
//       mIsProd = false;
//
//       if (!mIsCube) mEig.resize(mDim);
//       else mEig.resize(1);
//     };
//
//
//     // dpp_MR(double rho, int dim) : dpp_All(rho, dim) {
//     //   if (mDim == 1) mTau = mRho/2.;
//     //   else if (mDim == 2) mTau = mRho/M_PI;
//     //   else mTau = mRho*mDim*tgamma(mDim/2.)/(2.*pow(M_PI, mDim/2.));
//     //
//     //   // mAsprod = false;
//     //   mIsProj = true;
//     //   mIsProd = false;
//     //
//     //   if (!mIsCube) mEig.resize(mDim);
//     //   else mEig.resize(1);
//     // };
//
//     ~dpp_MR() { };
//
//
//     void computeIndex(const int k);
//
//     void computeEigenVec(const int k);
//
//
// };


// class dpp_Dir0 : public dpp_All {
//
//   private:
//     int mN;
//     bool mIsOdd;
//
//   public:
//
//     dpp_Dir0() : dpp_All() {
//       // mAsprod = false;
//       mIsProj = true;
//
//
//     };
//
//     dpp_Dir0(List args) : dpp_All(args),
//                           mN(as<int>(args["n0"])),
//                           mIsOdd(as<bool>(args["odd"])) {
//       // mTau = (pow(mRho, 1./mDim)-1.)/2.;
//       // std::cout<<"Tau = "<<mTau << std::endl;
//       // mAsprod = false;
//       // std::cout<<"mN = "<<mN<<", mIsOdd = "<<mIsOdd<<std::endl;
//       mIsProj = mIsOdd;
//       if (!mIsOdd & mIsCube) {
//         mIsCube = !mIsCube;
//         mEig.resize(mDim);
//       }
//
//     };
//
//
//     ~dpp_Dir0() { };
//
//
//     void computeIndex(const int k);
//
//     void computeEigenVec(const int k);
//
//
// };

class dpp_Dir : public dpp_All {

  private:
    IntegerVector mN;
    // IntegerVector mIsOdd;

  public:

    dpp_Dir() : dpp_All() {
      // mAsprod = false;
      mIsProj = true;
    };

    dpp_Dir(List args) : dpp_All(args) {
      // mTau = (pow(mRho, 1./mDim)-1.)/2.;
      // std::cout<<"Tau = "<<mTau << std::endl;
      // mAsprod = false;
      // std::cout<<"mN = "<<mN<<", mIsOdd = "<<mIsOdd<<std::endl;

      // std::cout<<"Read 'N' arguments"<<std::endl;
      mN = args["N"];
      // std::cout<<"Done"<<std::endl;
      // mIsOdd = args["odd"];

      // mIsCube = false;
      mEig.resize(mDim);

      // mIsProj = (bool)(std::accumulate(mIsOdd.begin(), mIsOdd.end(), 1, std::multiplies<int>()));

      // if (!mIsProj & mIsCube) {
      //   mIsCube = !mIsCube;
      //   mEig.resize(mDim);
      // }

    };


    ~dpp_Dir() { };


    void computeIndex(const int k);

    void computeEigenVec(const int k);


};



//
// class dpp_Eig : public dpp_All {
class dpp_Eig {

  private:

    int mDim;

    NumericVector mWscale;
    NumericVector mWcenter;

    // double mPws;        // prod(mWscale)

    int mProg;      // Progress of the i° sampling
    int mProgSim;   // Progress of the sampling

    bool mIsProj;      // Is it a Projection DPP or not?
    // bool mIsProd;      // Is it a d-dimensional DPP with a kernel which is a product of d identical one-deimensional kernels?

    // bool mIsCube;    // Is the domain a cube (i.e. product of d identical intervals)?

    std::vector<std::vector<int> > mIndex;   // Indexes in the Mercer's decomposition kept after Bernoulli drawns

    // std::vector<double> mEig;
    NumericVector mEig;

    std::vector< std::vector<int> > mIndextot;

    void setIndex(const std::vector<std::vector<int> >& index) {
      mIndex.resize(index.size());
      mIndex = index;
    };

    void resetIndex() {
      mIndex.clear();
    };

  public:

    // dpp_Eig() : dpp_All() { };
    dpp_Eig() { };
    // dpp_Eig(List args) : dpp_All(args)
    dpp_Eig(List args) : mDim(as<int>(args["dim"])),
                        mProg(as<int>(args["progress"])),
                        mProgSim(as<int>(args["simprogress"]))
                        // mIsCube(as<bool>(args["ic"]))
                        {
                          mWscale = args["Wscale"];
                          mWcenter = args["Wcenter"];

                          // mIsProd = false;
                          mIsProj = as<bool>(args["ip"]);
                          List temp;
                          // std::cout<<"Save 'param.eigenval'"<<std::endl;
                          // temp = args["param"];

                          // mEig.resize(1);
                          mEig = args["eigenval"];
                          // for (int i = 0; i < mEig.size(); ++i) {
                          //   std::cout<< (mEig[i]) << std::endl;
                          // }

                          // std::cout<<"Read 'index'"<<std::endl;

                          temp = args["index"];

                          // std::cout<<"Save 'index'"<<std::endl;

                          // std::vector<int> tp;
                          for (int i = 0; i < temp.size(); ++i) {
                            // tp = temp[i];
                            // for (int j = 0; j < mDim; ++j) std::cout<< (int)(tp[j]) << ";";
                            // std::cout<< std::endl;
                            mIndextot.push_back(temp[i]);
                          }

    };

    ~dpp_Eig() { };


    // void computeEigenVec(const int k) {  };


    void computeIndex();

    List computeListSamples(const int nsim);
    NumericMatrix computeSample();


    //
    // virtual void select(const std::vector<int>& coord, std::vector<double>& res) {
    //
    //   // int n = coord.size();
    //   if (mIsCube) {
    //     for (int i = 0; i < mDim; ++i)  res[i] = (mEig[0])[coord[i]];
    //   } else {
    //     for (int i = 0; i < mDim; ++i)  res[i] = (mEig[i])[coord[i]];
    //
    //   }
    //
    // } ;


};


//

// void print_vector(std::vector<std::complex<double> > v);
// void print_vector(std::vector<std::vector<int> > v);
// void print_vector(std::vector<double> v);

template <typename T> void print_vector(std::vector<T> v);


bool next_variation(std::vector<int>::iterator first, std::vector<int>::iterator last,  const int max) ;
std::complex<double> computeFourierbasis(const std::vector<int>& k, const NumericVector& x, const NumericVector& boxlengths) ;
std::complex<double> computeFourierbasis(const std::vector<int>& k, const NumericVector& x) ;
// void select(const std::vector<int>& coord, const std::vector<double>& v, std::vector<double>& res) ;
// void select(const std::vector<int>& coord, const std::vector<int>& v, std::vector<int>& res) ;
// template <typename T1, typename T2> void select(const std::vector<T1>& coord, const std::vector<T2>& v, std::vector<T2>& res) ;


#endif
