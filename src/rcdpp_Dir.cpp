#include "rcdpp_Dir.h"

using namespace Rcpp;

double dpp_Dir::methodComputationEigen (const std::vector<int>& coord, std::vector<double>& vec) {
  double res;
  if(mDim == 1) {
    select(coord, vec);
    res = vec[0];
  } else {
    if (mType.compare("prod") == 0) {
      select(coord, vec);
      res = std::accumulate(vec.begin(), vec.end(), 1., std::multiplies<double>());    // Eigenvalue associated to coord
    }

    if (mType.compare("sum") == 0) {
      // std::cout<<"Compute eigen value (type = 'sum')"<<std::endl;
      int count_nzeros = 0;
      // int i = 0;
      int ind_nzeros;
      for (int i = 0; (i < mDim) & (count_nzeros < 2); ++i) {
      // while (count_nzeros < 2) {
      // std::cout<<"coord["<<i<<"] = "<< coord[i]<<std::endl;
        if (coord[i] != 0) {
          count_nzeros++;
          ind_nzeros = i;
        }
      }
      if (count_nzeros == 0) {
        // std::cout<<"Que des 0 dans la grille"<<std::endl;
        select(coord, vec);
        res = std::accumulate(vec.begin(), vec.end(), 0.);
      } else if (count_nzeros == 1) {
        // std::cout<<"Une seule non-nulle dans la grille"<<std::endl;
        select(coord, vec);
        res = vec[ind_nzeros];
      } else {
        // std::cout<<"Plus d'une valeur nonnulle dans la grille"<<std::endl;
        res = 0;
      }
    }
  }
  return res;

};

void dpp_Dir::computeEigenDir() {


  // std::cout<<"In computeEigenVec"<<std::endl;
  std::vector<double> eig;
  // double tpInt = 0.;

  int ni;
  int size;
  if (mType.compare("prod") == 0) {
    size = mK;
  } else if (mType.compare("sum") == 0) {
    size = mK+1;
  }
  // int n1=0,n0=0;

  if (mIsCube) {

    std::vector<double> eig(size, 1.);

    if (mType.compare("sum") == 0) eig[0] = 0.;
    setEigenDir(eig, 0);

  } else {
    for (int i = 0; i < mDim; ++i) {
      ni = mN[i];
      if (mType.compare("sum") == 0) ni++;
      eig.clear();
      for (int tpk = 0; tpk < size; tpk++) {
        if (tpk < ni) {
          eig.push_back(1.);
          // n1++;
          // tpInt++;
        } else {
          eig.push_back(0.);
        }
      }
      if (mType.compare("sum") == 0) eig[0] = 0.;
      setEigenDir(eig, i);

      }
    }
    // std::cout<<"computeEigenVec done."<<std::endl;
}


void dpp_Dir::computeIndex() {

    // std::cout<<"In computeIndex: call computeEigenVec"<<std::endl;
    // computeEigenDir();
    // std::cout<<"In computeIndex: computeEigenVec done"<<std::endl;
    int init = 0;
    int end;

    if (mType.compare("prod") == 0) {
      // init = 0;
      end = mK-1;
    } else if (mType.compare("sum") == 0) {
      // init = 1;
      end = mK;
    }

    std::vector<int> coord (mDim, init);              // vector of a possible permutations of {0,...,k-1}^d initialized with 0^d
    std::vector<double> temp(mDim, 0.);                // subvector of K coresponding to coord
    std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d


    double tp;
    int curk = 0;
    // select(coord, mEig, temp);
    // select(coord, temp);
    // if (mDim == 1) tp = temp[0];
    // else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());
    tp = methodComputationEigen(coord, temp);

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
    // std::cout<<"totk = "<<totk<<std::endl;
    // if (!mIsOdd) --totk;
    std::vector<int>::iterator pos = coord.end();

    while (next_variation(coord.begin(), pos, init, end)) {
      ++curk;
    // while (next_variation(coord, totk)) {
      // select(coord, temp);
      tp = methodComputationEigen(coord, temp);
      // if (mDim == 1) tp = temp[0];
      // else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());

      if (mWithKernel & !mIsEigSet) {
        mEig[curk] = tp;
        // if (tp > 0) ++mInt;
      }

      if (tp == 1) {
        res.push_back(coord);
        // std::cout<<"Keep this coord: ";
        // print_vector(coord);
      }
      if (mType.compare("sum") == 0) {
        if (*(pos-1) == end) {
          *(pos-1) = init;
          --pos;
          // print_vector(coord);
        }
      }
    }

    // std::cout<<"computeIndex done: size of res ="<<res.size()<<std::endl;

    setIndex(res);
    if (mWithKernel & !mIsEigSet) mIsEigSet = true;
    // std::cout<<"computeIndex done."<<std::endl;

}


void dpp_Dir::computeEigenForKernel() {

    int init = 0;
    int end;

    if (mType.compare("prod") == 0) {
      // init = 0;
      end = mK-1;
    } else if (mType.compare("sum") == 0) {
      // init = 1;
      end = mK;
    }

    std::vector<int> coord (mDim, init);              // vector of a possible permutations of {0,...,k-1}^d initialized with 0^d
    std::vector<double> temp(mDim, 0.);                // subvector of K coresponding to coord
    std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d

    double tp;
    int curk = 0;
    // // select(coord, mEig, temp);
    // select(coord, temp);
    // if (mDim == 1) tp = temp[0];
    // else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());
    tp = methodComputationEigen(coord, temp);
    mEig[curk] = tp;


    std::vector<int>::iterator pos = coord.end();
    while (next_variation(coord.begin(), pos, init, end)) {
      ++curk;
    // while (next_variation(coord, totk)) {
      // select(coord, temp);
      //
      // if (mDim == 1) tp = temp[0];
      // else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());
      tp = methodComputationEigen(coord, temp);
      // if (mDim == 1) tp = temp[0];
      // else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());

      mEig[curk] = tp;

      if (mType.compare("sum") == 0) {
        if (*(pos-1) == end) {
          *(pos-1) = init;
          --pos;
          // print_vector(coord);
        }
      }

    }

    // std::cout<<"computeIndex done: size of res ="<<res.size()<<std::endl;

    if (!mIsEigSet) mIsEigSet = true;
    // std::cout<<"computeIndex done."<<std::endl;

}



std::complex<double> dpp_Dir::computeKernel(const NumericVector& X, const NumericVector& Y) {

  // int n = mIndextot.size();
  // std::vector<int> K;
  // int n = 2*mK;

  int minK = 0;
  int maxK;
  if (mType.compare("prod") == 0) {
    // minK = 0;
    maxK = mK-1;
  } else if (mType.compare("sum") == 0) {
    // minK = 1;
    maxK = mK;
  }
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
  std::vector<int>::iterator pos = coord.end();
  while (next_variation(coord.begin(), pos, minK, maxK)) {
    ++curK;
    tpEig = mEig[curK];
    res += tpEig*computeFourierbasis(coord, Z, mWscale);
    // res += tpcx;
    if (mType.compare("sum") == 0) {
      if (*(pos-1) == maxK) {
        *(pos-1) = minK;
        --pos;
        // print_vector(coord);
      }
    }
  }
  // res.r = tp.real(); res.i = tp.imag();
  return res;

}
