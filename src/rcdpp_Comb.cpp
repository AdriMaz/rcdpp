#include "rcdpp_Comb.h"

using namespace Rcpp;

double dpp_Comb::methodComputationEigen (const std::vector<int>& coord, std::vector<double>& vec) {
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
        if (coord[i] != mK) {
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
// compute eigen values for each direction
//
void dpp_Comb::computeEigenDir() {
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
      tpeig = mKernel1D->computeEigen1D(tpk, mWscale[0]);
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
          tpeig = mKernel1D->computeEigen1D(tpk, mWscale[i]);
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

void dpp_Comb::computeEigenForKernel() {

      int init;

      if (mType.compare("prod") == 0) {
        init = 0;
      } else if (mType.compare("sum") == 0) {
        init = mK;
      }

      std::vector<int> coord (mDim, init);               // vector of a possible permutations of {-k,...,k}

      int curk = init;


      double tp;

      std::vector<double> temp (mDim, 0.);          // vector of eigenvalues coresponding to coord
      // select(coord, temp);
      // if(mDim == 1) tp = temp[0];
      // else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());    // Eigenvalue associated to coord

      tp = methodComputationEigen(coord, temp);

      // if (mWithKernel & !mIsEigSet) {
      mEig[curk] = tp;
      if (tp > 0) mInt = tp;
      // }
      // std::cout << "Corresponding eigenvalue.: " << tp << std::endl;

      // if (mWithKernel) mEig.push_back(tp);
      // Repeat for all the permutations
      int end = 2*mK;
      std::vector<int>::iterator pos = coord.end();
      bool test = next_variation(coord.begin(), pos, init, end);
      while (test) {
        ++curk;
        // select(coord, temp);
        // if(mDim == 1) tp = temp[0];
        // else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());
        tp = methodComputationEigen(coord, temp);
        mEig[curk] = tp;
        if (tp > 0) mInt += tp;
        if (mType.compare("sum") == 0) {
          if (*(pos-1) == end) {
            *(pos-1) = init;
            --pos;
            // print_vector(coord);
          }
        }
        test = next_variation(coord.begin(), pos, init, end);
      }
      if (mType.compare("sum") == 0) {
        curk = init;
        end = 0;
        pos = coord.end();
        test = next_variation(coord.begin(), pos, init, end);
        while (test) {
          ++curk;
          // select(coord, temp);
          // if(mDim == 1) tp = temp[0];
          // else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());
          tp = methodComputationEigen(coord, temp);
          mEig[curk] = tp;
          if (tp > 0) mInt += tp;
          if (mType.compare("sum") == 0) {
            if (*(pos-1) == end) {
              *(pos-1) = init;
              --pos;
              // print_vector(coord);
            }
          }
          test = next_variation(coord.begin(), pos, init, end);
        }
      }

      if (!mIsEigSet) mIsEigSet = true;

}


void dpp_Comb::computeIndex() {

  // std::cout<<"In computeIndex"<<std::endl;
  // if (mEigDir.size() == 0) computeEigenDir();
    int init;

    if (mType.compare("prod") == 0) {
      init = 0;
    } else if (mType.compare("sum") == 0) {
      init = mK;
    }

    std::vector<int> coord (mDim, init);               // vector of a possible permutations of {-k,...,k}
    std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d
    // std::cout<<"current coord"<<print_vector(coord);
    int curk = init;
    // double tpInt =0. ;
      // int totk = 2*k+1;


    double tp;

    std::vector<double> temp (mDim, 0.);          // vector of eigenvalues coresponding to coord
    // select(coord, temp);
    // if(mDim == 1) tp = temp[0];
    // else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());    // Eigenvalue associated to  coord
    // std::cout<<"Proposed point:  ("<<temp[0]<<" , "<<temp[1]<<")"<<std::endl;
    // std::cout << "Point with coord.: ("<<coord[0]<<" , "<<coord[1]<<")"<<std::endl;
    tp = methodComputationEigen(coord, temp);

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
    int end = 2*mK;
    std::vector<int>::iterator pos = coord.end();
    bool test = next_variation(coord.begin(), pos, init, end);;
    // ++coord[mDim-1];
    while (test) {
      ++curk;
      tp = methodComputationEigen(coord, temp);
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
      if (mType.compare("sum") == 0) {
        if (*(pos-1) == end) {
          *(pos-1) = init;
          --pos;
          // print_vector(coord);
        }
      }
      test = next_variation(coord.begin(), pos, init, end);
    }
    if (mType.compare("sum") == 0) {
      // std::cout<<"2e étape : on part dans l'autre sens"<<std::endl;
      curk = init;
      end = 0;
      pos = coord.end();
      test = next_variation(coord.begin(), pos, init, end);
      while (test) {
        --curk;
        tp = methodComputationEigen(coord, temp);
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
        if (mType.compare("sum") == 0) {
          if (*(pos-1) == end) {
            *(pos-1) = init;
            --pos;
            // print_vector(coord);
          }
        }
        test = next_variation(coord.begin(), pos, init, end);
      }
    }
    // std::cout << "In computeIndex :  Size of the res = " << res.size() << std::endl;
    // print_vector(res);

    // std::cout << "In computeIndex :  Size of the mIndex = " << mIndex.size() << std::endl;


      // }
    setIndex(res);
    if (mWithKernel & !mIsEigSet) mIsEigSet = true;
    // setInt(tpInt);

}


ComplexMatrix dpp_Comb::computeKernelR(const NumericMatrix& PP) {

    if (!mIsEigSet & (mK >= 0)) computeEigenForKernel();

    int np = PP.nrow();
    ComplexMatrix res (np,np);
    NumericVector tpX (mDim);

    Rcomplex diag, tpR;
    // std::cout<<"mInt ="<<mInt<<std::endl;
    diag.r = mInt; diag.i = 0.;
    std::complex<double> tp;

    for (int i = 0; i < np; ++i) {
      // res(i,i) = diag;
      tpX = PP(i,_);
      for (int j = 0; j <= i; ++j) {
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

// ComplexMatrix dpp_Comb::computeExactKernelR(const NumericMatrix& PP) {
//
//     // if (!mIsEigSet & (mK > 0)) computeEigenForKernel();
//
//     int np = PP.nrow();
//     ComplexMatrix res (np,np);
//     NumericVector tpX (mDim);
//
//     Rcomplex diag, tpR;
//     diag.r = mInt; diag.i = 0.;
//     std::complex<double> tp;
//
//     for (int i = 0; i < np; ++i) {
//       res(i,i) = diag;
//       tpX = PP(i,_);
//       for (int j = 0; j < i; ++j) {
//         // if (j != i) {
//           // tpY = PP(j,_);
//           tp = computeExactKernel(tpX, PP(j,_));
//           tpR.r = tp.real(); tpR.i = tp.imag();
//           res(i,j) = tpR;
//           tpR.i *= -1;
//           res(j,i) = tpR;
//         // }
//       }
//     }
//
//
//   return res;
//
// }


// List dpp_Comb::computeOnlyExactKernelR(const List& PP) {
//
//   int np = PP.size();
//   List RES;
//   // if (!mIsEigSet) computeEigenForKernel();
//   for (int i=0; i < np; ++i) {
//     NumericMatrix pp = PP[i];
//     // Rcout<<"nrow(pp) ="<<pp.nrow()<<" ; ncol(pp) ="<<pp.ncol()<<"\n";
//     // for (int j =0; j < pp.nrow(); ++j) {
//     //   NumericVector tp = pp(j,_);
//     //   print_vector(tp);
//     // }
//
//     ComplexMatrix tpres = computeExactKernelR(pp);
//     RES.push_back(tpres);
//   }
//   return RES;
// }


std::complex<double> dpp_Comb::computeKernel(const NumericVector& X, const NumericVector& Y) {

  std::complex<double> res = 0.;
  if (mK >= 0){
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
    res = computeExactKernel(X, Y);
  }
  // res.r = tp.real(); res.i = tp.imag();
  return res;

}
