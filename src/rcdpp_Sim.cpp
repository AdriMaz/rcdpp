#include "rcdpp_Sim.h"
// void print_vector(std::vector<std::complex<double> > v) {
//   int n = v.size();
//   for (int i = 0; i < n; ++i) std::cout << std::real(v[i]) << " + " << std::imag(v[i]) << "i" << "  ";
//   std::cout << std::endl;
// }
//
// void print_vector(std::vector<std::vector<int> > v) {
//   int n = v.size();
//   std::vector<int> K;
//   for (int i = 0; i < n; ++i) {
//     K = v[i];
//     for (int j = 0; j < K.size() ; ++j) {
//       std::cout << K[j] << "  ";
//     }
//     std::cout << std::endl;
//   }
//   std::cout << std::endl;
// }


template <typename T>
void print_vector(std::vector<T> v) {
  int n = v.size();
  for (int i = 0; i < n; ++i) std::cout << v[i] << " ";
  std::cout << std::endl;
}


bool next_variation(std::vector<int>::iterator first, std::vector<int>::iterator last,  const int max) {

    if (first == last) return false; // empty sequence (n==0)

    std::vector<int>::iterator i(last); --i;   // Point to the rightmost element

    // Check if I can just increase it
    if (*i < max) { ++(*i); return true; } // Increase this element and return

    // Find the rightmost element to increase
    while( i != first )
       {
        *i = 0; // reset the right-hand element
        --i; // point to the left adjacent
        if (*i < max) { ++(*i); return true; } // Increase this element and return
       }

    // If here all elements are the maximum symbol (max=k-1), so there are no more variations
    //for (i=first; i!=last; ++i) *i = 0; // Should reset to the lowest sequence (0)?
    return false;
}

std::complex<double> computeFourierbasis(const std::vector<int>& k, const NumericVector& x, const NumericVector& boxlengths) {

  std::complex<double> dpii (0.,2.*M_PI);
  std::complex<double> res;

  double scal = 0;
  double bli, bl = 1;
  int i;
  int s = x.size();
  for (i = 0; i < s; ++i) {
    bli = boxlengths[i];
    scal += k[i]*x[i]/bli;
    bl *=bli;
  }

  res = exp(dpii*scal)*pow(bl, -.5);

  return res;

}

std::complex<double> computeFourierbasis(const std::vector<int>& k, const NumericVector& x) {

  std::complex<double> dpii (0.,2.*M_PI);
  std::complex<double> res;

  double scal=0;
  int i;
  int s = x.size();
  for (i = 0; i < s; ++i) scal += k[i]*x[i];

  res = exp(dpii*scal);

  return res;

}




List dpp_All::computeListSamples(const int k, const int nsim) {

  List RES;
  NumericMatrix res;
  NumericVector x;
  int i,j,l;

  // Set seed
  RNGScope scope;

  // Rescaling index if necessary

  if (mIsProj) {    // If projection DPP: compute index only once

    // std::cout<<"Call computeIndex (Proj. case)"<<std::endl;
    this->computeIndex(k);
    // std::cout<<"computeIndex done"<<std::endl;

    for (i = 0; i < nsim; ++i) {
      // std::cout << "Sim n° " << i+1 << std::endl;
      if(mProgSim > 0) {
        if ((i+1)%mProgSim == 0) std::cout << "Compute sample n° " << i+1 << std::endl;
      }

      res = this->computeSample(k);   // Compute sample

      for (l = 0; l < res.nrow(); ++l){
        x = res(l,_);
        for(j = 0; j < mDim; ++j) {
          // std::cout<<"x[j] = "<<x[j]<<std::endl;
          // std::cout<<"WSC[j] = "<<mWscale[j]<<std::endl;
          // std::cout<<"WC[j] = "<<mWcenter[j]<<std::endl;
          x[j] = x[j]*mWscale[j]+mWcenter[j];
          // std::cout<<"Affine transf."<<std::endl;
          // std::cout<<"x[j] = "<<x[j]<<std::endl;
        }
        res(l,_) = x;
      }

      RES.push_back(res);       // Stock new point process
    }
  } else {
  // std::cout<< "Compute Eigen" << std::endl;
    this->computeEigenVec(k);
    // std::cout<< "Done -> mEig.size = " << mEig.size() << std::endl;
    // print_vector(mEig);

    for (i = 0; i < nsim; ++i) {

      this->computeIndex(k);

      if(mProgSim > 0) {
        if ((i+1)%mProgSim == 0) std::cout << "Compute sample n° " << i+1 << std::endl;
      }
      // std::cout << "Sim n° " << i+1 << std::endl;
      res = this->computeSample(k);   // Compute sample

      for (l = 0; l < res.nrow(); ++l){
        x = res(l,_);
        for(j = 0; j < mDim; ++j) {
          x[j] = x[j]*mWscale[j]+mWcenter[j];
        }
        res(l,_) = x;
      }

      RES.push_back(res);       // Stock new point process
      // std::cout << "Before reset  mIndex.size() = " << mIndex.size() << std::endl;
      this->resetIndex();       // reset mIndex
      // std::cout << "Reset Index -> mIndex.size() = " << mIndex.size() << std::endl;
    }
  }
  return RES;
}


NumericMatrix dpp_All::computeSample(const int k) {

  // Remark: Sample is computed on [-1/2 ; 1/2]^d
  // to be rescaled on domain W

  int rejectmax = 1e4;          // Max number of repetition for reject algorithm
  // std::cout<< "Compute Index" << std::endl;
  // if (!mIsProj) this->computeIndex(k);
  int n = mIndex.size();
  // std::cout<< "Number of points to be computed " << n << std::endl;
  // print_vector(mIndex);
  int i, j, it;

  if (n == 0) {
    NumericMatrix res (1, 1);
    res(0, 0) = 0;
    return res;
  }
  // std::vector<std::vector<double> > res;
  NumericMatrix res (n, mDim);    // matrix of coordinates

  // First point: uniformaly distributed
  NumericVector x;

  // if (mIsCube) {
  //   x = runif(mDim, mBinfs[0], mBsups[0]);
  // } else {
  //   for (i = 0; i < mDim; ++i) x.push_back(runif(1, mBinfs[i], mBsups[i])[0]);
  // }
  x = runif(mDim, -0.5, 0.5);
  // res.push_back(x);
  // std::cout<< "Adding first point" << std::endl;
  // for(i = 0; i < mDim; ++i) {
  //   x[i] = x[i]*mWscale[i]+mWcenter[i];
  // }
  res(0,_) = x;

  // std::cout<<"Done"<<std::endl;
  if (n == 1) return res;

  std::vector<int> K(mDim);                   // Element of vector index

  std::vector<std::complex<double> > v(n);      // Vector of eigenvector evaluated at a given point
  std::vector<std::complex<double> > w(n);      // Temp vector for each computation
  std::vector<std::complex<double> > wei;            // Temp vector of weight;
  std::vector<std::complex<double> > etp(n);    // Temp vector e_i
  std::vector<std::complex<double> > etpstar(n);   // Conj of e_i
  std::vector<std::vector<std::complex<double> > > e;    // All vectors e_i
  std::vector<std::vector<std::complex<double> > > estar;   // All conj of vectors e_i

  std::complex<double> tpcplx;

  double tp = 0, scal = 0, nv = 0;

  double accept;    // Accept probability

  // std::cout<<"Basis-functions evaluated at first point"<<std::endl;
  for (i = 0; i < n; ++i) {       // Basis-functions evaluated at first point
    K = mIndex[i];
    tpcplx = computeFourierbasis(K, x);
    v[i] = tpcplx;
    nv += pow(abs(tpcplx), 2.);
  }

  nv = pow(nv, 0.5);
  // std::cout<<"Record normalized version in the Gram-Schmidt matrices"<<std::endl;

  for (int i = 0; i < n; ++i ) {       // Record normalized version in the Gram-Schmidt matrices:
    etp[i] = v[i]/nv;
    etpstar[i] = std::conj(etp[i]);
  }
  e.push_back(etp); estar.push_back(etpstar);

  int tries;
  // std::cout<< "Main for loop over number of points"<< std::endl;
  for (it = n-1; it > 0; --it) {    // Main for loop over number of points:

    if (mProg > 0) {
      if ((n-it+1)%mProg == 0) std::cout << "Compute point n° " << n-it+1 << std::endl;
    }
    tries = 1;
    bool stop = false;
    wei.resize(n-it);  // Weight vector at it-th step

    // std::cout << "Acceptance algo" << std::endl;
    while(tries < rejectmax && !stop) {

      // std::cout << "Try n° " << tries << std::endl;

      // Proposed point
      // x = runif (mDim, 0., 1.);
      // if (mIsCube) {
      //   x = runif(mDim, mBinfs[0], mBsups[0]);
      // } else {
      //   for (i = 0; i < mDim; ++i) x.push_back(runif(1, mBinfs[i], mBsups[i])[0]);
      // }
      x = runif(mDim, -0.5, 0.5);
      // for (i = 0; i < mDim; ++i)  x[i] = dU(gen);

      nv = 0;
      // Basis functions eval. at proposed point:
      // std::cout<<"Basis functions eval. at proposed point"<<std::endl;
      for (i = 0; i < n; ++i) {
        K = mIndex[i];
        tpcplx = computeFourierbasis(K, x);
        v[i] = tpcplx;
      }



      scal = 0;
      // Vector of projection weights (has length n-it)
      //       wei <- t(v)%*%estar

      // std::cout << "Compute vector of projection weights" << std::endl;
      // itwei = wei.begin();
      for (i = 0; i < n-it; ++i) {
        tpcplx = 0;
        etpstar = estar[i];
        for (j = 0; j < n; ++j) {
          tpcplx += v[j]*etpstar[j];
        }
        wei[i] = tpcplx;
        nv += pow(abs(tpcplx), 2.);
      }

      // accept = 1-mPws*nv/n;     // Accept probability
      accept = 1-nv/n;


      tp = runif (1, 0., 1.)[0];
      if (tp < accept) stop = true;

      ++tries;
      if (tries > rejectmax) std::cerr << "Rejection sampling failed reject_max =" << rejectmax << "times in a row" << std::endl;

    }

    // res.push_back(x);      // Record the accepted point
    // std::cout << "Record the accepted point" << std::endl;

    // std::cout << "Initialy: x = ("<< x[0] << ", "<< x[1]<<")"<<std::endl;
    // for(i = 0; i < mDim; ++i) {
    //   x[i] = x[i]*mWscale[i]+mWcenter[i];
    // }
    // std::cout << "After rescloing: x = ("<< x[0] << ", "<< x[1]<<")"<<std::endl;
    res(n-it,_) = x;

    if (it > 1) {      // while it is not the last point

      // Calculate orthogonal vector for Gram-Schmidt procedure:
      //  w <- v - rowSums(matrix(wei,n,n-i,byrow=TRUE)*e[,1:(n-i)])
      double nw = 0;
      // std::cout<<"Calculate orthogonal vector for Gram-Schmidt procedure 1"<<std::endl;
      for (i = 0; i < n; ++i) {
        tpcplx = 0;
        j = 0;
        for (j = 0; j < n-it; ++j) {
          tpcplx += wei[j]*e[j][i];
        }
        tpcplx = v[i]-tpcplx;
        w[i] = tpcplx;
        nw += pow(abs(tpcplx), 2);
      }

      nw = pow(nw, 0.5);
      // std::cout<<"Calculate orthogonal vector for Gram-Schmidt procedure 2"<<std::endl;
      for (i = 0; i < n ; ++i) {
        tpcplx = w[i]/nw;
        etp[i] = tpcplx;
        etpstar[i] = std::conj(tpcplx);
      }

      e.push_back(etp); estar.push_back(etpstar);
    }

  }

  return res;

}

void dpp_Prod::computeEigenVec(const int k) {
  std::vector<double> eig;

  double wsc;
  // double sumEig;
  if (mIsCube) {
    wsc = mWscale[0];
    // std::cout<<"Wscale = " << wsc << std::endl;
    for (int tpk = -k; tpk < k+1; tpk++)  eig.push_back(computeEigen(tpk, wsc));
    this->setEigen(eig, 0);
    // sumEig = std::accumulate(eig.begin(), eig.end(), 0.);
    // std::cout<<"sum(eig) = " << sumEig << std::endl;
  } else {
    double wscm = 0.;
    // double wsc = 1.;
    for (int i = 0; i < mDim; ++i) {
      wsc = mWscale[i];
      if (wsc != wscm) {   // Compute only if current wsc is different from previous one
        eig.clear();
        for (int tpk = -k; tpk < k+1; tpk++)  eig.push_back(computeEigen(tpk, wsc));
        wscm = wsc;
      }
      this->setEigen(eig, i);
      // this->setEigen(eig);
    }
  }
}


void dpp_Prod::computeIndex(const int k) {

  std::vector<int> coord (mDim, 0);               // vector of a possible permutations of {-k,...,k}
  std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d

    // int totk = 2*k+1;


    double tp;

    // if (mAsprod) {    // Until now: for Gaussian and L1Exponential kernels

      // int totk = 2*k+1;
      // std::vector<double> eig = computeEigen(k);    // Compute eigenvalues on [-k:k]^d
      // std::cout << "In computeIndex: mEig.size() = " << mEig.size() << std::endl;

      // if (mEig.size() == 0) {
        // std::cout << "Computation of eigenvalues required" << std::endl;
        // this->computeEigenVec(k);
      // }

      // else std::cout << "Computation of eigenvalues not required anymore" << std::endl;
      // std::cout << "In computeIndex2: mEig.size() = " << mEig.size() << std::endl;
      // std::vector<int> coord (mDim, 0);               // vector of a possible permutations of {-k,...,k}
      std::vector<double> temp (mDim, 0.);          // vector of eigenvalues coresponding to coord
      // std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d
      // std::cout << "Lambda_K of" << ' ';
      // print_vector(coord);
      // select(coord, mEig, temp);
      this->select(coord, temp);
      // std::cout<<"Proposed point:  ("<<temp[0]<<" , "<<temp[1]<<")"<<std::endl;
      // std::cout << "With coord.: ("<<coord[0]<<" , "<<coord[1]<<")"<<std::endl;
      // int cpt = 1;
      if(mDim == 1) tp = temp[0];
      else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());    // Eigenvalue associated to  coord
      // std::cout << "Corresponding eigenvalue.: " << tp << std::endl;
      double tpbool;

      tpbool = rbinom(1, 1, tp)[0];  // Bernoulli draw
      if (tpbool == 1) res.push_back(coord);     // Keep the elements with proba tp
      // std::cout<< "Result of the draw: " << tpbool <<std::endl;

      // Repeat for all the permutations
      int totk = 2*k;
      while (next_variation(coord.begin(), coord.end(), totk)) {
        // std::cout << "Next permutation" <<std::endl;
        // ++cpt;
        // select(coord, mEig, temp);
        this->select(coord, temp);
        // std::cout<<"Proposed point:  ("<<temp[0]<<" , "<<temp[1]<<")"<<std::endl;
        // std::cout << "With coord.: ("<<coord[0]<<" , "<<coord[1]<<")"<<std::endl;
        if(mDim == 1) tp = temp[0];
        else tp = std::accumulate(temp.begin(), temp.end(), 1., std::multiplies<double>());
        // std::cout << "Corresponding eigenvalue.: " << tp << std::endl;
        // tp = rbinom(1, 1, tp)[0];
        // if (tp == 1.) res.push_back(coord);
        if (tp > 1e-8) {
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
  this->setIndex(res);

}


// Simulate a projection kernel
// std::vector< std::vector<double> > dpp_Gauss::computeSample(std::vector<std::vector<int> > index) {


// void dpp_MR::computeEigenVec(const int k) {
//
//
//   // std::cout<<"In ComputeEigen"<<std::endl;
//   std::vector<double> eig;
//
//   double wsc;
//   if (mIsCube) {
//     wsc = mWscale[0];
//     if (mDim == 1) {
//         for (int tpk = -k; tpk < k+1; tpk++)  eig.push_back(fabs(tpk/wsc));
//     } else {
//       for (int tpk = -k; tpk < k+1; tpk++)  eig.push_back(pow(tpk/wsc, 2));
//     }
//     this->setEigen(eig, 0);
//   } else {
//     double wscm = 0.;
//     // double wsc = 1.;
//     for (int i = 0; i < mDim; ++i) {
//       wsc = mWscale[i];
//       if (wsc != wscm) {   // Compute only if current wsc is different from previous one
//         eig.clear();
//         if (mDim == 1) {
//             for (int tpk = -k; tpk < k+1; tpk++)  eig.push_back(fabs(tpk/wsc));
//         } else {
//           for (int tpk = -k; tpk < k+1; tpk++)  eig.push_back(pow(tpk/wsc, 2));
//         }
//         wscm = wsc;
//       }
//       this->setEigen(eig, i);
//       // this->setEigen(eig);
//     }
//   }
//   // std::cout<<"ComputeEigen Done"<<std::endl;
//
// }



// void dpp_MR::computeIndex(const int k) {
//
//     // int totk = 2*k+1;
//     // int tpk = -k;
//     // std::cout<<"In ComputeIndex"<<std::endl;
//
//     // if (mEig.size() == 0) {
//       // std::cout << "Computation of eigenvalues required" << std::endl;
//     this->computeEigenVec(k);
//     // }
//     std::vector<int> coord (mDim, 0);              // vector of a possible permutations of {-k,...,k} initialized with -k
//     std::vector<double> temp(mDim, 0.);                // subvector of K coresponding to coord
//     std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d
//
//     // select(coord, mEig, temp);
//     this->select(coord, temp);
//     // std::cout<<"Proposed point:  ("<<temp[0]<<" , "<<temp[1]<<")"<<std::endl;
//     // std::cout << "With coord.: ("<<coord[0]<<" , "<<coord[1]<<")"<<std::endl;
//     double tp = 0;
//     // for (i = 0; i < mDim; ++i)  tp += pow(temp[i], 2.);
//     if (mDim == 1) tp = temp[0];
//     else tp = std::accumulate(temp.begin(), temp.end(), 0.);
//
//     if (tp <= mTau) {
//       // std::cout<<"Accepted point:  ("<<coord[0]<<" , "<<coord[1]<<")"<<std::endl;
//       res.push_back(coord);
//     }
//
//     int totk = 2*k;
//     while (next_variation(coord.begin(), coord.end(), totk)) {
//
//       // select(coord, mEig, temp);
//       this->select(coord, temp);
//       // std::cout<<"Proposed point:  ("<<temp[0]<<" , "<<temp[1]<<")"<<std::endl;
//       // std::cout << "With coord.: ("<<coord[0]<<" , "<<coord[1]<<")"<<std::endl;
//       // tp = 0;
//       // for (i = 0; i < mDim; ++i)  tp += pow(temp[i], 2.);
//       if(mDim == 1) tp = temp[0];
//       else tp = std::accumulate(temp.begin(), temp.end(), 0.);
//
//       if (tp <= mTau) {
//
//         res.push_back(coord);
//       }
//     }
//
//     // std::cout<< "Size of res ="<< res.size();
//     this->setIndex(res);
//     // std::cout<<"ComputeIndex Done"<<std::endl;
//
// }


// void dpp_Dir0::computeEigenVec(const int k) {
//
//   std::vector<double> eig;
//
//   if (mIsOdd) {
//     // std::cout<<"Case N = k^d with k odd number."<<std::endl;
//     if (mIsCube) {
//       // std::cout<<"Window is a cube ."<<std::endl;
//
//       for (int tpk = -k; tpk < k+1; tpk++)  eig.push_back(1.);
//       // std::cout<<"Store the vector eig (eig.size() = "<<eig.size()<<")"<<std::endl;
//
//       this->setEigen(eig, 0);
//       // std::cout<<"Done"<<std::endl;
//     } else {
//       for (int i = 0; i < mDim; ++i) {
//         eig.clear();
//         for (int tpk = -k; tpk < k+1; tpk++)  eig.push_back(1.);
//         this->setEigen(eig, i);
//       }
//     }
//   }
//   else {
//     double tpbool;
//     // if (mIsCube) {
//     //   tpbool = rbinom(1, 1, 0.5)[0]; // Bernoulli draw
//     //   if (tpbool == 0) {
//     //     for (int tpk = -k; tpk < k; tpk++)  eig.push_back(1.);
//     //     eig.push_back(0.);
//     //   } else {
//     //     eig.push_back(0.);
//     //     for (int tpk = -k+1; tpk < k+1; tpk++)  eig.push_back(1.);
//     //   }
//     //   this->setEigen(eig, 0);
//     //
//     // } else {
//       for (int i = 0; i < mDim; ++i) {
//         eig.clear();
//         tpbool = rbinom(1, 1, 0.5)[0]; // Bernoulli draw
//         if (tpbool == 0) {
//           for (int tpk = -k; tpk < k; tpk++)  eig.push_back(1.);
//           eig.push_back(0.);
//         } else {
//           eig.push_back(0.);
//           for (int tpk = -k+1; tpk < k+1; tpk++)  eig.push_back(1.);
//         }
//         this->setEigen(eig, i);
//       }
//     // }
//   }
// }
//
// void dpp_Dir0::computeIndex(const int k) {
//
//     // std::cout<<"In computeIndex: call computeEigenVec"<<std::endl;
//     this->computeEigenVec(k);
//     // std::cout<<"In computeIndex: computeEigenVec done"<<std::endl;
//
//     std::vector<int> coord (mDim, 0);              // vector of a possible permutations of {-k,...,k} initialized with -k
//     std::vector<double> temp(mDim, 0.);                // subvector of K coresponding to coord
//     std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d
//
//     double tp = 1;
//     // select(coord, mEig, temp);
//     this->select(coord, temp);
//     if (mDim == 1) tp = temp[0];
//     else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());
//
//     if (tp == 1) res.push_back(coord);
//     // }
//
//     int totk = 2*k;
//
//     // if (!mIsOdd) --totk;
//
//     while (next_variation(coord.begin(), coord.end(), totk)) {
//       this->select(coord, temp);
//
//       if (mDim == 1) tp = temp[0];
//       else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());
//
//       if (tp == 1) res.push_back(coord);
//     }
//
//     this->setIndex(res);
//
// }



void dpp_Dir::computeEigenVec(const int k) {


  // std::cout<<"In computeEigenVec"<<std::endl;
  std::vector<double> eig;

  int ni;
  int n1=0,n0=0;

  if (mIsCube) {
    // std::cout<<"Window is a cube ."<<std::endl;
    ni = mN[0];
    for (int tpk = 0; tpk < k; tpk++)  {
      if (tpk < ni) {
        eig.push_back(1.);
        n1++;
      } else {
        eig.push_back(0.);
        n0++;
      }
    }
    // std::cout<<"Store the vector eig (eig.size() = "<<eig.size()<<")"<<std::endl;
    // std::cout<<"eig.size() = "<<eig.size()<<std::endl;
    // std::cout<<"# of 1 = "<<n1<<std::endl;
    // std::cout<<"# of 0 = "<<n0<<std::endl;
    this->setEigen(eig, 0);
    // std::cout<<"Done"<<std::endl;
  } else {
    for (int i = 0; i < mDim; ++i) {
      ni = mN[i];
      eig.clear();
      for (int tpk = 0; tpk < k; tpk++) {
        if (tpk <= ni) {
          eig.push_back(1.);
          n1++;
        } else {
          eig.push_back(0.);
          n0++;
        }
      }
      // std::cout<<"eig.size() = "<<eig.size()<<std::endl;
      this->setEigen(eig, i);

      }
    }
}


void dpp_Dir::computeIndex(const int k) {

    // std::cout<<"In computeIndex: call computeEigenVec"<<std::endl;
    this->computeEigenVec(k);
    // std::cout<<"In computeIndex: computeEigenVec done"<<std::endl;

    std::vector<int> coord (mDim, 1);              // vector of a possible permutations of {-k,...,k} initialized with -k
    std::vector<double> temp(mDim, 0.);                // subvector of K coresponding to coord
    std::vector< std::vector<int> > res;          // Vector of kept elements in {-k,...,k}^d

    double tp = 1;
    // select(coord, mEig, temp);
    this->select(coord, temp);
    if (mDim == 1) tp = temp[0];
    else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());

    if (tp == 1) {
      res.push_back(coord);
      // std::cout<<"Keep this coord: (";
      // for (int i = 0; i < mDim; i++) std::cout<<coord[i]<<",";
      // std::cout<<")"<<std::endl;
    }
    // }

    int totk = k;

    // if (!mIsOdd) --totk;

    while (next_variation(coord.begin(), coord.end(), totk)) {
      this->select(coord, temp);

      if (mDim == 1) tp = temp[0];
      else tp = std::accumulate(temp.begin(), temp.end(), 1, std::multiplies<double>());

      if (tp == 1) {
        res.push_back(coord);
        // std::cout<<"Keep this coord: (";
        // for (int i = 0; i < mDim; i++) std::cout<<coord[i]<<",";
        // std::cout<<")"<<std::endl;
      }
    }

    this->setIndex(res);

}





////// dpp_Eig class

List dpp_Eig::computeListSamples(const int nsim) {

  List RES;
  NumericMatrix res;
  NumericVector x;
  int i,j,l;

  // Set seed
  RNGScope scope;

  // Rescaling index if necessary

  if (mIsProj) {    // If projection DPP: compute index only once

    this->computeIndex();

    for (i = 0; i < nsim; ++i) {
      // std::cout << "Sim n° " << i+1 << std::endl;
      if(mProgSim > 0) {
        if ((i+1)%mProgSim == 0) std::cout << "Compute sample n° " << i+1 << std::endl;
      }
      res = this->computeSample();   // Compute sample

      for (l = 0; l < res.nrow(); ++l){
        x = res(l,_);
        for(j = 0; j < mDim; ++j) {
          x[j] = x[j]*mWscale[j]+mWcenter[j];
        }
        res(l,_) = x;
      }

      RES.push_back(res);       // Stock new point process
    }
  } else {
  // std::cout<< "Compute Eigen" << std::endl;
    // this->computeEigenVec(k);
    // std::cout<< "Done -> mEig.size = " << mEig.size() << std::endl;
    // print_vector(mEig);

    for (i = 0; i < nsim; ++i) {

      this->computeIndex();

      if(mProgSim > 0) {
        if ((i+1)%mProgSim == 0) std::cout << "Compute sample n° " << i+1 << std::endl;
      }
      // std::cout << "Sim n° " << i+1 << std::endl;
      res = this->computeSample();   // Compute sample

      for (l = 0; l < res.nrow(); ++l){
        x = res(l,_);
        for(j = 0; j < mDim; ++j) {
          x[j] = x[j]*mWscale[j]+mWcenter[j];
        }
        res(l,_) = x;
      }

      RES.push_back(res);       // Stock new point process
      // std::cout << "Before reset  mIndex.size() = " << mIndex.size() << std::endl;
      this->resetIndex();       // reset mIndex
      // std::cout << "Reset Index -> mIndex.size() = " << mIndex.size() << std::endl;
    }
  }
  return RES;
}


NumericMatrix dpp_Eig::computeSample() {


  int rejectmax = 1e4;          // Max number of repetition for reject algorithm
  // std::cout<< "Compute Index" << std::endl;
  // if (!mIsProj) this->computeIndex();
  int n = mIndex.size();
  // std::cout<< "Number of points to be computed " << n << std::endl;
  // print_vector(mIndex);
  int i, j, it;

  if (n == 0) {
    NumericMatrix res (1, 1);
    res(0, 0) = 0;
    return res;
  }
  // std::vector<std::vector<double> > res;
  NumericMatrix res (n, mDim);    // matrix of coordinates

  // First point: uniformaly distributed
  NumericVector x (mDim);

  // if (mIsCube) {
  //   x = runif(mDim, mBinfs[0], mBsups[0]);
  // } else {
  //   for (i = 0; i < mDim; ++i) x.push_back(runif(1, mBinfs[i], mBsups[i])[0]);
  // }

  x = runif(mDim, -0.5, 0.5);

  // x.resize(mDim);
  // x(0) = 0.5173693;
  // res.push_back(x);
  // std::cout<< "Adding first point" << std::endl;
  // for(i = 0; i < mDim; ++i) {
  //   x[i] = x[i]*mWscale[i]+mWcenter[i];
  // }
  res(0,_) = x;

  // std::cout<<"Done"<<std::endl;
  if (n == 1) return res;

  std::vector<int> K(mDim);                   // Element of vector index

  std::vector<std::complex<double> > v(n);      // Vector of eigenvector evaluated at a given point
  std::vector<std::complex<double> > w(n);      // Temp vector for each computation
  std::vector<std::complex<double> > wei;            // Temp vector of weight;
  std::vector<std::complex<double> > etp(n);    // Temp vector e_i
  std::vector<std::complex<double> > etpstar(n);   // Conj of e_i
  std::vector<std::vector<std::complex<double> > > e;    // All vectors e_i
  std::vector<std::vector<std::complex<double> > > estar;   // All conj of vectors e_i

  std::complex<double> tpcplx;

  double tp = 0, scal = 0, nv = 0;

  double accept;    // Accept probability

  // std::cout<<"Basis-functions evaluated at first point"<<std::endl;
  for (i = 0; i < n; ++i) {       // Basis-functions evaluated at first point
    K = mIndex[i];
    tpcplx = computeFourierbasis(K, x);
    v[i] = tpcplx;
    nv += pow(abs(tpcplx), 2.);
  }

  // std::cout<<"v =";
  // print_vector(v);
  // std::cout<<std::endl;

  nv = pow(nv, 0.5);
  // std::cout<<"nv ="<<nv<<std::endl;
  // std::cout<<"Record normalized version in the Gram-Schmidt matrices"<<std::endl;

  for (int i = 0; i < n; ++i ) {       // Record normalized version in the Gram-Schmidt matrices:
    etp[i] = 1/nv*v[i];
    // etp[i] = std::complex<double>(v[i].real()/nv, v[i].imag()/nv);
    etpstar[i] = std::conj(etp[i]);
  }
  //
  // std::cout<<"etp =";
  // print_vector(etp);
  // std::cout<<std::endl;

  e.push_back(etp); estar.push_back(etpstar);

  int tries;
  // std::cout<< "Main for loop over number of points"<< std::endl;
  for (it = n-1; it > 0; --it) {    // Main for loop over number of points:

    if (mProg > 0) {
      if ((n-it+1)%mProg == 0) std::cout << "Compute point n° " << n-it+1 << std::endl;
    }
    tries = 1;
    bool stop = false;
    wei.resize(n-it);  // Weight vector at it-th step

    // std::cout << "Acceptance algo" << std::endl;
    while(tries < rejectmax && !stop) {

      // std::cout << "Try n° " << tries << std::endl;

      // Proposed point
      // x = runif (mDim, 0., 1.);
      // if (mIsCube) {
      //   x = runif(mDim, mBinfs[0], mBsups[0]);
      // } else {
      //   for (i = 0; i < mDim; ++i) x.push_back(runif(1, mBinfs[i], mBsups[i])[0]);
      // }
      x = runif(mDim, -0.5, 0.5);
      // for (i = 0; i < mDim; ++i)  x[i] = dU(gen);

      nv = 0;
      // Basis functions eval. at proposed point:
      // std::cout<<"Basis functions eval. at proposed point"<<std::endl;
      for (i = 0; i < n; ++i) {
        K = mIndex[i];
        tpcplx = computeFourierbasis(K, x);
        v[i] = tpcplx;
      }



      scal = 0;
      // Vector of projection weights (has length n-it)
      //       wei <- t(v)%*%estar

      // std::cout << "Compute vector of projection weights" << std::endl;
      // itwei = wei.begin();
      for (i = 0; i < n-it; ++i) {
        tpcplx = 0;
        etpstar = estar[i];
        for (j = 0; j < n; ++j) {
          tpcplx += v[j]*etpstar[j];
        }
        wei[i] = tpcplx;
        nv += pow(abs(tpcplx), 2.);
      }

      // accept = 1-mPws*nv/n;     // Accept probability
      accept = 1-nv/n;


      tp = runif (1, 0., 1.)[0];
      if (tp < accept) stop = true;

      ++tries;
      if (tries > rejectmax) std::cerr << "Rejection sampling failed reject_max =" << rejectmax << "times in a row" << std::endl;

    }

    // res.push_back(x);      // Record the accepted point
    // std::cout << "Record the accepted point" << std::endl;

    // std::cout << "Initialy: x = ("<< x[0] << ", "<< x[1]<<")"<<std::endl;
    // for(i = 0; i < mDim; ++i) {
    //   x[i] = x[i]*mWscale[i]+mWcenter[i];
    // }
    // std::cout << "After rescloing: x = ("<< x[0] << ", "<< x[1]<<")"<<std::endl;
    res(n-it,_) = x;

    if (it > 1) {      // while it is not the last point

      // Calculate orthogonal vector for Gram-Schmidt procedure:
      //  w <- v - rowSums(matrix(wei,n,n-i,byrow=TRUE)*e[,1:(n-i)])
      double nw = 0;
      // std::cout<<"Calculate orthogonal vector for Gram-Schmidt procedure 1"<<std::endl;
      for (i = 0; i < n; ++i) {
        tpcplx = 0;
        j = 0;
        for (j = 0; j < n-it; ++j) {
          tpcplx += wei[j]*e[j][i];
        }
        tpcplx = v[i]-tpcplx;
        w[i] = tpcplx;
        nw += pow(abs(tpcplx), 2);
      }

      nw = pow(nw, 0.5);
      // std::cout<<"Calculate orthogonal vector for Gram-Schmidt procedure 2"<<std::endl;
      for (i = 0; i < n ; ++i) {
        tpcplx = w[i]/nw;
        etp[i] = tpcplx;
        etpstar[i] = std::conj(tpcplx);
      }

      e.push_back(etp); estar.push_back(etpstar);
    }

  }

  return res;

}


void dpp_Eig::computeIndex() {

  // std::vector<int> coord (mDim, 0);              // vector of a possible permutations of {-k,...,k} initialized with -k
  // std::vector<double> temp(mDim, 0.);                // subvector of K coresponding to coord

  // std::cout<<"In ComputeIndex"<<std::endl;
  if (mIsProj) {
    this->setIndex(mIndextot);
  } else {

    std::vector< std::vector<int> > res;

    // std::vector<double> temp = mEig;

    double tp, tpbool;

    for (int i = 0; i < mIndextot.size(); ++i) {
      // tp = temp[i];
      // std::cout<<"In ComputeIndex"<<std::endl;
      tp = mEig[i];
      // std::cout<<"mEig[" <<i<< "]=" << tp <<std::endl;
      tpbool = rbinom(1, 1., tp)[0];  // Bernoulli draw
      // std::cout<<"tpbool=" << tpbool <<std::endl;
      if (tpbool == 1) {
        // coord = mIndextot[i];
        res.push_back(mIndextot[i]);
      }
    }

    this->setIndex(res);
  }
}














//
// NumericMatrix dpp_Gauss::testSample() {
//
//     // int rejectmax = 1e4;          // Max number of repetition for reject algorithm
//     // std::cout<< "Compute Index" << std::endl;
//     this->initIndex();
//     // this->computeIndex(k);
//     int n = mIndex.size();
//     // std::cout<< "Number of points to be computed -> mIndex.size() = " << n << std::endl;
//     int i, j, it;
//
//     // std::vector<std::vector<double> > res;
//     NumericMatrix res (n, mDim);    // matrix of coordinates
//
//     // First point: uniformaly distributed
//     // NumericVector x = runif (mDim, 0., 1.);
//     NumericVector x(3);
//
//     x[0] = 0.9478397; x[1] = 0.8304194; x[2] = 0.6271219;
//
//     // res.push_back(x);
//     // std::cout<< "Adding first point" << std::endl;
//     res(0,_) = x;
//
//     // std::cout<<"Done"<<std::endl;
//
//     std::vector<int> K(mDim);                   // Element of vector index
//
//     std::vector<std::complex<double> > v(n);      // Vector of eigenvector evaluated at a given point
//     std::vector<std::complex<double> > w(n);      // Temp vector for each computation
//     std::vector<std::complex<double> > wei;            // Temp vector of weight;
//     std::vector<std::complex<double> > etp(n);    // Temp vector e_i
//     std::vector<std::complex<double> > etpstar(n);   // Conj of e_i
//     std::vector<std::vector<std::complex<double> > > e;    // All vectors e_i
//     std::vector<std::vector<std::complex<double> > > estar;   // All conj of vectors e_i
//
//     std::complex<double> tpcplx;
//
//     double scal = 0, nv = 0;
//
//     // double accept;    // Accept probability
//
//     // std::cout<<"Basis-functions evaluated at first point"<<std::endl;
//     for (i = 0; i < n; ++i) {       // Basis-functions evaluated at first point
//       K = mIndex[i];
//       tpcplx = computeFourierbasis(K, x);
//       v[i] = tpcplx;
//       nv += pow(abs(tpcplx), 2.);
//     }
//     std::cout << "First v" << std::endl;
//     print_vector(v);
//     std::cout << std::endl;
//
//     nv = pow(nv, 0.5);
//     // std::cout<<"Record normalized version in the Gram-Schmidt matrices"<<std::endl;
//
//     for (int i = 0; i < mDim; ++i ) {       // Record normalized version in the Gram-Schmidt matrices:
//       etp[i] = v[i]/nv;
//       etpstar[i] = std::conj(etp[i]);
//     }
//     std::cout << "First e" << std::endl;
//     print_vector(etp);
//     std::cout << std::endl;
//     e.push_back(etp); estar.push_back(etpstar);
//
//     // int tries;
//     // std::cout<< "Main for loop over number of points"<< std::endl;
//     for (it = n-1; it > 0; --it) {    // Main for loop over number of points:
//       // std::cout << "Compute point n° " << it << std::endl;
//       // tries = 0;
//       // bool stop = false;
//       std::cout << "--------------" << std::endl;
//       std::cout << "Point n°" << n-i+1 << std::endl;
//       wei.resize(n-it);  // Weight vector at it-th step
//
//       // std::cout << "Acceptance algo" << std::endl;
//       // while(tries < rejectmax && !stop) {
//         // ++tries;
//         // std::cout << "Try n° " << tries << std::endl;
//         // Proposed point
//         // x = runif (mDim, 0., 1.);
//         if (it == n-1) {
//           x[0] = 0.9590675; x[1] = 0.5749953; x[2] = 0.7375211;
//         } else if (it == n-2) {
//           x[0] = 0.5207400; x[1] = 0.8601489; x[2] = 0.5395330;
//         }
//         // for (i = 0; i < mDim; ++i)  x[i] = dU(gen);
//
//         nv = 0;
//         // Basis functions eval. at proposed point:
//         // std::cout<<"Basis functions eval. at proposed point"<<std::endl;
//         for (i = 0; i < n; ++i) {
//           K = mIndex[i];
//           tpcplx = computeFourierbasis(K, x);
//           v[i] = tpcplx;
//         }
//         std::cout << "Current v" << std::endl;
//         print_vector(v);
//         std::cout << std::endl;
//
//
//
//         scal = 0;
//         // Vector of projection weights (has length n-it)
//         //       wei <- t(v)%*%estar
//
//         // std::cout << "Compute vector of projection weights" << std::endl;
//         // itwei = wei.begin();
//         for (i = 0; i < n-it; ++i) {
//           tpcplx = 0;
//           etpstar = estar[i];
//           for (j = 0; j < n; ++j) {
//             tpcplx += v[j]*etpstar[j];
//           }
//           wei[i] = tpcplx;
//         }
//         std::cout << "Current wei" << std::endl;
//         print_vector(wei);
//         std::cout << std::endl;
//
//
//         // accept = 1-nv/n;     // Accept probability
//
//         // if (runif (1, 0., 1.)[0] < accept) stop = true;
//         // stop = true;
//       // }
//
//       // res.push_back(x);      // Record the accepted point
//       // std::cout << "Record the accepted point" << std::endl;
//       res(n-it,_) = x;
//
//       if (it > 1) {      // while it is not the last point
//
//         // Calculate orthogonal vector for Gram-Schmidt procedure:
//         //  w <- v - rowSums(matrix(wei,n,n-i,byrow=TRUE)*e[,1:(n-i)])
//         double nw = 0;
//         // std::cout<<"Calculate orthogonal vector for Gram-Schmidt procedure 1"<<std::endl;
//         for (i = 0; i < n; ++i) {
//           tpcplx = 0;
//           j = 0;
//           for (j = 0; j < n-it; ++j) {
//             tpcplx += wei[j]*e[j][i];
//           }
//           tpcplx = v[i]-tpcplx;
//           w[i] = tpcplx;
//           nw += pow(abs(tpcplx), 2);
//         }
//
//         std::cout << "Current ortho vector w" << std::endl;
//         print_vector(w);
//         std::cout << std::endl;
//
//         nw = pow(nw, 0.5);
//         // std::cout<<"Calculate orthogonal vector for Gram-Schmidt procedure 2"<<std::endl;
//         for (i = 0; i < n ; ++i) {
//           tpcplx = w[i]/nw;
//           etp[i] = tpcplx;
//           etpstar[i] = std::conj(tpcplx);
//         }
//         std::cout << "Current e" << std::endl;
//         print_vector(etp);
//         std::cout << std::endl;
//
//         e.push_back(etp); estar.push_back(etpstar);
//       }
//
//     }
//
//     return res;
//
//
//
// }
