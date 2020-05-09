#include "rcdpp_All.h"

using namespace Rcpp;

template <typename T>
void print_vector(std::vector<T> v) {
  int n = v.size();
  for (int i = 0; i < n; ++i) Rcout << v[i] << " ";
  Rcout << std::endl;
}

void print_vector(NumericVector v) {
  int n = v.size();
  for (int i = 0; i < n; ++i) Rcout << v[i] << " ";
  Rcout << std::endl;
}



bool next_variation(std::vector<int>::iterator first, std::vector<int>::iterator last, const int min, const int max) {
// bool next_variation(std::vector<int> vect,  const int max) {
    // std::vector<int>::iterator first (vect.begin());
    // std::vector<int>::iterator last (vect.end());
    // double eps=1e-8;

    // std::cout<<"In next_variation"<<std::endl;
    // std::cout<<"first = "<<*first<<"; last = "<<*last<<std::endl;
    if (first == last) return false; // empty sequence (n==0)

    std::vector<int>::iterator i(last);
    --i;   // Point to the rightmost element

    // Check if I can just increase it
    // std::cout<<max<<"-"<<*i<<" = "<<max-(*i)<<std::endl;
    if (*i < max) {
    // if (max-(*i) > eps) {
      // std::cout<<"First 'if' satisfied"<<std::endl;
       ++(*i);   // Increase this element and return
       return true;
     }

    // Find the rightmost element to increase
    while (i != first) {
      // std::cout<<"In 'while' loop"<<std::endl;
      *i = min; // reset the right-hand element
      --i; // point to the left adjacent
      // std::cout<<max<<"-"<<*i<<" = "<<max-(*i)<<std::endl;
      if (*i < max) {
        // std::cout<<"'if' in the 'while' loop is satisfied"<<std::endl;
        ++(*i);   // Increase this element and return
        return true;
      }
    }
       // std::cout<<"Skip the 'while' loop"<<std::endl;
    // If here all elements are the maximum symbol (max=k-1), so there are no more variations
    //for (i=first; i!=last; ++i) *i = 0; // Should reset to the lowest sequence (0)?
    return false;
}

std::complex<double> computeFourierbasis(const std::vector<int>& k, const NumericVector& x, const NumericVector& boxlengths) {

  std::complex<double> dpii (0., 2.*M_PI);
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




List dpp_All::computeListSamples(const int nsim) {

  List RES;
  // List tpres;
  NumericMatrix res;
  NumericVector x;
  int i,j,l;

  // Set seed
  RNGScope scope;

  // Rescaling index if necessary

  if (mIsProj) {    // If projection DPP: compute index only once

    // std::cout<<"Call computeIndex (Proj. case)"<<std::endl;
    computeIndex();
    // std::cout<<"computeIndex done"<<std::endl;

    for (i = 0; i < nsim; ++i) {
      // std::cout << "Sim n° " << i+1 << std::endl;
      if(mProgSim > 0) {
        // if ((i+1)%mProgSim == 0) std::cout << "Compute sample n° " << i+1 << std::endl;
        if ((i+1)%mProgSim == 0) Rcout << "Compute sample n° " << i+1 << "\n";
      }

      res = computeSample();   // Compute sample

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
      if (mWithKernel) {
        List tpres;
        tpres.push_back(res);
        tpres.push_back(computeKernelR(res));

        RES.push_back(tpres);
      } else {
        RES.push_back(res);       // Stock new point process
      }
    }
  } else {
  // std::cout<< "Compute Eigen" << std::endl;
    // computeEigenVec();
    // std::cout<< "Done -> mEig.size = " << mEig.size() << std::endl;
    // print_vector(mEig);

    for (i = 0; i < nsim; ++i) {

      computeIndex();

      if(mProgSim > 0) {
        // if ((i+1)%mProgSim == 0) std::cout << "Compute sample n° " << i+1 << std::endl;
        if ((i+1)%mProgSim == 0) Rcout << "Compute sample n° " << i+1 << "\n";
      }
      // std::cout << "Sim n° " << i+1 << std::endl;
      res = computeSample();   // Compute sample

      for (l = 0; l < res.nrow(); ++l){
        x = res(l,_);
        for(j = 0; j < mDim; ++j) {
          x[j] = x[j]*mWscale[j]+mWcenter[j];
        }
        res(l,_) = x;
      }

      if (mWithKernel) {
        List tpres;
        tpres.push_back(res);
        tpres.push_back(computeKernelR(res));

        RES.push_back(tpres);
      } else {
        RES.push_back(res);       // Stock new point process
      }
      // std::cout << "Before reset  mIndex.size() = " << mIndex.size() << std::endl;
      resetIndex();       // reset mIndex
      // std::cout << "Reset Index -> mIndex.size() = " << mIndex.size() << std::endl;
    }
  }
  return RES;
}


NumericMatrix dpp_All::computeSample() {

  // Remark: Sample is computed on [-1/2 ; 1/2]^d
  // to be rescaled on domain W

  int rejectmax = 1e4;          // Max number of repetition for reject algorithm
  // std::cout<< "Compute Index" << std::endl;
  // if (!mIsProj) computeIndex(k);
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
      // if ((n-it+1)%mProg == 0) std::cout << "Compute point n° " << n-it+1 << std::endl;
      if ((n-it+1)%mProg == 0) Rcout << "Compute point n° " << n-it+1 << "\n";
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


List dpp_All::computeOnlyKernelR(const List& PP) {

  int np = PP.size();
  List RES;
  // if (!mIsEigSet) computeEigenForKernel();
  for (int i=0; i < np; ++i) {
    NumericMatrix pp = PP[i];
    // Rcout<<"nrow(pp) ="<<pp.nrow()<<" ; ncol(pp) ="<<pp.ncol()<<"\n";
    // for (int j =0; j < pp.nrow(); ++j) {
    //   NumericVector tp = pp(j,_);
    //   print_vector(tp);
    // }

    ComplexMatrix tpres = computeKernelR(pp);
    RES.push_back(tpres);
  }
  return RES;
}

List dpp_All::computePCFR(const List& PP) {

  int np = PP.size();
  List RES;
  // if (!mIsEigSet) computeEigenForKernel();
  for (int i=0; i < np; ++i) {
    NumericMatrix pp = PP[i];
    // Rcout<<"nrow(pp) ="<<pp.nrow()<<" ; ncol(pp) ="<<pp.ncol()<<"\n";
    // for (int j =0; j < pp.nrow(); ++j) {
      // NumericVector tp = pp(j,_);
      // print_vector(tp);
    // }

    NumericMatrix tpres = computePCF(pp);
    RES.push_back(tpres);
  }
  return RES;
}


NumericMatrix dpp_All::computePCF(const NumericMatrix& PP) {

  int np = PP.nrow();
  NumericMatrix res (np,np);
  ComplexMatrix K = computeKernelR(PP);
  // NumericVector tpX (mDim);
  // NumericVector tpY (mDim);
  Rcomplex tpc;
  double tp;
  // diag.r = mInt; diag.i = 0.;


  for (int i = 0; i < np; ++i) {
    res(i,i) = 0.;
    NumericVector tpX = PP(i,_);
    for (int j = 0; j < i; ++j) {
      // if (j != i) {
        // tpY = PP(j,_);
        tpc = K(i,j);
        tp = pow(tpc.r,2)+pow(tpc.i,2);
        tp /= pow(mInt,2);
        res(i,j) = tp; res(j,i) = tp;
      // }
    }
  }
  return res;

}
