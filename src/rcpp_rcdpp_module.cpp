#include "rcpp_rcdpp_module.h"

RCPP_MODULE(rcdpp_module) {

  class_<dpp_All>("dppAll")
    // .constructor<List>()
    .method("simulate", &dpp_All::computeListSamples, "Simulate DPPs.")
    .method("kernel", &dpp_All::computeOnlyKernelR, "Evaluate kernel on point configuration.")
    .method("pcf", &dpp_All::computePCFR, "Evaluate pcf on point configuration.")
    // .field("eig", &dppAll::mEig, "eigen values")
    // .method("eigen", &dpp_All::getEigen, "get eigen values")
    // .method("index", &dpp_All::getIndex, "get index lattice")
    // .method("binom", &dpp_All::test_binom, "Test Rcpp::rbinom")
  ;

  class_<dpp_Prod>("dppProd")
    .derives<dpp_All>("dppAll")
    // .constructor<List>()
    // .method("eigen", &dpp_Prod::getEigen, "get eigen values")
  ;

  class_<dpp_Gauss>("dppGauss")
    .derives<dpp_Prod>("dppProd")
    .constructor<List>()
    // .method("simulate", &dppGauss::computeListSamples, "Simulate Gaussian DPPs.")
    // .method("test", &dppGauss::testSample, "Test sample function.")
    // .method("eigen", &dppGauss::computeEigen, "Compute Eigenvalues.")
  ;

  class_<dpp_L1Exp>("dppL1Exp")
    .derives<dpp_Prod>("dppProd")
    .constructor<List>()
  ;

  class_<dpp_Dir>("dppDir")
    .derives<dpp_Prod>("dppProd")
    .constructor<List>()
  ;


  class_<dpp_Eig>("dppEig")
    .derives<dpp_All>("dppAll")
    .constructor<List>()
    // .method("simulate", &dpp_Eig::computeListSamples, "Simulate DPPs given eigenvalues in Fourier basis.")
    // .method("kernel", &dpp_Eig::computeKernelR, "Evaluate kernels on point configuration.")
    // .method("rhok", &dpp_Eig::computekInt, "Evaluate kème intensity function of kxk matrix.")
    // .method("pcf", &dpp_Eig::computePCFR, "Evaluate pcfs on point configurations.")
  ;
}
