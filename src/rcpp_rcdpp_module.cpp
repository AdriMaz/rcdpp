#include "rcpp_rcdpp_module.h"

RCPP_MODULE(rcdpp_module) {

  class_<dpp_All>("dppAll")
    // .constructor<List>()
    .method("simulate", &dpp_All::computeListSamples, "Simulate DPPs.")
    // .field("eig", &dppAll::mEig, "eigen values")
    .method("eigen", &dpp_All::getEigen, "get eigen values")
    // .method("binom", &dpp_All::test_binom, "Test Rcpp::rbinom")
  ;

  class_<dpp_Prod>("dppProd")
    .derives<dpp_All>("dppAll")
    // .constructor<List>()
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

  // class_<dpp_MR>("dppMR")
  //   .derives<dpp_All>("dppAll")
  //   .constructor<List>()
  // ;

  class_<dpp_Dir0>("dppDir0")
    .derives<dpp_All>("dppAll")
    .constructor<List>()
  ;

  class_<dpp_Dir>("dppDir")
    .derives<dpp_All>("dppAll")
    .constructor<List>()
  ;


  //
  class_<dpp_Eig>("dppEig")
    // .derives<dpp_All>("dppAll")
    .constructor<List>()
    .method("simulate", &dpp_Eig::computeListSamples, "Simulate DPPs given eigenvalues in Fourier basis.")
  ;
}
