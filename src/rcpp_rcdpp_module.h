#ifndef RCPP_RCDPP_MODULE_H
#define RCPP_RCDPP_MODULE_H

#include "rcdpp_Dir.h"
#include "rcdpp_Eig.h"
#include "rcdpp_Prod.h"

RCPP_EXPOSED_AS(dpp_All)
RCPP_EXPOSED_WRAP(dpp_All)

RCPP_EXPOSED_AS(dpp_Prod)
RCPP_EXPOSED_WRAP(dpp_Prod)

RCPP_EXPOSED_AS(dpp_Gauss)
RCPP_EXPOSED_WRAP(dpp_Gauss)

RCPP_EXPOSED_AS(dpp_L1Exp)
RCPP_EXPOSED_WRAP(dpp_L1Exp)

// RCPP_EXPOSED_AS(dpp_MR)
// RCPP_EXPOSED_WRAP(dpp_MR)

// RCPP_EXPOSED_AS(dpp_Dir0)
// RCPP_EXPOSED_WRAP(dpp_Dir0)

RCPP_EXPOSED_AS(dpp_Dir)
RCPP_EXPOSED_WRAP(dpp_Dir)


RCPP_EXPOSED_AS(dpp_Eig)
RCPP_EXPOSED_WRAP(dpp_Eig)

#endif
