# require(Rcpp)
# require(spatstat)

# sourceCpp("~/Documents/GitHub/rcdpp/src/rcpp_rcdpp_module.cpp", verbose = TRUE)

######### Sampling DPPs given a model.
### - nsim: number of samples
### - d: dimension of the DPP
### - param: list of parameters. Depends on model:
###         * If model = "G", 'param' must have the following form
###           param = list(rho, alpha)
###           where 'rho' is the intensity of the DPP,
###           and alpha is a scale parameter which must satisfy
###           rho > (sqrt(pi)*alpha)^(-d)
###         * If model = "L1E", 'param' must have the following form
###           param = list(rho, alpha)
###           where 'rho' is the intensity of the DPP,
###           and alpha is a scale parameter which must satisfy
###           rho > (2*alpha)^(-d)
###         * If model = "D", 'param' must have the following form
###           param = list(N)
###           where 'N' is a vector of length 'd'
###           such that prod(N) is the deterministic number of points of the DPP.
###         * If model = "Eig", 'param' must have the following form
###           param = list(eigen, index)
###           where 'eigen' is a vector of eigenvalues (numbers non-greater than 1)
###           and 'index' is a matrix where index[i, ] is the vector in Z^d
###           where eigen[i] is located
### - domain: domain of observation
### - progress: if > 0, print a message each 'progress' simulated points
### - progress.sim: if > 0, print a message each 'progress.sim' simulation
### - with.kernel: boolean. If TRUE (default), the kernel is evaluated for each
###                couple of points in the simulated point configuration(s).
###                is returned as a symmetric matrix K, with indexes corresponding to those
###                of point configuration: K[i,j] is the kernel evaluated at
###                i-th and jth points.
rdppC <- function(nsim = 1, d, param = NULL, model = c("D", "L1E", "G", "Eig"),
                  type = c("prod", "sum"),
                  domain = boxx(rep(list(0:1), d)),
                  progress = 0, progress.sim = 0,
                  with.kernel = FALSE) {

if (missing(type)) type <- "prod"
type <- match.arg(type)
model <- match.arg(model)
# cat("Calling buil function \n")
  dpp <- build(d = d, param = param, model = model, type = type, domain = domain,
              progress = progress, progress.sim = progress.sim,
              with.kernel = with.kernel, sim = TRUE)
# cat("Calling cpp-simulate function \n")
  res <-  dpp$simulate(nsim)

  # cat("Done: convert Matrix into ppx \n")

  res <- lapply(res, function (pp) {
    if (!with.kernel) {
      tp <- if (length(pp) == 1) {
        if (d == 2) {
          ppp(numeric(0), numeric(0), domain = as.owin(domain))
        } else if (d == 3){
          pp3(numeric(0), numeric(0), numeric(0), domain)
        } else ppx(replicate(d, numeric(0), simplify = TRUE), domain = domain)
      } else {
        ppx(pp, domain = domain, simplify = TRUE)
      }
    } else {
      tp <- if (length(pp[[1]]) == 1) {
        if (d == 2) {
          ppp(numeric(0), numeric(0), domain = as.owin(domain))
        } else if (d == 3){
          pp3(numeric(0), numeric(0), numeric(0), domain)
        } else ppx(replicate(d, numeric(0), simplify = TRUE), domain = domain)
      } else {
        ppx(pp[[1]], domain = domain, simplify = TRUE)
      }
      tp <- c(PP = list(tp), K = list(pp[[2]]))
    }
    tp
  })

  if (nsim == 1) res <- res[[1]]

  # cat("Done \n")
  res

}


####### Kernel computation
computeKernel <- function(pp, param = NULL, model = c("D", "G", "L1E", "Eig"), type = c("prod", "sum"), domain = NULL) {

  if (is.ppx(pp)) {
    domain <- pp$domain
    pp <- as.matrix(pp$data[[6]])
    d <- ncol(pp)
    r <- domain$ranges
    binfs <- as.numeric(r[1, ])
    bsups <- as.numeric(r[2, ])
    pp <- list(pp)
    np <- 1
  } else if (is.ppp(pp)) {
    domain <- as.boxx(pp$window)
    d <- 2
    pp <- cbind(pp$x, pp$y)
    r <- domain$ranges
    binfs <- as.numeric(r[1, ])
    bsups <- as.numeric(r[2, ])
    pp <- list(pp)
    np <- 1
  } else if (is.matrix(pp)){
    if (is.null(domain)) stop("If 'pp' is a matrix, 'domain' must be non-empty.")
    else if (is.owin(domain)) domain <- as.boxx(domain)
    r <- domain$ranges
    d <- ncol(r)
    binfs <- as.numeric(r[1, ])
    bsups <- as.numeric(r[2, ])
    if (ncol(pp) != d) stop("'ncol(pp)' is not equal to dimension of 'domain'.")
    if (any(apply(pp, 2, min) < binfs)) stop("Some 'pp' points are out of 'domain'.")
    if (any(apply(pp, 2, max) > bsups)) stop("Some 'pp' points are out of 'domain'.")
    pp <- list(pp)
    np <- 1
  } else if (is.list(pp)) {
    np <- length(pp)
    if (is.null(domain)) stop("If 'pp' is a list of matrices, 'domain' must be non-empty.")
    else if (is.owin(domain)) domain <- as.boxx(domain)
    r <- domain$ranges
    binfs <- as.numeric(r[1, ])
    bsups <- as.numeric(r[2, ])
    d <- ncol(r)
    tests = sapply(pp, function(tpp) {
      cbind(!is.matrix(tpp),
            ncol(tpp) != d,
            any(apply(tpp, 2, min) < as.numeric(r[1, ])),
            any(apply(tpp, 2, max) > as.numeric(r[2, ]))
            )
      })
    if (any(tests[1,])) stop("If 'pp' is a list, each component must be a matrix.")
    if (any(tests[2,])) stop("Some 'pp' ncols do not fit to 'domain'.")
    if (any(tests[3,]) | any(tests[4,])) stop("Points of some 'pp' elements are out of 'domain'.")
  } else stop("'pp' must be a matrix, a ppp/ppx object or a list of matrices.")

  dpp <- build(d = d, param = param, model = model, type = type, domain = domain,
                progress = 0, progress.sim = 0,
                with.kernel = TRUE, sim = FALSE)

  # res <- lapply(pp, dpp$kernel)
  res <-  dpp$kernel(pp)
  # if (any(model == c("G", "L1E")) & !is.finite(k)) res <- lapply(res,  Re)

  if (np == 1) res <- res[[1]]

  res

}

####### PCF computation
computePCF <- function(pp, param = NULL, model = c("D", "G", "L1E", "Eig"), type = c("prod", "sum"), domain = NULL) {

    if (is.ppx(pp)) {
      domain <- pp$domain
      pp <- as.matrix(pp$data[[6]])
      d <- ncol(pp)
      r <- domain$ranges
      binfs <- as.numeric(r[1, ])
      bsups <- as.numeric(r[2, ])
      pp <- list(pp)
      np <- 1
    } else if (is.ppp(pp)) {
      domain <- as.boxx(pp$window)
      d <- 2
      pp <- cbind(pp$x, pp$y)
      r <- domain$ranges
      binfs <- as.numeric(r[1, ])
      bsups <- as.numeric(r[2, ])
      pp <- list(pp)
      np <- 1
    } else if (is.matrix(pp)){
      if (is.null(domain)) stop("If 'pp' is a matrix, 'domain' must be non-empty.")
      else if (is.owin(domain)) domain <- as.boxx(domain)
      r <- domain$ranges
      d <- ncol(r)
      binfs <- as.numeric(r[1, ])
      bsups <- as.numeric(r[2, ])
      if (ncol(pp) != d) stop("'ncol(pp)' is not equal to dimension of 'domain'.")
      if (any(apply(pp, 2, min) < binfs)) stop("Some 'pp' points are out of 'domain'.")
      if (any(apply(pp, 2, max) > bsups)) stop("Some 'pp' points are out of 'domain'.")
      pp <- list(pp)
      np <- 1
    } else if (is.list(pp)) {
      np <- length(pp)
      if (is.null(domain)) stop("If 'pp' is a list of matrices, 'domain' must be non-empty.")
      else if (is.owin(domain)) domain <- as.boxx(domain)
      r <- domain$ranges
      binfs <- as.numeric(r[1, ])
      bsups <- as.numeric(r[2, ])
      d <- ncol(r)
      tests = sapply(pp, function(tpp) {
        cbind(!is.matrix(tpp),
              ncol(tpp) != d,
              any(apply(tpp, 2, min) < as.numeric(r[1, ])),
              any(apply(tpp, 2, max) > as.numeric(r[2, ]))
              )
        })
      if (any(tests[1,])) stop("If 'pp' is a list, each component must be a matrix.")
      if (any(tests[2,])) stop("Some 'pp' ncols do not fit to 'domain'.")
      if (any(tests[3,]) | any(tests[4,])) stop("Points of some 'pp' elements are out of 'domain'.")
    } else stop("'pp' must be a matrix, a ppp/ppx object or a list of matrices.")

  dpp <- build(d = d, param = param, model = model, type = type, domain = domain,
                progress = 0, progress.sim = 0,
                with.kernel = TRUE, sim = FALSE)

  # res <- lapply(pp, dpp$kernel)
  res <-  dpp$pcf(pp)
  # if (any(model == c("G", "L1E")) & !is.finite(k)) res <- lapply(res,  Re)

  if (np == 1) res <- res[[1]]

  res

}


#### Construction modèles
build <- function(d, param = NULL, model = c("D", "L1E", "G", "Eig"),
                  type = c("prod", "sum"),
                  domain = boxx(rep(list(0:1), d)),
                  progress = 0, progress.sim = 0,
                  with.kernel = FALSE, sim = TRUE) {

  # cat("In build function \n")
  if (d <= 0) stop("'d' must be a positive integer.")

  if (missing(type)) type <- "prod"

  model <- match.arg(model)

  if (is.owin(domain)) domain <- as.boxx(domain)
  r <- domain$ranges
  # d <- ncol(r)
  if (d != ncol(r)) stop("Dimension of 'domain' must be equal to 'd'.")
  # }

  if (progress < 0) stop("'progress' must be a non-negative integer.")
  if (progress.sim < 0) stop("'progress.sim' must be a non-negative integer.")

  binfs <- as.numeric(r[1, ])
  bsups <- as.numeric(r[2, ])
  wsc <- bsups-binfs
  wc <- as.numeric(colMeans(r))
  # d <- ncol(r)
  args <- list(dim = d,
               # binfs = as.numeric(r[1, ]), bsups = as.numeric(r[2, ]),
               Wscale = wsc,
               Wcenter = wc,
               ic = (min(binfs) == max(binfs)) & (min(bsups) == max(bsups)),
               progress = progress, simprogress = progress.sim,
               wk = with.kernel
  )

  switch(model,
    'G' = {
      if (is.null(param)) stop("'G' model requires 3 parameters ('rho', 'alpha' and 'k').")
      if (length(param) != 3) stop("'G' model requires 3 parameters ('rho', 'alpha' and 'k').")
      npar <- names(param)
      if (!is.null(npar)) {
        if (!any(npar == "rho") | !any(npar == "alpha") | !any(npar == "k")) stop("If non-empty, 'names(param)' must contains 'rho', 'alpha' and 'k'.")
        rho <- param$rho
        alpha <- param$alpha
        k <- param$k
      } else {
        rho <- param[[1]]
        alpha <- param[[2]]
        k <- param[[3]]
      }
      if (rho < 0) stop("'rho' must be non-negative.")
      if (alpha < 0) stop("'alpha' must be non-negative.")
      if (sim) {
        if (floor(k) != k | k <= 0 | !is.finite(k)) stop("'k' must be a finite positive integer (simulation).")
      } else {
        if (floor(k) != k | k <= 0) stop("'k' must be a finite positive integer or infinite (exact kernel computation).")
        if (!is.finite(k)) k <- -1
      }

      if (rho > (sqrt(pi)*alpha)^(-d)) stop("'G' model is not valid (see  'help(rdppC)'.")

      # if (length(param) > 3) stop("'G' model requires only 3 parameters ('rho' and 'alpha').")
      args <- c(args, model = "G", type = type, rho = rho, alpha = alpha, k = k)
      # args <- c(args, model = "G", type = type, param = param)

      dpp <- new(dppComb, args)
    }
    ,
    'L1E' = {
      if (is.null(param)) stop("'L1E' model requires 3 parameters ('rho', 'alpha' and 'k').")
      if (length(param) != 3) stop("'L1E' model requires 3 parameters ('rho', 'alpha' and 'k').")
      npar <- names(param)
      if (!is.null(npar)) {
        if (!any(npar == "rho") | !any(npar == "alpha") | !any(npar == "k")) stop("If non-empty, 'names(param)' must contains 'rho', 'alpha' and 'k'.")
        rho <- param$rho
        alpha <- param$alpha
        k <- param$k
      } else {
        rho <- param[[1]]
        alpha <- param[[2]]
        k <- param[[3]]
      }
      if (rho < 0) stop("'rho' must be non-negative.")
      if (alpha < 0) stop("'alpha' must be non-negative.")
      if (sim) {
        if (floor(k) != k | k <= 0 | !is.finite(k)) stop("'k' must be a finite positive integer (simulation).")
      } else {
        if (floor(k) != k | k <= 0) stop("'k' must be a finite positive integer or infinite (exact kernel computation).")
        if (!is.finite(k)) k <- -1
      }

      if (rho > (2*alpha)^(-d)) stop("'L1E' model is not valid(see 'help(rdppC)').")
      args <- c(args, model = "L1E", type = type, rho = rho, alpha = alpha, k = k)
      # args <- c(args, model = "L1E", type = type, param = param)

      dpp <- new(dppComb, args)
    },
    'D' = {
      if (is.null(param)) stop("'D' model requires 1 parameter ('N').")
      if (length(param) != 1) stop("'D' model requires 1 parameter ('N').")
      npar <- names(param)
      if (!is.null(npar)) {
        if (!any(npar == "N")) stop("If non-empty, 'names(param)' must contains 'N'.")
        N <- param$N
      } else N <- param[[1]]
      if(any(N <= 0)) stop("'N' must be a vector of positive integers.")
      if (length(N) > d) stop("Size of 'N' must be non-greater than 'd'.")

      if (length(N) < d) {
        warning("Size of 'N' is smaller than 'd': 'N' is completed with '1'.")
        N <- c(N, rep(1, d-length(N)))
      }
      if (d > 1) args$ic <- args$ic & (sd(N) == 0)
      args <- c(args, model = model, type = type, N = 0, k = max(N))
      args$N <- N
      dpp <- new(dppDir, args)
    },

    'Eig' = {
      if (is.null(param)) stop("'Eig' model requires 2 parameters ('eigen', and 'index').")
      if (length(param) != 2) stop("'Eig' model requires 2 parameters ('eigen', and 'index').")
      npar <- names(param)
      if (!is.null(npar)) {
        if (!any(npar == "eigen") | !any(npar =="index")) stop("If non-empty, 'names(param)' must contains 'eigen' and 'index'.")
        eigen <- param$eigen
        index <- param$index
      } else {
        eigen = param[[1]]
        index = param[[2]]
      }
      if(!is.matrix(index)) stop("'index' must be a matrix with 'd' columns.")
      if (ncol(index) != d) stop("'index' must be a matrix with 'd' columns.")
      if (any(eigen > 1)) stop("'eigen' must be a vector of numbers non-greater than 1.")
      if (length(eigen) != nrow(index)) stop("'eigen' length must be equal to 'index' number of rows.")
      w <- which(eigen == 0)
      if (length(w) > 0) {
        eigen <- eigen[-w]
        index <- index[-w, ]
      }
      args <- c(args, eigenval = 0, index = 0,
                ip = sum(eigen == 1) == length(eigen))
      args$eigenval = eigen
      args$index = index
      dpp <- new(dppEig, args)
    }
  )


}




###### K-th intensity function computation
# computekInt <- function(pp, eigen, index, domain = boxx(rep(list(0:1), ncol(index)))) {
#
#   if (!is.ppx(pp) & !is.ppp(pp) & !is.matrix(pp)) stop("'pp' must be a matrix or a ppp/ppx object.")
#   if (!is.matrix(index)) stop("'index' must be a matrix.")
#   # else {
#   #   nd <- sapply(index, length)
#   #   if (var(nd) != 0) stop("'index' must be a list of numeric vectors with same size.")
#   # }
#   if (length(eigen) != nrow(index)) stop("'eigen' and 'index' must have same length.")
#
#   w <- which(eigen == 0)
#   if (length(w) > 0) {
#     eigen <- eigen[-w]
#     index <- index[-w]
#   }
#   # names(param) <- c("eigenval", "index")
#   # if (is.null(domain)) {
#   #   d <- ncol(index)
#   #   domain <- boxx(rep(list(0:1), d))
#   #   r <- domain$ranges
#   # } else {
#   domain <- as.boxx(domain)
#   r <- domain$ranges
#   d <- ncol(r)
#   if (d != ncol(index)) stop("Dimension of 'domain' must be equal to size of numeric vectors of 'index'.")
#   # }
#
#   binfs <- as.numeric(r[1, ])
#   bsups <- as.numeric(r[2, ])
#   wsc <- bsups-binfs
#   wc <- as.numeric(colMeans(r))
#   # d <- ncol(r)
#   args <- list(dim = as.integer(d),
#                # binfs = as.numeric(r[1, ]), bsups = as.numeric(r[2, ]),
#                Wscale = wsc,
#                Wcenter = wc,
#                # ic = (min(binfs) == max(binfs)) & (min(bsups) == max(bsups)),
#                progress = 0, simprogress = 0,
#                eigenval = eigen, index = index, ip = sum(eigen == 1) == length(eigen))
#   # cat('length(args) = ', length(args), '\n')
#   # cat('names(args) = ', names(args), '\n')
#
#   dpp <- new(dppEig, args)
#
#   res <- if (is.matrix(pp)) {
#     dpp$rhok(pp)
#   } else if (is.ppp(pp)) {
#     dpp$rhok(cbind(pp$x, pp$y))
#   } else if (is.ppx(pp)){
#     dpp$rhok(as.matrix(pp$data[[6]]))
#   }
#
#   res
# }
#
#
# computeEigen <- function(d, param = NULL, model = c("G", "L1E", "D"),
#                          domain = boxx(rep(list(0:1), d))) {
#
#
#     # if (rho < 0) stop("'rho' must be a non-negative number.")
#     if (d <= 0) stop("'d' must be a positive integer.")
#     # if (!is.null(k)) {
#     #   if (k < 0) {
#     #   # warning("'k' is a negative integer: its absolute value is taken.")
#     #     k <- abs(k)
#     #   }
#     # }
#
#
#     # if (!is.null(param) & !is.list(param)) stop("'param' must be a list.")
#
#     if (missing(model)) model <- "D"
#     model <- match.arg(model)
#
#     domain <- as.boxx(domain)
#     r <- domain$ranges
#     if (ncol(r) != d) stop("Dimension of 'domain' is different from 'd'.")
#
#     binfs <- as.numeric(r[1, ])
#     bsups <- as.numeric(r[2, ])
#     wsc <- bsups-binfs
#     wc <- as.numeric(colMeans(r))
#     args <- list( # rho = rho,
#       dim = d,
#       Wscale = wsc,
#       Wcenter = wc,
#       ic = (min(binfs) == max(binfs)) & (min(bsups) == max(bsups)),
#       progress = 0, simprogress = 0,
#       wk = TRUE
#     )
#
#     switch(model,
#       'G' = {
#         if(is.null(param)) {
#           stop("'G' model requires 3 parameters ('rho', 'alpha'and 'k').")
#           # args <- c(args, alpha = 1/(sqrt(pi)*rho^(1/d)))
#         } else {
#           if (length(param) > 3) stop("'G' model requires only 3 parameters ('rho', 'alpha'and 'k').")
#           if (param[[1]] > (sqrt(pi)*param[[2]])^(-d)) stop("'G' model is not valid.")
#           args <- c(args, rho = param[[1]], alpha = param[[2]], k = param[[3]])
#         }
#         dpp <- new(dppGauss, args)
#       }
#       ,
#       'L1E' = {
#         if (is.null(param)) {
#           # args <- c(args, alpha = 1/(2*rho^(1/d)))
#           stop("'L1E' model requires 3 parameters ('rho', 'alpha'and 'k').")
#         } else {
#           if (length(param) > 3) stop("'L1E' model requires only 3 parameters ('rho', 'alpha'and 'k').")
#           if (param[[1]] > (2*param[[2]])^(-d)) stop("'L1E' model is not valid.")
#           args <- c(args, rho = param[[1]], alpha = param[[2]], k = param[[3]])
#         }
#         dpp <- new(dppL1Exp, args)
#       },
#       'D' = {
#         if (is.null(param)) {
#           stop("'D' model requires 1 parameter ('N').")
#         } else {
#           N <- param[[1]]
#           if(sum(N <= 0) > 0) stop("'N' must be a vector of positive number.")
#           if (length(N) > d) stop("Size of 'N' must be non-greater than 'd'.")
#
#
#           if (length(N) < d) {
#             warning("Size of 'N' is smaller than 'd': 'N' is completed with '1'.")
#             N <- c(N, rep(1, d-length(N)))
#           }
#           if (d > 1) args$ic <- args$ic & (sd(N) == 0)
#           args <- c(args, N = 0, k = max(N))
#           args$N <- N
#           dpp <- new(dppDir, args)
#         }
#       }
#     )
#
#     res <- dpp$eigen()
#
#     res
#
# }
#
# computeIndex <- function(d, param = NULL, model = c("G", "L1E", "D"), domain = boxx(rep(list(0:1), d))) {
#
#
#     # if (rho < 0) stop("'rho' must be a non-negative number.")
#     if (d <= 0) stop("'d' must be a positive integer.")
#
#
#     # if (!is.null(param) & !is.list(param)) stop("'param' must be a list.")
#
#     if (missing(model)) model <- "G"
#     model <- match.arg(model)
#
#     domain <- as.boxx(domain)
#     r <- domain$ranges
#     if (ncol(r) != d) stop("Dimension of 'domain' is different from 'd'.")
#
#     # binfs <- as.numeric(r[1, ])
#     # bsups <- as.numeric(r[2, ])
#     wsc <- as.numeric(r[2,]-r[1,])
#     wc <- as.numeric(colMeans(r))
#     args <- list(dim = d,
#                  Wscale = wsc,
#                  Wcenter = wc,
#                  ic = TRUE,
#                  progress = 0, simprogress = 0
#                 )
#
#     switch(model,
#       'G' = {
#         if(is.null(param)) {
#           stop("'G' model requires 3 parameters ('rho', 'alpha'and 'k').")
#           # args <- c(args, alpha = 1/(sqrt(pi)*rho^(1/d)))
#         } else {
#           if (length(param) > 3) stop("'G' model requires only 3 parameters ('rho', 'alpha'and 'k').")
#           if (param[[1]] > (sqrt(pi)*param[[2]])^(-d)) stop("'G' model is not valid.")
#           args <- c(args, rho = param[[1]], alpha = param[[2]], k = param[[3]])
#         }
#         dpp <- new(dppGauss, args)
#       }
#       ,
#       'L1E' = {
#         if (is.null(param)) {
#           # args <- c(args, alpha = 1/(2*rho^(1/d)))
#           stop("'L1E' model requires 3 parameters ('rho', 'alpha'and 'k').")
#         } else {
#           if (length(param) > 3) stop("'L1E' model requires only 3 parameters ('rho', 'alpha'and 'k').")
#           if (param[[1]] > (2*param[[2]])^(-d)) stop("'L1E' model is not valid.")
#           args <- c(args, rho = param[[1]], alpha = param[[2]], k = param[[3]])
#         }
#         dpp <- new(dppL1Exp, args)
#       },
#       'D' = {
#         if (is.null(param)) {
#           stop("'D' model requires 1 parameter ('N').")
#         } else {
#           N <- param[[1]]
#           if(sum(N <= 0) > 0) stop("'N' must be a vector of positive number.")
#           if (length(N) > d) stop("Size of 'N' must be non-greater than 'd'.")
#
#
#           if (length(N) < d) {
#             warning("Size of 'N' is smaller than 'd': 'N' is completed with '1'.")
#             N <- c(N, rep(1, d-length(N)))
#           }
#           if (d > 1) args$ic <- args$ic & (sd(N) == 0)
#           args <- c(args, N = 0)
#           # cat("N =", args$N,"\n")
#           args$N <- N
#           # cat("N =", args$N,"\n")
#           # warning("Parameter 'k' is ignored and set to max(N) when 'model' is 'D'.")
#           # k <- max(N)-1
#
#           dpp <- new(dppDir, args)
#         }
#       }
#     )
#
#     # res <- sapply(0:(d-1), function(i) dpp$index(k,i))
#     res <- dpp$index()
#
#     res
#
# }
