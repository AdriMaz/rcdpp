# require(spatstat, quiet = TRUE)

#

rdppEigC <- function(nsim = 1, eigen, index, window = NULL, progress = 0, progress.sim = 0) {

    # Note: not necessarily to enter zero eigenvalues

    if (!is.list(index)) stop("'index' must be a list of numeric vectors with same size.")
    else {
      nd <- sapply(index, length)
      if (var(nd) != 0) stop("'index' must be a list of numeric vectors with same size.")
    }
    if (length(eigen) != length(index)) stop("'eigenval' and 'index' must have same length.")

    w <- which(eigen == 0)
    if (length(w) > 0) {
      eigen <- eigen[-w]
      index <- index[-w]
    }
    # names(param) <- c("eigenval", "index")
    if (is.null(window)) {
      d <- length(index[[1]])
      window <- boxx(rep(list(0:1), d))
      r <- window$ranges
    } else {
      window <- as.boxx(window)
      r <- window$ranges
      if (ncol(r) != length(index[[1]])) stop("Dimension of 'window' must be equal to size of numeric vectors of 'index'.")
    }


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
                 # ic = (min(binfs) == max(binfs)) & (min(bsups) == max(bsups)),
                 progress = progress, simprogress = progress.sim,
                 eigenval = eigen, index = index, ip = sum(eigen == 1) == length(eigen))
    # cat('length(args) = ', length(args), '\n')
    # cat('names(args) = ', names(args), '\n')

    dpp <- new(dppEig, args)

    res <-  dpp$simulate(nsim)
    res <- lapply(res, function (pp) {
      if (length(pp) == 1) {
        if (d == 2) {
          ppp(numeric(0), numeric(0), window = as.owin(window))
        } else if (d == 3){
          pp3(numeric(0), numeric(0), numeric(0), window)
        } else ppx(replicate(d, numeric(0), simplify=FALSE), domain = window)
      }
      else {
        # pp <- pp%*%diag(wsc)+wc
        if (d == 2) {
          ppp(pp[ ,1], pp[ ,2], window = as.owin(window))
        } else if (d == 3) {
          pp3(pp[ ,1], pp[ ,2], pp[ ,3], domain = window)
        } else ppx(pp, domain = window)
      }


    })
    if (nsim == 1) res <- res[[1]]

    # cat("Done \n")
    res

}

rdppC <- function(rho, k, d, nsim = 1, param = NULL, model = c("G", "L1G", "MR", "MRProd"), window = boxx(rep(list(0:1), d)), progress = 0, progress.sim = 0) {


  if (rho < 0) stop("'rho' must be a non-negative number.")
  if (d <= 0) stop("'d' must be a positive integer.")
  if (k < 0) {
    # warning("'k' is a negative integer: its absolute value is taken.")
    k <- abs(k)
  }


  # if (!is.null(param) & !is.list(param)) stop("'param' must be a list.")

  if (missing(model)) model <- "G"
  model <- match.arg(model)

  window <- as.boxx(window)
  r <- window$ranges
  if (ncol(r) != d) stop("Dimension of 'window' is different from 'd'.")

  if (progress < 0) stop("'progress' must be a non-negative integer.")
  if (progress.sim < 0) stop("'progress.sim' must be a non-negative integer.")

  binfs <- as.numeric(r[1, ])
  bsups <- as.numeric(r[2, ])
  wsc <- bsups-binfs
  wc <- as.numeric(colMeans(r))
  args <- list(rho = rho, dim = d,
               # binfs = as.numeric(r[1, ]), bsups = as.numeric(r[2, ]),
               Wscale = wsc,
               Wcenter = wc,
               ic = (min(binfs) == max(binfs)) & (min(bsups) == max(bsups)),
               progress = progress, simprogress = progress.sim
              )
  # cat('length(args) = ', length(args), '\n')
  # cat('names(args) = ', names(args), '\n')
  # if(!is.null(param)) {
  #
  #   args <- c(args, param)
  # }
  switch(model,
    'G' = {
      if(is.null(param)) {
        args <- c(args, alpha = 1/(sqrt(pi)*rho^(1/d)))
      } else {
        if (length(param) > 1) stop("'G' model requires only 1 parameter ('alpha').")
        args <- c(args, alpha = param[[1]])
      }
      dpp <- new(dppGauss, args)
    }
    ,
    'L1G' = {
      if (is.null(param)) {
        args <- c(args, alpha = 1/(2*rho^(1/d)))
      } else {
        if (length(param) > 1) stop("'param' is missing: 'L1G' model requires only 1 parameter ('alpha').")
        args <- c(args, alpha = param[[1]])
      }
      dpp <- new(dppL1Gauss, args)
    },
    'MR' = {
      if (is.null(param)) {
        # stop("'param' is missing: 'MR' model requires 1 parameter ('tau').")
        tau <- if (d == 1) rho/2 else if (d == 2) rho/pi else (rho*gamma(d/2+1)/(pi^(d/2)))^(2/d)
        args <- c(args, tau = tau)
      } else {
        if (length(param) > 1) stop("'MR' model requires only 1 parameter ('tau').")
        args <- c(args, tau = param[[1]])
      }
      dpp <- new(dppMR, args)
    },
    'MRProd' = {
      if (!is.null(param)) warning("'MRProd' model does not require parameter: 'param' is ignored." )
      dpp <- new(dppMRProd, args)
      tp <- ((rho*prod(wsc))^(1/d)-1)/2
      kmax <- floor(round(tp, 8))
      # if (tp%%2 != 1) warning("'rho' is not a 'd'-power of an odd number: the number of points will not be equal to 'int'.")
      if (abs(tp - kmax) > 1e8) warning("'rho' is not a 'd'-power of an odd number: the number of points will not be equal to 'int'.")
      k <- min(k, kmax)
    }
  )

  res <-  dpp$simulate(k, nsim)

  # cat("Done: convert Matrix into ppx \n")

  res <- lapply(res, function (pp) {
    if (length(pp) == 1) {
      if (d == 2) {
        ppp(numeric(0), numeric(0), window = as.owin(window))
      } else if (d == 3){
        pp3(numeric(0), numeric(0), numeric(0), window)
      } else ppx(replicate(d, numeric(0), simplify=FALSE), domain = window)
    }
    else {
      # pp <- pp%*%diag(wsc)+wc
      if (d == 2) {
        ppp(pp[ ,1], pp[ ,2], window = as.owin(window))
      } else if (d == 3) {
        pp3(pp[ ,1], pp[ ,2], pp[ ,3], domain = window)
      } else ppx(pp, domain = window)
    }


  })
  if (nsim == 1) res <- res[[1]]

  # cat("Done \n")
  res


}

computeEigen <- function(rho, k, d, param = NULL, model = c("G", "L1G", "MR", "MRProd")) {


    if (rho < 0) stop("'rho' must be a non-negative number.")
    if (d <= 0) stop("'d' must be a positive integer.")
    if (k < 0) {
      # warning("'k' is a negative integer: its absolute value is taken.")
      k <- abs(k)
    }


    # if (!is.null(param) & !is.list(param)) stop("'param' must be a list.")

    if (missing(model)) model <- "G"
    model <- match.arg(model)

    window <- as.boxx(window)
    r <- window$ranges
    if (ncol(r) != d) stop("Dimension of 'window' is different from 'd'.")

    # binfs <- as.numeric(r[1, ])
    # bsups <- as.numeric(r[2, ])
    wsc <- as.numeric(r[2,]-r[1,])
    wc <- as.numeric(colMeans(r))
    args <- list(rho = rho, dim = d,
                 Wscale = wsc,
                 Wcenter = wc,
                 ic = TRUE,
                 progress = 0, simprogress = 0
                )

    switch(model,
      'G' = {
        if(is.null(param)) {
          args <- c(args, alpha = 1/(sqrt(pi)*rho^(1/d)))
        } else {
          if (length(param) > 1) stop("'G' model requires only 1 parameter ('alpha').")
          args <- c(args, alpha = param[[1]])
        }
        dpp <- new(dppGauss, args)
      }
      ,
      'L1G' = {
        if (is.null(param)) {
          args <- c(args, alpha = 1/(2*rho^(1/d)))
        } else {
          if (length(param) > 1) stop("'param' is missing: 'L1G' model requires only 1 parameter ('alpha').")
          args <- c(args, alpha = param[[1]])
        }
        dpp <- new(dppL1Gauss, args)
      },
      'MR' = {
        if (is.null(param)) {
          stop("'param' is missing: 'MR' model requires 1 parameter ('tau').")
        } else {
          if (length(param) > 1) stop("'MR' model requires only 1 parameter ('tau').")
          args <- c(args, tau = param[[1]])
        }
        dpp <- new(dppMR, args)
      },
      'MRProd' = {
        if (!is.null(param)) warning("'MRProd' model does not require parameter: 'param' is ignored." )
        dpp <- new(dppMRProd, args)
        tp <- ((rho*prod(wsc))^(1/d)-1)/2
        kmax <- floor(round(tp, 8))
        # if (tp%%2 != 1) warning("'rho' is not a 'd'-power of an odd number: the number of points will not be equal to 'int'.")
        if (abs(tp - kmax) > 1e8) warning("'rho' is not a 'd'-power of an odd number: the number of points will not be equal to 'int'.")
        k <- min(k, kmax)
      }
    )

    res <- dpp$eigen(k)

    res

}

#
# rdppGauss <- function(rho, d, alpha = 1/(sqrt(pi)*rho^(1/d)), k, nsim = 1
#   # , window = boxx(rep(list(0:1), d))
#               ) {
#
#    dppG <- new(dppGauss, list(rho = rho, alpha = alpha, dim = d
#      # , binfs = unlist(ranges[1, ]), bsups = unlist(ranges[2, ]))
#    ))
#
#   # if (nsim > 1)  {
#     res <- dppG$simulate(k, nsim)
#
#     res <- lapply(res, function (pp) {
#       if (length(pp) == 1) {
#         if (d == 2) {
#           ppp(numeric(0), numeric(0), window=owin(c(0, 1), c(0, 1)))
#         } else if (d == 3){
#           pp3(numeric(0), numeric(0), numeric(0), box3(rep(list(0:1), d)))
#         } else ppx(replicate(d, numeric(0), simplify=FALSE), domain = boxx(rep(list(0:1), d)))
#       }
#       else {
#         if (d == 2) {
#           ppp(pp[,1], pp[,2], window=owin(c(0, 1), c(0, 1)))
#         } else if (d == 3) {
#           pp3(pp[,1], pp[,2], pp[,3], domain = box3(rep(list(0:1), 3)))
#         } else ppx(pp, domain = boxx(rep(list(0:1), d)))
#       }
#
#
#     })
#     if (nsim == 1) res <- res[[1]]
#
#     res
# }
#
#
# rdppL1Gauss <- function(rho, d, alpha = 1/(2*rho^(1/d)), k, nsim = 1
#   # , window = boxx(rep(list(0:1), d))
#               ) {
#
#    dppE <- new(dppL1Gauss, list(rho = rho, alpha = alpha, dim = d
#      # , binfs = unlist(ranges[1, ]), bsups = unlist(ranges[2, ]))
#    ))
#
#   # if (nsim > 1)  {
#     res <- dppE$simulate(k, nsim)
#
#     res <- lapply(res, function (pp) {
#       if (length(pp) == 1) {
#         if (d == 2) {
#           ppp(numeric(0), numeric(0), window=owin(c(0, 1), c(0, 1)))
#         } else if (d == 3){
#           pp3(numeric(0), numeric(0), numeric(0), box3(rep(list(0:1), d)))
#         } else ppx(replicate(d, numeric(0), simplify=FALSE), domain = boxx(rep(list(0:1), d)))
#       }
#       else {
#         if (d == 2) {
#           ppp(pp[,1], pp[,2], window=owin(c(0, 1), c(0, 1)))
#         } else if (d == 3) {
#           pp3(pp[,1], pp[,2], pp[,3], domain = box3(rep(list(0:1), 3)))
#         } else ppx(pp, domain = boxx(rep(list(0:1), d)))
#       }
#
#
#     })
#     if (nsim == 1) res <- res[[1]]
#
#     res
# }
#
#
# rdppMR <- function(rho, d, k, nsim = 1) {
#
#   dppM <- new(dppMR, list(rho = rho, dim = d))
#
#
#   res <- dppE$simulate(k, nsim)
#
#   res <- lapply(res, function (pp) {
#     if (length(pp) == 1) {
#       if (d == 2) {
#         ppp(numeric(0), numeric(0), window=owin(c(0, 1), c(0, 1)))
#       } else if (d == 3){
#         pp3(numeric(0), numeric(0), numeric(0), box3(rep(list(0:1), d)))
#       } else ppx(replicate(d, numeric(0), simplify=FALSE), domain = boxx(rep(list(0:1), d)))
#     }
#     else {
#       if (d == 2) {
#         ppp(pp[,1], pp[,2], window=owin(c(0, 1), c(0, 1)))
#       } else if (d == 3) {
#         pp3(pp[,1], pp[,2], pp[,3], domain = box3(rep(list(0:1), 3)))
#       } else ppx(pp, domain = boxx(rep(list(0:1), d)))
#     }
#
#
#   })
#   if (nsim == 1) res <- res[[1]]
#
#   res
#
# }



## Compute eigen of -k:k
# eigGauss <- function(k, rho, d , alpha = 1/(sqrt(pi)*rho^(1/d))) {
#
#
#   dppG <- new(dppGauss, list(rho = rho, alpha = alpha, dim = d))
#
#   dppG$eigen(k)
#
# }
#

#
# testsim <- function() {
#   rho <- 200
#   d <- 3
#   alpha <- 1/(sqrt(pi)*rho^(1/d))
#   dppG <- new(dppGauss, list(rho = rho, alpha = alpha, dim = d))
#
#   dppG$test()
#
# }
#
#
