\name{rdppC}
\alias{rdppC}
\encoding{UTF-8}

\title{Simulation of a Determinantal Point Process}

\description{
Generates simulated realisations from a determinantal point process model.
}

\usage{

  rdppC(nsim = 1, d, param = NULL,
        model = c("D", "L1E", "G", "Eig"),
        domain = boxx(rep(list(0:1), d)),
        progress = 0, progress.sim = 0,
        with.kernel = FALSE)
}

\arguments{
 \item{nsim}{ number of simulated realisations.}
 \item{d}{ dimension of the DPP.}
 \item{param}{ list of parameters. See Details.}
 \item{domain}{ domain of the DPP. Must be of class \code{"ownin"}, \code{"box3"} or \code{"boxx"}.}
 \item{model}{ Class of the DPP. See Details.}
 \item{progress}{ if > 0, print a message each \code{progress} simulated points.}
 \item{progress.sim}{ if > 0, print a message each \code{progress.sim} simulated realization.}
 \item{with.kernel}{ boolean. If \code{TRUE}, the kernel associated to the simulated DPP is computed for each couple of points. }
}
\details{
  A DPP is usually characterized by a kernel. The exact simulation of the DPP models requires the eigendecomposition
  of the associated kernel, which generaly involves an infinite serie with no analytical solution.

  The list \code{param} depends on \code{model} (the order of appearance of the components must be respected, or the user must set the names as below)

    \code{"G"} model stands for Gaussian DPP, for which \code{param} must have the following form
     \code{param = list(rho, alpha, k)}. A Gaussian DPP is a DPP with kernel
      \deqn{K(x,y) = \rho exp(-||(x-y)/\alpha||^2)}.
     where \eqn{\rho} is the intensity of the considered DPP and \eqn{alpha} a scale parameter which must satisfy
      \eqn{\rho \le (\sqrt(\pi)*\alpha)^(-d)}. The exact simulation of the Gaussian DPPs involves the eigendecomposition of the kernel K
      with respect to the Fourier basis, i.e. an infinite serie with no analytical solution. A truncation is therefore required.
      The above parameter \code{'k'} corresponds to this truncation : the serie will be computed only on the lattice \eqn{\{-k,...,k\}^d}.

    \code{"L1E"} model stands for Exponential DPP, for which \code{param} must have the following form
     \code{param = list(rho, alpha, k)}. An Exponential DPP is a DPP with kernel
         \deqn{K(x,y) = \rho exp(-|(x-y)/\alpha|)}.
        where \eqn{\rho} is the intensity of the considered DPP and \eqn{alpha} a scale parameter which must satisfy
         \eqn{\rho \le (2*\alpha)^{-d}}. Here, \eqn{|x|} denotes the \eqn{L^1}-norm of the vector \eqn{x}. As Gaussian DPP, truncation
         of the eigen decomposition of the kernel is required. The above parameter \code{'k'} corresponds to this truncation : the serie will be computed only on the lattice \eqn{\{-k,...,k\}^d}.

    \code{"D"} model stands for Dirichlet DPP, for which \code{param} must have the following form
     \code{param = list(N)}. A one-dimensionnal Dirichlet kernel is given by
     \deqn{C(x,y;n) = \sum exp(2i\pi k(x-y))}
     where the index \eqn{k} takes value in \eqn{\{0,...,n-1\}}.
     A Dirichlet DPP is a DPP where its kernel \eqn{K} can be expressed as the product of \code{d} one-dimensionnal Dirichlet kernels
     \deqn{K(x,y)=\prod C(x_i,y_i;n_i)}
     The vector \code{N} contains all the \eqn{n_i} in the above equation. A realization of a DPP with such kernel will contain a deterministic number of points, equal to \eqn{\prod n_i}. Unlike the previous models, no truncation is needed in that case.

     \code{"Eig"} model is a specific case in which the user directly specifies the choosen eigenvalues (with respect to the Fourier basis) of the kernel.
     \code{param} must have the following form \code{param = list(eigen, index)}. \code{eigen} is the vector of eigenvalues (which must be non-greater than 1), and \code{index} is a matrix such that \code{eigen[i]} is located at \code{index[i,]} in \eqn{Z^d}.
 }

\value{
  If \code{nsim = 1}, returns a \code{"ppp"} (or \code{"ppx"}) class object if \code{with.kernel = FALSE}, else a list with 2 components:
    \item{PP}{ a \code{"ppp"} (or \code{"ppx"}) class object}
    \item{K}{ a symmetric matrix which represents the evaluation of the kernel for each couple of points in the object \code{pp}. The \code{(i,j)} element in K is the kernel evaluated at the \code{i}-th and \code{j}-th points in \code{pp}.}

  If \code{nsim > 1}, returns a list with size \code{nsim}, where each component is described as above.
}

\seealso{computeKernel, computePCF}

\references{

  F. Lavancier, J. Moller and E. Rubak: Determinantal point process models and statistical inference: Extended version.
  J. Roy. Stat. Soc. B 77, 853–877 (2015).

  A. Mazoyer, J.-F. Coeurjolly and P.-O. Amblard: Projection of determinantal point processes.
  to appear in Spatial Statistics (2020).

}

\examples{
  ### 2-dimensional samples with intensity 100
  rho <- 100; d <- 2
  ## Gaussian DPP
  alphaG <- rho^(-1/d)/sqrt(pi)*0.9  ## ensures that the model is well defined
  dppG <- rdppC(d = d, param = list(rho = rho, alpha = alphaG, k = 20), model = "G")

  ## Exponential DPP
  alphaL1E <- rho^(-1/d)/2*0.9  ## ensures that the model is well defined
  dppL1E <- rdppC(d = d, param = list(rho = rho, alpha = alphaL1E, k = 50), model = "L1E")

  ## Dirichlet DPP
  N <- c(10, 10)
  dppD <- rdppC(d = d, param = list(N = N), model = "D")

  ## 'Eig' model
  eigen <- runif(100)
  index <- as.matrix(expand.grid(rep(list(0:9), d)))
  dppEig <- rdppC(d = d, param = list(eigen = eigen, index = index), model = "Eig")


  ## Compute the kernel
  dppG <- rdppC(d = d, param = list(rho, alphaG, 20), model = "G", with.kernel = TRUE)
  dppG$PP  ### Simulated realisation
  dppG$K  ### Evaluation of the kernel for each couple of points in dppG$PP

  ### Influence of the parameter 'k'
set.seed(1)
dppG5 <- rdppC(d = d, param = list(rho, alphaG, 5), model = "G")
set.seed(1)
dppG20 <- rdppC(d = d, param = list(rho, alphaG, 20), model = "G")
c(dppG5$n, dppG20$n)

}
