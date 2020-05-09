\name{computePCF}
\alias{computePCF}
\encoding{UTF-8}

\title{Computation of the pair correlation function of a point configuration.}

\description{
Evaluates a pair correlation function on a point configuration according to a given kernel class.
}

\usage{

  computePCF(pp, param = NULL,
        model = c("D", "L1E", "G", "Eig"),
        domain = NULL)
        )
}

\arguments{
 \item{pp}{ ppp/ppx object or matrix of points. Can also be a list of matrices.}
 \item{param}{ list of parameters. See Details.}
 \item{domain}{ domain of observation. Must be of class \code{"ownin"}, \code{"box3"} or \code{"boxx"}. Ignored if \code{pp} is a ppp/ppx object.}
 \item{model}{ Class of the kernel. See Details.}
}
\details{

  The pair correlation function associated to a kernel \eqn{K} is given by
  \deqn{
    g(x,y) = 1-|K(x,y)|^2/(K(x,x)*K(x,y))
  }

  The list \code{param} depends on \code{model} (the order of appearance of the components must be respected, or the user must set the names as below)

    \code{"G"} model stands for Gaussian kernel, for which \code{param} must have the following form
     \code{param = list(rho, alpha, k)}. A Gaussian kernel is given by
      \deqn{K(x,y) = \rho exp(-||(x-y)/\alpha||^2)}.
     where \eqn{\rho} is the intensity of the considered DPP and \eqn{alpha} a scale parameter which must satisfy
      \eqn{\rho \le (\sqrt(\pi)*\alpha)^(-d)}.
      The exact simulation of the Gaussian DPPs involves the eigendecomposition of the kernel K
      with respect to the Fourier basis, i.e. an infinite serie with no analytical solution. A truncation is therefore required.
      The above parameter \code{'k'} corresponds to this truncation : the serie will be computed only on the lattice \eqn{\{-k,...,k\}^d}.
      In \code{computeKernel} it is possible to set \code{k = +Inf}. In that case, the exact kernel (i.e. givn by the above equation) will be computed.

    \code{"L1E"} model stands for Exponential kernel, for which \code{param} must have the following form
     \code{param = list(rho, alpha, k)}. An Exponential kernel is given by
         \deqn{K(x,y) = \rho exp(-|(x-y)/\alpha|)}.
        where \eqn{\rho} is the intensity of the considered DPP and \eqn{alpha} a scale parameter which must satisfy
         \eqn{\rho \le (2*\alpha)^{-d}}. Here, \eqn{|x|} denotes the \eqn{L^1}-norm of the vector \eqn{x}. As Gaussian DPP, truncation
         of the eigen decomposition of the kernel is required. The above parameter \code{'k'} corresponds to this truncation : the serie will be computed only on the lattice \eqn{\{-k,...,k\}^d}.
         In \code{computeKernel} it is possible to set \code{k = +Inf}. In that case, the exact kernel (i.e. givn by the above equation) will be computed.

    \code{"D"} model stands for Dirichlet kernel, for which \code{param} must have the following form
     \code{param = list(N)}. A one-dimensionnal Dirichlet kernel is given by
     \deqn{C(x,y;n) = \sum exp(2i\pi k(x-y))}
     where the index \eqn{k} takes value in \eqn{\{0,...,n-1\}}.
     A Dirichlet kernel \eqn{K} can be expressed as the product of \code{d} one-dimensionnal Dirichlet kernels
     \deqn{K(x,y)=\prod C(x_i,y_i;n_i)}
     The vector \code{N} contains all the \eqn{n_i} in the above equation. A realization of a DPP with such kernel will contain a deterministic number of points, equal to \eqn{\prod n_i}. Unlike the previous models, no truncation is needed in that case.

     \code{"Eig"} model is a specific case in which the user directly specifies the choosen eigenvalues (with respect to the Fourier basis) of the kernel.
     \code{param} must have the following form \code{param = list(eigen, index)}. \code{eigen} is the vector of eigenvalues (which must be non-greater than 1), and \code{index} is a matrix such that \code{eigen[i]} is located at \code{index[i,]} in \eqn{Z^d}.
 }

\value{
  If \code{pp} is a matrix or a ppp/ppx object, returns a symmetric matrix which represents the evaluation of the pcf for each couple of points in the object \code{pp}. The \code{(i,j)} element in of the matrix is the pcf evaluated at the \code{i}-th and \code{j}-th points in \code{pp}.

  If \code{pp} is a list of matrices, returns a a list of symmetric matrices as above.


}

\seealso{rdppC,computeKernel}

\references{

  F. Lavancier, J. Moller and E. Rubak: Determinantal point process models and statistical inference: Extended version.
  J. Roy. Stat. Soc. B 77, 853–877 (2015).

  A. Mazoyer, J.-F. Coeurjolly and P.-O. Amblard: Projection of determinantal point processes.
  to appear in Spatial Statistics (2020).

}

\examples{
### 2-dimensional Gaussian DPP simulation with intensity 100
rho <- 100; d <- 2
alphaG <- rho^(-1/d)/sqrt(pi)*0.9  ## ensures that the model is well defined
pp <- rdppC(d = d, param = list(rho, alphaG, 20), model = "G")

### Influence of the parameter 'k'
KG5 <- computePCF(pp, param = list(rho, alphaG, 5), model = "G")
KG20 <- computePCF(pp, param = list(rho, alphaG, 20), model = "G")
KGExact <- computePCF(pp, param = list(rho, alphaG, +Inf), model = "G")
}
