\name{computeKernel}
\alias{computeKernel}
\encoding{UTF-8}

\title{Computation of kernel of a point configuration.}

\description{
Evaluates a kernel on a point configuration according to a given kernel class.
}

\usage{

  computeKernel(pp, param = NULL,
        type = c("prod", "sum"),
        model = c("D", "L1E", "G", "Eig"),
        domain = NULL)
        )
}

\arguments{
 \item{pp}{ ppp/ppx object or matrix of points. Can also be a list of matrices.}
 \item{param}{ list of parameters. See Details.}
 \item{type}{ Construction of the DPP. See Details. }
 \item{model}{ Class of the kernel. See Details.}
 \item{domain}{ domain of observation. Must be of class \code{"ownin"}, \code{"box3"} or \code{"boxx"}. Ignored if \code{pp} is a ppp/ppx object.}
}
\details{

Two kinds of kernels can be considered. When \code{type} is set to \code{'prod'}, the kernel is defined as
\deqn{K(x,y) = \sum K_j(x_j,y_j)}
If \code{type} is set to \code{'prod'}, the kernel is defined as
\deqn{K(x,y) = \prod K_j(x_j,y_j)}
where the index \eqn{j} takes value in \eqn{{1,...,d}}.
The \eqn{K_j}'s are 1-dimensional kernels defining according to \code{model}.

The list \code{param} depends also on \code{model} (the order of appearance of the components must be respected, or the user must set the names as below)

  \code{"G"} model stands for Gaussian DPP, for which \code{param} must have the following form
   \code{param = list(rho, alpha, k)}. The parameter \eqn{\rho} is the intensity of the considered DPP and \eqn{alpha} a scale parameter. If \code{type = 'prod'} then \eqn{K_j} is defined by
    \deqn{K_j(x,y) = \rho^(1/d) exp(-((x-y)/\alpha)^2)}.
    where \eqn{\alpha} must satisfy \eqn{\rho \le (\sqrt(\pi)*\alpha)^(-d)}. If \code{type = 'sum'} then \eqn{K_j} is defined by
     \deqn{K_j(x,y) = \rho/d exp(-((x-y)/\alpha)^2)}.
     where \eqn{\alpha} must satisfy \eqn{\rho \le (\sqrt(\pi)*\alpha)^(-1)}. The exact simulation of the Gaussian DPPs involves the eigendecomposition of the kernel \eqn{K} with respect to the Fourier basis, i.e. an infinite serie with no analytical solution. A truncation is therefore required.
    The above parameter \code{'k'} corresponds to this truncation : the serie will be computed only on the lattice \eqn{{-k,...,k}^d}.

  \code{"L1E"} model stands for Exponential DPP, for which \code{param} must have the following form
   \code{param = list(rho, alpha, k)}. The parameter \eqn{\rho} is the intensity of the considered DPP and \eqn{alpha} a scale parameter. If \code{type = 'prod'} then \eqn{K_j} is defined by
       \deqn{K_j(x,y) = \rho^(1/d) exp(-|(x-y)/\alpha|)}.
       where \eqn{alpha} must satisfy \eqn{\rho \le (2*\alpha)^(-d)}. If \code{type = 'sum'} then \eqn{K_j} is defined by
           \deqn{K_j(x,y) = \rho/d exp(-|(x-y)/\alpha|)}.
           where \eqn{alpha} must satisfy \eqn{\rho \le (2*\alpha)^(-1)}. As Gaussian DPP, truncation
       of the eigen decomposition of the kernel is required. The above parameter \code{'k'} corresponds to this truncation : the serie will be computed only on the lattice \eqn{{-k,...,k}^d}.

  \code{"D"} model stands for Dirichlet DPP, for which \code{param} must have the following form
   \code{param = list(N)}. If \code{type = 'prod'} then \eqn{K_j} is defined by
   \deqn{K_j(x,y; n_j) = \sum exp(2i\pi k(x-y))}
   where the index \eqn{k} takes value in \eqn{{0,...,n_j-1}}. If \code{type = 'sum'} then \eqn{K_j} is defined by
   \deqn{K_j(x,y; n_j) = \sum exp(2i\pi k(x-y))}
   where the index \eqn{k} takes value in \eqn{{1,...,n_j}}.

   The vector \code{N} contains all the \eqn{n_j} in the above equations. A realization of a DPP with such kernel will contain a deterministic number of points, equal to \eqn{\prod n_j} if \code{type = 'prod'}, \eqn{\sum n_j} if \code{type = 'sum'}. Unlike the previous models, no truncation is needed in that case.

   \code{"Eig"} model is a specific case in which the user directly specifies the choosen eigenvalues (with respect to the Fourier basis) of the kernel.
   \code{param} must have the following form \code{param = list(eigen, index)}. \code{eigen} is the vector of eigenvalues (which must be non-greater than 1), and \code{index} is a matrix such that \code{eigen[i]} is located at \code{index[i,]} in \eqn{Z^d}.
}

\value{
  If \code{pp} is a matrix or a ppp/ppx object, returns a symmetric matrix which represents the evaluation of the kernel for each couple of points in the object \code{pp}. The \code{(i,j)} element in of the matrix is the kernel evaluated at the \code{i}-th and \code{j}-th points in \code{pp}.

  If \code{pp} is a list of matrices, returns a a list of symmetric matrices as above.


}

\seealso{rdppC,computePCF}

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
KG5 <- computeKernel(pp, param = list(rho, alphaG, 5), model = "G")
KG20 <- computeKernel(pp, param = list(rho, alphaG, 20), model = "G")
KGExact <- computeKernel(pp, param = list(rho, alphaG, +Inf), model = "G")
}
