\name{calculate_bic}
\alias{calculate_bic}
\title{
BIC calculation}
\description{
This function is to calculate BIC score for LOCUS output following the BIC formula defined in Wang, Y. and Guo, Y. (2022+). 
}
\usage{
calculate_bic(Y, A, S)
}

\arguments{
  \item{Y}{Group-level connectivity data from N subjects, which is of dimension N x p, where p is number of edges. Each row of Y represents a subject's vectorized connectivity matrix by \code{Ltrans} function. }
  \item{A}{Mixing Matrix.}
  \item{S}{Latent Source.}
}

\details{
  In Wang, Y. and Guo, Y. (2022+), the BIC-type criterion is proposed to select parameters \eqn{\phi, \rho}. The criterion balances between model fitting and model sparsity. 
}
\value{A BIC value will be returned.}