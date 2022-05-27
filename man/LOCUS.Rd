\name{LOCUS}
\alias{LOCUS}
\title{
LOCUS: Low-rank decomposition of brain connectivity matrices with universal sparsity
}
\description{
This function is to conduct the LOCUS approach for decomposing brain connectivity data into subnetworks. 
}
\usage{
LOCUS(Y, q, V, MaxIteration=100, penalty="SCAD", phi=0.9, approximation=T, preprocess=T, espli1=0.001, espli2=0.001, rho=0.95, silent=F, demean=T)
}

\arguments{
  \item{Y}{Group-level connectivity data from N subjects, which is of dimension N x p, where p is number of edges. Each row of Y represents a subject's vectorized connectivity matrix by \code{Ltrans} function. }
  \item{q}{Number of ICs/subnetworks to extract.}
  \item{V}{Number of nodes in the network. Note: p should be equal to V(V-1)/2. }
  \item{MaxIteration}{Maximum number of iteractions.}
  \item{penalty}{The penalization approach for universal sparsity, which can be \code{NULL}, \code{SCAD}, \code{L1}, \code{Hardthreshold}. }
  \item{phi}{Tuning parameter for uniform sparse penalty.}
  \item{approximation}{Whether to use an approximated algorithm to speed up the algorithm.}
  \item{preprocess}{Whether to preprocess the data, which reduces the data dimension to \code{q} and whiten the data.}
  \item{espli1}{Toleration for convergence on mixing coefficient matrix, i.e. A.}
  \item{espli2}{Toleration for convergence on latent sources, i.e. S.}
  \item{rho}{Tuning parameter for selecting number of ranks in each subnetwork's decomposition. }
  \item{silent}{Whether to print intermediate steps.}
  \item{demean}{Whether to demean Y.}
}

\details{
  This is the main function for LOCUS decomposition of brain connectivity matrices, which is to minimize the following objective function:
    \deqn{\sum_{i=1}^N\|Y_i - \sum_{l=1}^q a_{il} S_l\|_2^2 + \phi \sum_{l=1}^q\|S_l\|_*,}
  where \eqn{Y_i} is the transpose of \eqn{i}th row in \eqn{Y}, \eqn{S_l = L(X_l D_l X_l')} represents the \eqn{l}th vectorized latent source/subnetwork with low-rank decomposition, \eqn{L} is \code{Ltrans} function, \eqn{\|\cdot\|_*} represents the penalty which can either be L1, SCAD or hard-thresholded. Further details can be found in the LOCUS paper.  
}
\value{
  An R list from Locus containing the following terms:
   \item{Conver}{Whether the algorithm is converaged.}
   \item{A}{Mixing matrix \eqn{\{a_{il}\}} of dimension N by q.}
   \item{S}{Subnetworks of dimension q by p, where each row represents a vectorized subnetwork based on \code{Ltrans} function.}
   \item{theta}{A list of length q, where \code{theta[[i]]} contains the symmetric low-rank decomposition of \code{i}th subnetwork. }
   }


\references{
Wang, Y. and Guo, Y. (2020+).   \emph{LOCUS: A novel signal decomposition method for brain network connectivity matrices using low-rank structure with universal sparsity. } In preparation. 
}

\author{Yikai Wang\cr
  Maintainer: Yikai Wang \email{johnzon.wyk@gmail.com} }
  
\examples{
## Simulated the data to use. 
V = 50
S1 = S2 = S3 = matrix(0,ncol = V,nrow = V)
S1[5:20,5:20] = 4;S1[23:37,23:37] = 3;S1[40:48,40:48] = 3
S2[15:20,] = -3;S2[,15:20] = -3
S3[15:25,36:45] = 3; S3[36:45,15:25] = 3
Struth = rbind(Ltrans(S1,FALSE) , Ltrans(S2,FALSE), Ltrans(S3,FALSE))
set.seed(100)
Atruth = matrix(rnorm(100*3),nrow=100,ncol=3)
Residual = matrix(rnorm(100*dim(Struth)[2]),nrow=100)
Yraw = Atruth\%*\%Struth + Residual

##### Run Locus on the data ##### 
Locus_result = LOCUS(Yraw,3,V)

par(mfrow=c(2,3))
for(i in 1:dim(Struth)[1]){image(Ltrinv(Struth[i,],V,FALSE))}
for(i in 1:dim(Locus_result$S)[1]){image(Ltrinv(Locus_result$S[i,],V,FALSE))}

}