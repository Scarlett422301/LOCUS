# LOCUS: Low-rank Decomposition of Brain Connectivity Matrices with Universal Sparsity Method

LOCUS (Low-rank Decomposition of Brain Connectivity Matrices with Universal Sparsity Method) is an innovative blind source separation method that incorporates regularization and low-rank structure to investigate brain connectivity.

## Installation

We assume you are running R 4.1.0 or newer. There is no guarantee for backward or forward comparability. Please raise the issue on GitHub if something breaks.

You can install LOCUS from CRAN with:

``` r
install.packages("LOCUS")
```

If you want to install LOCUS from github,

the following R packages are required:

-   MASS (\>= 7.3.60)
-   ica (\>= 1.0.3)
-   far (\>= 0.6.6)
-   devtools

You can install them by running this code:

``` r
if(!require(c("MASS","ica","far", "devtools"))){
    install.packages(c("MASS","ica","far","devtools"))
}
```

Then you can install LOCUS from github with:

``` r
library(devtools)
install_github("Emory-CBIS/LOCUS")
# Load the package
library(LOCUS)
```

## Tutorial

The `LOCUS()` is the main function of our LOCUS algorithm. The `LOCUS_BIC_selection()` selects the tuning parameter phi and rho based on our proposed BIC-like criterion.

### Explanation of Arguments

#### 1. LOCUS function

```         
LOCUS(Y, q, V, MaxIteration=100, penalty="SCAD", phi = 0.9, approximation=TRUE, 
preprocess=TRUE, espli1=0.001, espli2=0.001, rho=0.95, demean = TRUE, silent=FALSE)
```

-   `Y`: Group-level connectivity data from N subjects, which is of dimension N x p, where p is number of edges. Each row of Y represents a subject's vectorized connectivity matrix by `Ltrans` function.
-   `q`: Number of ICs/subnetworks to extract.
-   `V`: Number of nodes in the network. Note: p should be equal to V(V-1)/2.
-   `MaxIteration`: Maximum number of iterations. The default number is 100.
-   `penalty`: The penalization approach for uniform sparsity, which can be "NULL"","SCAD", "L1" and "Hardthreshold". Defaults to "SCAD".
-   `phi`: The tuning parameter $\phi$ for uniform sparse penalty. The default is 0.9.
-   `approximation`: Whether to use an approximated algorithm to speed up the algorithm. The default is "TRUE". It is suggested to be used.
-   `preprocess`: Whether to preprocess the data, which reduces the data dimension to q and whiten the data. The default is "TRUE".
-   `espli1`: Toleration for convergence on mixing coefficient matrix, i.e. A. The default is 0.001.
-   `espli2`: Toleration for convergence on latent sources, i.e. S. The default is 0.001.
-   `rho`: The tuning parameter $\rho$ for selecting number of ranks in each subnetwork's decomposition. The default is 0.95.
-   `demean`: Whether to subtract the mean from each column of Y. The default is "TRUE".
-   `silent`: Whether to print intermediate steps. The default is "FALSE".

#### 2. LOCUS_BIC_selection function

```         
LOCUS_BIC_selection(Y, q, V, MaxIteration=50, penalty="SCAD", 
phi_grid_search=seq(0.2, 1, 0.2), rho_grid_search=c(0.95), 
espli1=0.001, espli2=0.001, save_LOCUS_output=TRUE, 
preprocess=TRUE, demean = TRUE)
```

-   `Y`: Group-level connectivity data from N subjects, which is of dimension N x p, where p is number of edges. Each row of Y represents a subject's vectorized connectivity matrix by `Ltrans` function.
-   `q`: Number of ICs/subnetworks to extract.
-   `V`: Number of nodes in the network. Note: p should be equal to V(V-1)/2.
-   `MaxIteration`: Maximum number of iterations. The default number is 100.
-   `penalty`: The penalization approach for uniform sparsity, which can be "NULL"","SCAD", "L1" and "Hardthreshold". Defaults to "SCAD".
-   `phi_grid_search`: Grid search candidates for tuning parameter of uniform sparse penalty.
-   `espli1`: Toleration for convergence on mixing coefficient matrix, i.e. A. The default is 0.001.
-   `espli2`: Toleration for convergence on latent sources, i.e. S. The default is 0.001.
-   `rho_grid_search`: Grid search candidates for tuning parameter for selecting number of ranks in each subnetwork's decomposition.
-   `save_LOCUS_output`: Whether to save LOCUS output from each grid search. The default is TRUE.
-   `preprocess`: Whether to preprocess the data, which reduces the data dimension to q and whiten the data. The default is TRUE.
-   `demean`: Whether to subtract the mean from each column of Y. The default is TRUE.

### Explanation of Output

#### 1. LOCUS function

The output will be a list with 4 components as such: - `Conver`: Whether the algorithm is converaged. - `A`: Mixing matrix ${\{a_{il}\}$ of dimension N x q. - `S`: Subnetworks of dimension q x p, where each row represents a vectorized subnetwork based on `Ltrans` function. - `theta`: A list of length q, where `theta[[i]]` contains the symmetric low-rank decomposition of \code{i}th subnetwork.

#### 2. LOCUS_BIC_selection function

Function outputs a list including the following: - `bic_tab`: BIC values per phi and rho. - `LOCUS_results`: LOCUS output, if save_LOCUS_output is TRUE.

## Example

Load the example data:

``` r
## Simulated the data to use
V = 50
S1 = S2 = S3 = matrix(0,ncol = V,nrow = V)
S1[5:20,5:20] = 4;S1[23:37,23:37] = 3;S1[40:48,40:48] = 3
S2[15:20,] = -3;S2[,15:20] = -3
S3[15:25,36:45] = 3; S3[36:45,15:25] = 3
Struth = rbind(Ltrans(S1,FALSE) , Ltrans(S2,FALSE), Ltrans(S3,FALSE))
set.seed(100)
Atruth = matrix(rnorm(100*3),nrow=100,ncol=3)
Residual = matrix(rnorm(100*dim(Struth)[2]),nrow=100)
Yraw = Atruth%*%Struth + Residual

## Run Locus on the data 
Locus_result = LOCUS(Yraw,3,V)

## Visualize the result
par(mfrow=c(2,3))
for(i in 1:dim(Struth)[1]){image(Ltrinv(Struth[i,],V,FALSE))}
for(i in 1:dim(Locus_result$S)[1]){image(Ltrinv(Locus_result$S[i,],V,FALSE))}
```

## References

<div id="refs" class="references">
