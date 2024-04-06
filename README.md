<img src="fig/LOCUS.png" width="900" align="center"/>

# LOCUS: Low-rank Decomposition of Brain Connectivity Matrices with Universal Sparsity Method

LOCUS (Low-rank Decomposition of Brain Connectivity Matrices with Universal Sparsity Method) is an innovative blind source separation method that incorporates regularization and low-rank structure to investigate brain connectivity. This package is as an implementation of the LOCUS method proposed by [Wang and Guo, 2023](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-17/issue-2/LOCUS--A-regularized-blind-source-separation-method-with-low/10.1214/22-AOAS1670.full).

## Installation

You can install LOCUS from CRAN with:

``` r
install.packages("LOCUS")
```

## Method

LOCUS is a fully data-driven blind signal separation method for decomposing imaging-based brain connectivity data to reveal underlying neural circuits. Specifically, LOCUS decomposes subjects' connectivity data, $\mathbf Y$, into a linear combination of latent connectivity traits or source signals, $\{\mathbf S_{\ell}\}$, weighted by mixing coefficients $\{\mathbf{a}_{\ell}\}_{\ell=1}^{q}$, i.e. $$\mathbf{Y} = \sum_{\ell=1}^{q} {\mathbf{a}_\ell} \mathbf{S}_{\ell} + \text{error}.$$

Here, each of the connectivity source signals $\mathbf S_\ell$ represents an underlying neural circuit and the mixing coefficients $\{\mathbf{a}_{\ell}\}$ represent subject-specific loadings on the trait. The method has the following highlights:

-   LOCUS models the source signals $\mathbf S_\ell$ using a low-rank structure $\mathbf X_{\ell}\mathbf D_{\ell}\mathbf X_{\ell}'$ where $\mathbf D_{\ell}$ is a diagonal matrix. This is well motivated by the observation that brain connectivity matrices often include block-diagonal or banded structure that can be efficiently captured with a low-rank factorization. The low-rank structure leads to a significant reduction in the number of parameters, hence improving accuracy and reliability in the recovery of underlying connectivity traits. LOCUS also incorporates an adaptive rank selection approach to flexibly choose the rank for each of the connectivity source signals. The source-specific rank selection allows better accommodation of varying spatial patterns in the neural circuits across the brain.

-   The subject-specific trait loadings $\{\mathbf{a}_{\ell}\}$, generated by LOCUS, quantify the prominence or presence of each trait in a subject's connectivity. These subject-specific trait loadings capture between-subject heterogeneity in each of the connectivity traits and also allow identifying which traits are associated with clinical or demographical variables.

-   To reduce spurious scientific findings, LOCUS proposes a more intuitive way of sparsity regularization by directly penalizing on latent sources $\mathbf S_{\ell}$ via the low-rank structure $\mathbf X_{\ell}\mathbf D_{\ell}\mathbf X_{\ell}'$. This is equivalent as angle-based penalization.

-   The learning of LOCUS is formulated as a non-convex optimization problem. The optimization function as below has a block multi-convex structure. LOCUS incorporates an efficient node-rotation algorithm with closed-form solutions at each iteration for estimating the parameters.

$$\text{min} \sum_{i=1}^N\| \mathbf{y}_i - \sum_{\ell=1}^q a_{i\ell} \mathcal{L}(\mathbf{X}_\ell\mathbf{D}_\ell\mathbf{X}_{\ell}') \|_2^2 + \phi\sum_{\ell=1}^q \sum_{u<v} |\mathbf{x}_{\ell}(u)'\mathbf{D}_\ell\mathbf{x}_{\ell}(v)|$$

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

The output will be a list with 4 components as such:

-   `Conver`: Whether the algorithm is converaged.
-   `A`: Mixing matrix $\{a_{il}\}$ of dimension N x q.
-   `S`: Subnetworks of dimension q x p, where each row represents a vectorized subnetwork based on `Ltrans` function.
-   `theta`: A list of length q, where `theta[[i]]` contains the symmetric low-rank decomposition of ith subnetwork.

#### 2. LOCUS_BIC_selection function

Function outputs a list including the following:

-   `bic_tab`: BIC values per phi and rho.
-   `LOCUS_results`: LOCUS output, if save_LOCUS_output is TRUE.

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

If you find this package useful, please cite:

```         
@article{wang2023locus,
  title={LOCUS: A regularized blind source separation method with low-rank structure for investigating brain connectivity},
  author={Wang, Yikai and Guo, Ying},
  journal={The Annals of Applied Statistics},
  volume={17},
  number={2},
  pages={1307--1332},
  year={2023},
  publisher={Institute of Mathematical Statistics}
}
```
