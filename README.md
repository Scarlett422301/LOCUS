# R/`LOCUS`

> Low-rank decomposition of brain connectivity matrices with universal sparsity

**Author:** [Yikai Wang](https://sites.google.com/view/yikaiw/), Jialu Ran, Ying Guo

-----

## Description

`LOCUS` is an R package that conducts the low-rank decomposition of brain connectivity matrices 
with universal sparsity to reveal sub brain networks.

-----

## Installation

-----

## Usage

Below is an illustration of the the main function on simulated data.

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
