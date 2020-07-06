[![](https://www.r-pkg.org/badges/version/svmadmm?color=green)](https://cran.r-project.org/package=svmadmm)
[![](http://cranlogs.r-pkg.org/badges/grand-total/svmadmm?color=red)](https://cran.r-project.org/package=svmadmm)
[![](http://cranlogs.r-pkg.org/badges/last-month/svmadmm?color=blue)](https://cran.r-project.org/package=svmadmm)
[![](http://cranlogs.r-pkg.org/badges/last-week/svmadmm?color=green)](https://cran.r-project.org/package=svmadmm)

# svmadmm
Linear/Nonlinear SVM Classification Solver Based on ADMM and IADMM Algorithms

## Description

`svm.admm` is a simple function for solving large-scale regularized linear/nonlinear
 classification by using ADMM and IADMM algorithms. This function provides
 linear L2-regularized primal classification (both ADMM and IADMM are available),
 kernel L2-regularized dual classification (IADMM) as well as L1-regularized primal
 classification (both ADMM and IADMM are available). The training of the models perform well
 practice.

## Installation

You can install `svmadmm` by 
```r
devtools::install_github("statmlben/svmadmm")
```

## Usage

```r
svm.admm(x.tr, y.tr, type = 0, kernel = "radial", sigma = 1/ncol(x.tr),
  degree = 1, scale = 1, offset = 1, algo = "iadmm", lambda = 1,
  rho = 1, eps = 0.01)
```


## Arguments

Argument      |Description
------------- |----------------
```x.tr```     |     a n*p data matrix. Each row stands for an example (sample, point) and each column stands for a dimension (feature, variable).
```y.tr```     |     a n-length vector. The values correspond to class labels.
```type```     |     `svmadmm` can provide 3 types of linear/kernel models. The default value for `type` is 0. See details below. Valid options are:   

*  0 -- L2-regularized kernel svm (dual)  

*  1 -- L2-regularized linear svm (primal)  

*  2 -- L1-regularized linear svm (primal) 
```kernel```     |     the kernel used in training and predicting when `type` is 0. The default value for `kernel` is radial. See details below. Valid options are:   

*   `radial` -- The Gaussian RBF kernel k(x,x') = exp(-sigma \|x - x'\|^2)  

*   `linear` -- The Linear kernel k(x,x') = <x, x'>  

*   `polynomial` -- The Polynomial kernel k(x,x') = (scale <x, x'> + offset)^ degree  
```sigma```     |     The inverse kernel width used by the Gaussian.
```degree```     |     The degree of the polynomial kernel function. This has to be an positive integer.
```scale```     |     The scaling parameter of the polynomial kernel is a convenient way of normalizing patterns without the need to modify the data itself.
```offset```     |     The offset used in a polynomial kernel.
```algo```     |     the algorithm to solve the problem w.r.t. `type` .
```lambda```     |     regularization constant (default: 1). Rules the trade-off between regularization and correct classification on `data` .
```rho```     |     regularization constant (default: 1).
```eps```     |     epsilon in the termination condition.

## Details


 `svmadmm` internally computing kernel matrix when `type` is 0, which is based by the package kernlab .


## Examples

```r 
 library(svmadmm)
 n = 100
 p = 10
 x = matrix(runif(2 * n * p, -1, 1), nrow = 2 * n)
 y = sign(x[, 1])
 y.ind = sample(1 : (2 * n), n / 10, replace = FALSE)
 y[y.ind] = - y[y.ind]
 x.tr = x[1 : n, ]
 y.tr = y[1 : n]
 x.te = x[-(1 : n), ]
 y.te = y[-(1 : n)]
 model = svm.admm(x.tr, y.tr)
 fit = svm.predict(x.te, model)
``` 
