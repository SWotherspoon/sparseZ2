
# sparseZ2

<!-- badges: start -->
<!-- badges: end -->

The sparseZ2 package provides facilities for solving sparse linear
systems $A x = b$ over the field $\mathbb{Z}_{2}$.

## Installation
The package is easily installed from GitHub, using the devtools package. 

```r
devtools::install_github("SWotherspoon/sparseZ2")
```

If you don't have `devtools` installed already, install it first. 

```r
install.packages("devtools")
```

(sparseZ2 otherwise does not need devtools for normal use.)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sparseZ2)
# Create a sparse matrix
A <- sparseZ2(4, 4)
A[1,c(1,3,4)] <- 1
A[2,c(2,3)] <- 1
A[3,c(1,2,4)] <- 1
A[4,c(1,4)] <- 1

b <- c(1, 0, 1, 1)

# Perform Gaussian elimination
sys <- gaussJordanZ2(A, b)
x <- solution(sys$A,sys$b)
x
```

