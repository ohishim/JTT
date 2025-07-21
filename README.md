
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JTT (v0.1.0)

<!-- badges: start -->

<!-- badges: end -->

This package implements consistent clustering for group-wise linear
regression with graph structure.

**cite this package:**  
Ohishi, M. (2025). JTT: Consistent clustering for group-wise linear
regression. R package version 0.1.0. <https://github.com/ohishim/JTT>

## Installation

You can install the R package `JTT` like so:

``` r
# install.packages("devtools")
devtools::install_github("ohishim/JTT")
library(JTT)
```

## Example

You can get simulation data by `genData()` and implement JTT method by
`JTT()`. The simulation data is generated as

``` r
n0 <- 50; p <- 20
Data <- genData(n0, p, SNR=3)
```

You can obtain true clusters as

``` r
Data$cluster
#> $c1
#> [1] 1 2 3
#> 
#> $c2
#> [1]  4  5  6  9 10
#> 
#> $c3
#> [1] 7 8
```

Then, JTT method is implemented as

``` r
y <- Data$y; X <- Data$X; group <- Data$group; adj <- Data$adj
res <- JTT(y, X, group, adj)
```

Estimated clusters are

``` r
res$cluster
#> $c1
#> [1] 1 2 3
#> 
#> $c2
#> [1]  4  5  6  9 10
#> 
#> $c3
#> [1] 7 8
```
