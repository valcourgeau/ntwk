
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ntwk

<!-- badges: start -->

[![codecov](https://codecov.io/gh/valcourgeau/ntwk/branch/main/graph/badge.svg?token=4MEZZBRBML)](https://codecov.io/gh/valcourgeau/ntwk)
[![R-CMD-check](https://github.com/valcourgeau/ntwk/workflows/R-CMD-check/badge.svg)](https://github.com/valcourgeau/ntwk/actions)
<!-- badges: end -->

The goal of ntwk is to provide functions for the statistical modelling
of network time series, especially with irregularly-spaced data and a
Lévy-type driving noise.

## Installation

*SOON* You can install the released version of ntwk from
[CRAN](https://CRAN.R-project.org) with:

``` r
# install.packages("ntwk")
```

And the development version from
[GitHub](https://github.com/valcourgeau/ntwk) with:

``` r
# install.packages("devtools")
devtools::install_github("valcourgeau/ntwk")
```

## Example

This is a basic example of the GrOU 2-parameter MLE after constructing a
sample path driven by a Brownian motion:

``` r
library(ntwk)
set.seed(1)
n <- 10000
d <- 5
theta_1 <- 1
theta_2 <- 2
mesh_size <- 0.01

# Generate a graph adjacency matrix
adj <- polymer_network(d = d, theta_1 = theta_1, theta_2 = theta_2)

# Generate a path
times <- seq(from = 0, by = mesh_size, length.out = n) # time grid
noise <- matrix(rnorm(d * n, mean = 0, sd = sqrt(mesh_size)), ncol = d)
path <- construct_path(
  nw_topo = adj, noise = noise, y_init = rep(0, d), delta_time = mesh_size
)

# Fits the Theta-GrOU (2 parameters) process
fit <- grou_mle(times = times, data = path, mode = "network")
print(fit) # should be c(theta_1, theta_2) = c(1, 2)
#> [1] 1.077150 2.137579
```
