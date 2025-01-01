
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMORE

<!-- badges: start -->
<!-- badges: end -->

scMORE is designed to identify disease-relevant regulons by leveraging GWAS summary statistics and multimodal single-cell measurements, where gene expression and chromatin accessibility are either measured for each individual cell or integrated into metacell or clusters to capture both modalities.

## Installation

We recommend installing scMORE via Github using devtools:

``` r
# install.packages("devtools")
devtools::install_github("mayunlong89/scMORE")
```
You can install the old version of scMORE from
[GitHub](https://github.com/mayunlong89/ctDRTF)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scMORE)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
