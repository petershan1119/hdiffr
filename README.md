
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hdiffr

<!-- badges: start -->

<!-- badges: end -->

**hdiffr** performs an event history diffusion model that is modified in
order to account for heterogeneity in the diffusion process. In addition
to modeling characteristics that influnece the direct propensity of a
focal actor to adopt, **hdiffr** can also model characteristics of the
focal actor that make it more susceptible to the prior adoptions of
other actors, as well as characteristics of prior-adopting actors that
make them more or less influential on others, and characteristics of the
social structure among actors that influence the transmission of
influence between actors.

## Author

[Sang Won Han](https://sociology.columbia.edu/content/sang-won-han),
Ph.D. candidate in [Sociology](https://sociology.columbia.edu/) at
[Columbia University](https://www.columbia.edu)

## Reference

[Strang, David and Nancy Brandon Tuma. 1993. “Spatial and Temporal
Heterogeneity in Diffusion.” *American Journal of Sociology* 99(3):
614-639.](https://www.journals.uchicago.edu/doi/abs/10.1086/230318)

## Installation

The experimental version can be obtained via:

``` r
devtools::install_github("petershan1119/hdiffr")
```

## Example

This is a basic example which performs heterogeneous diffusion model
with multispell data:

``` r
library(hdiffr)
data(exampleData)
result <- hdiffr(data = exampleData, xvars = c('lnsale', 'roa', 'vote'), vvars = 'activistgrp', wvars = 'lnsale', 
                 zgvars = 'sic', idvar = 'gvkey', multispell = 1, vintercept = 1, hrno = 1)
#> Warning in hdiffr(data = exampleData, xvars = c("lnsale", "roa", "vote"), :
#> Warning: More than one observation per ID detected.
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

If you want to see the results, use `summary`:

``` r
summary(result)
#> 
#> Call:
#> survival::survreg(formula = as.formula(formula), data = data_wf_gi, 
#>     dist = dist)
#>                 Value Std. Error     z       p
#> (Intercept)    5.6198     1.3413  4.19 2.8e-05
#> lnsale        -0.3842     0.1492 -2.58   0.010
#> roa           -2.2903     2.1201 -1.08   0.280
#> vote          -0.0241     0.0230 -1.05   0.296
#> v_activistgrp -0.1493     0.0571 -2.61   0.009
#> v_intercept    0.1355     0.0826  1.64   0.101
#> w_lnsale      -0.0417     0.0377 -1.11   0.268
#> zg_sic        -0.3177     0.1340 -2.37   0.018
#> 
#> Scale fixed at 1 
#> 
#> Exponential distribution
#> Loglik(model)= -225.4   Loglik(intercept only)= -245.3
#>  Chisq= 39.72 on 7 degrees of freedom, p= 1.4e-06 
#> Number of Newton-Raphson Iterations: 6 
#> n= 575
```
