
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Raquifer

<!-- badges: start -->

<!-- badges: end -->

`Raquifer` estimates the cumulative water influx into hydrocarbon
reservoirs using un-steady and pseudo-steady state modeling approaches.
It generates a data frame of cumulative water influx over time for edge-
and bottom-drive aquifers.

## Installation

You can install the released version of Raquifer from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Raquifer")
```

## Example

``` r
library(Raquifer)

aqu_time <- aquifer_time(x = c(0,0.368,2.439,4.957,7.732,11.926,18.126,30.044) * 365, unit = "day")

parameters <- aquifer_param(input_unit = "Field", output_unit = "Field", model = "uss", 
                            flow_type = "radial", water_drive = "edge", phi = 0.27, perm_h = 64.2, 
                            h_a = 20, r_a = 5 * 14892, r_R = 14892, tetha = 180,
                            mu_water = 0.485, c_water = 3.88e-6, c_rock = 2e-6, 
                            pressure = c(1640,1600,1400,1200,1000,800,600,400))

pred_veh <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)
```
