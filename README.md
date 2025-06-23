
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PsyR

<!-- badges: start -->

<!-- badges: end -->

The goal of PsyR is to …

## Installation

You can install the development version of PsyR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("garner-code/PsyR")
#> Downloading GitHub repo garner-code/PsyR@HEAD
#> stringi      (1.8.4  -> 1.8.7 ) [CRAN]
#> rlang        (1.1.5  -> 1.1.6 ) [CRAN]
#> cli          (3.6.4  -> 3.6.5 ) [CRAN]
#> utf8         (1.2.4  -> 1.2.6 ) [CRAN]
#> Rdpack       (2.6.2  -> 2.6.4 ) [CRAN]
#> pillar       (1.10.1 -> 1.10.2) [CRAN]
#> ggplot2      (3.5.1  -> 3.5.2 ) [CRAN]
#> Deriv        (4.1.6  -> 4.2.0 ) [CRAN]
#> tibble       (3.2.1  -> 3.3.0 ) [CRAN]
#> generics     (0.1.3  -> 0.1.4 ) [CRAN]
#> MatrixModels (0.5-3  -> 0.5-4 ) [CRAN]
#> doBy         (4.6.26 -> 4.6.27) [CRAN]
#> broom        (1.0.7  -> 1.0.8 ) [CRAN]
#> reformulas   (0.4.0  -> 0.4.1 ) [CRAN]
#> nloptr       (2.1.1  -> 2.2.1 ) [CRAN]
#> scales       (1.3.0  -> 1.4.0 ) [CRAN]
#> lme4         (1.1-36 -> 1.1-37) [CRAN]
#> quantreg     (6.00   -> 6.1   ) [CRAN]
#> emmeans      (1.10.7 -> 1.11.1) [CRAN]
#> Installing 19 packages: stringi, rlang, cli, utf8, Rdpack, pillar, ggplot2, Deriv, tibble, generics, MatrixModels, doBy, broom, reformulas, nloptr, scales, lme4, quantreg, emmeans
#> Installing packages into '/tmp/RtmpBkYDYw/temp_libpathc3ed7e1a584a'
#> (as 'lib' is unspecified)
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/tmp/RtmpIbwuBI/remotes102697bae37d2/garner-code-PsyR-6932a1c/DESCRIPTION’ ... OK
#> * preparing ‘PsyR’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘PsyR_0.0.0.9000.tar.gz’
#> Installing package into '/tmp/RtmpBkYDYw/temp_libpathc3ed7e1a584a'
#> (as 'lib' is unspecified)
```

## Example

This is a basic example which shows you how to obtain between subject
confidence intervals. All of the following examples will use the
post-hoc method (see Bird, 2002, for details on the post-hoc method):

``` r

## load relevant packages
library(emmeans)
#> Welcome to emmeans.
#> Caution: You lose important information if you filter this package's results.
#> See '? untidy'
```

``` r
library(afex)
#> Loading required package: lme4
#> Loading required package: Matrix
#> ************
#> Welcome to afex. For support visit: http://afex.singmann.science/
#> - Functions for ANOVAs: aov_car(), aov_ez(), and aov_4()
#> - Methods for calculating p-values with mixed(): 'S', 'KR', 'LRT', and 'PB'
#> - 'afex_aov' and 'mixed' objects can be passed to emmeans() for follow-up tests
#> - Get and set global package options with: afex_options()
#> - Set sum-to-zero contrasts globally: set_sum_contrasts()
#> - For example analyses see: browseVignettes("afex")
#> ************
#> 
#> Attaching package: 'afex'
#> The following object is masked from 'package:lme4':
#> 
#>     lmer
```

``` r
library(readxl)
library(PsyR)

# load some data
dat <- spacing
# load(spacing)?

# set emmeans option to multivariate
afex_options(emmeans_model = "multivariate")

# perform the statistical model
mod <- aov_ez("subj", "yield", dat, within = "spacing", between = "group")
#> Converting to factor: group
#> Contrasts set to contr.sum for the following variables: group
```

``` r

# define some between group contrasts
emm_btwn <- emmeans(mod, "group")
con_b <- list(
  "12vs34" = c(0.5, 0.5, -0.5, -0.5), # groups 1 & 2 vs groups 3 & 4
  "1vs2" = c(1, -1, 0, 0), # and so on...
  "3vs4" = c(0, 0, 1, -1)
)
btwn_con <- contrast(emm_btwn, con_b)

# feed the contrast table into psyci() with the chosen method, family, and factor
# names
# extra documentation needed re: mapping of method to specific. Potentially
# add Sidak
psyci(model=mod, contrast_table = btwn_con, method="ph", family="b", 
      between_factors = list("group"))
#>  contrast estimate   SE df t.ratio p.value  lower upper
#>  12vs34       2.96 2.12 12   1.398  0.1873  -3.89  9.80
#>  1vs2       -14.08 2.99 12  -4.707  0.0005 -23.77 -4.40
#>  3vs4        -2.50 2.99 12  -0.836  0.4197 -12.18  7.18
#> 
#> Results are averaged over the levels of: spacing
```

``` r
# add attributes to the table about the method used to compute CIs
# think about adding the critical constant as an attribute
# specify df also as attributes
```

Here is an example, that carries on from the last, for obtaining within
subject confidence intervals:

``` r

# define some within group contrasts.
# and get emm table of contrast effects
emm_win <- emmeans(mod, "spacing")
con_w <- list(
  "20vs40" = c(1, -1, 0),
  "20vs60" = c(1, 0, -1),
  "Quad" = c(0.5, -1, 0.5)
)
con_win <- contrast(emm_win, con_w)

# generate 95% CIs for the within subjects contrasts
psyci(model=mod, contrast_table = con_win, method="ph", family="w", 
      within_factors = list("spacing"))
#>  contrast estimate    SE df t.ratio p.value lower   upper
#>  20vs40      -1.75 0.832 12  -2.103  0.0573 -4.20  0.7033
#>  20vs60      -3.00 1.010 12  -2.958  0.0120 -5.99 -0.0104
#>  Quad        -0.25 0.459 12  -0.545  0.5956 -1.60  1.1017
#> 
#> Results are averaged over the levels of: group
```

And finally, here is an example of how to generate CIs for between x
within contrasts, using the same post-hoc method:

``` r

# get emms for each cell from the between x within design
emm_int <- emmeans(mod, c("group", "spacing"))

# the handy thing about emmeans is that you can use the already defined between 
# and within contrasts to generate your interaction contrasts. The extra 
# delightful thing is that, when used this way, emmeans will scale the contrasts
# appropriately so that you can interpret the estimated effect as the size of 
# the effect (aka it is scaled appropriately).

con_int <- contrast(emm_int, interaction=list(con_b, con_w))
# check mean difference option in Psy (go back to win app and doc)
# potentially add a scheffe function 
# add p values from gcr etc procedures to the table
# generate 95% CIs for the between x within subjects contrasts
psyci(model=mod, contrast_table = con_int, method="ph", family="bw", 
      within_factors = list("spacing"), between_factors=list("group"))
#>  group_custom spacing_custom estimate    SE df t.ratio p.value  lower  upper
#>  12vs34       20vs40           -8.000 1.660 12  -4.806  0.0004 -14.98  -1.02
#>  1vs2         20vs40            5.000 2.350 12   2.124  0.0551  -4.87  14.87
#>  3vs4         20vs40           -3.000 2.350 12  -1.274  0.2267 -12.87   6.87
#>  12vs34       20vs60          -19.250 2.030 12  -9.490  <.0001 -27.75 -10.75
#>  1vs2         20vs60            8.750 2.870 12   3.050  0.0101  -3.28  20.78
#>  3vs4         20vs60           -6.750 2.870 12  -2.353  0.0365 -18.78   5.28
#>  12vs34       Quad              1.625 0.917 12   1.772  0.1018  -2.22   5.47
#>  1vs2         Quad              0.625 1.300 12   0.482  0.6386  -4.81   6.06
#>  3vs4         Quad              0.375 1.300 12   0.289  0.7774  -5.06   5.81
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
