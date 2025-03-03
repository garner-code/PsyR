
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
psyci(model=mod, contrast_table = btwn_con, method="ph", family="b", 
      between_factors = list("group"))
#>  contrast estimate   SE df t.ratio p.value  lower upper
#>  12vs34       2.96 2.12 12   1.398  0.1873  -3.89  9.80
#>  1vs2       -14.08 2.99 12  -4.707  0.0005 -23.77 -4.40
#>  3vs4        -2.50 2.99 12  -0.836  0.4197 -12.18  7.18
#> 
#> Results are averaged over the levels of: spacing
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
#>  20vs60      -3.00 1.014 12  -2.958  0.0120 -5.99 -0.0104
#>  Quad        -0.25 0.459 12  -0.545  0.5956 -1.60  1.1017
#> 
#> Results are averaged over the levels of: group
```

And finally, here is an example of how to generate CIs for between x
within contrasts, using the same post-hoc method:

``` r

# get emms for each cell from the between x within design
emm_int <- emmeans(mod, c("group", "spacing"))

# sometimes defining the between x within contrast is a little more involved,
# but here is one way. You need to wind up with a contrast vector that has
# as many elements as there are rows in your table of emms for each cell in 
# the between x within design
n_within = 3
n_group = 4
win_full <- lapply(con_w, function(x) rep(x, each = n_group))
btwn_full <- lapply(con_b, function(x) rep(x, times = n_within))
con_i <- lapply(btwn_full, function(x) lapply(win_full, function(y) x * y))
con_int <- contrast(emm_int, con_i)

# generate 95% CIs for the between x within subjects contrasts
psyci(model=mod, contrast_table = con_win, method="ph", family="bw", 
      within_factors = list("spacing"), between_factors=list("group"))
#>  contrast estimate    SE df t.ratio p.value lower upper
#>  20vs40      -1.75 0.832 12  -2.103  0.0573 -5.24  1.74
#>  20vs60      -3.00 1.014 12  -2.958  0.0120 -7.25  1.25
#>  Quad        -0.25 0.459 12  -0.545  0.5956 -2.17  1.67
#> 
#> Results are averaged over the levels of: group
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
