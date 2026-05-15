
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PsyR

<!-- badges: start -->
<!-- badges: end -->

PsyR computes confidence intervals for planned and post-hoc contrasts
from ANOVA models, including procedures described by Bird (2002) for
psychology research.

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
library(afex)
library(PsyR)

# load some data
dat <- spacing

# set emmeans option to multivariate
afex_options(emmeans_model = "multivariate")

# perform the statistical model
mod <- aov_ez("subj", "yield", dat, within = "spacing", between = "group")

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
psyci(model=mod, contrast_tables = list(btwn_con), method="ph", family_list=list("b"),
      between_factors = list("group"))
#> $b
#>  contrast estimate   SE df t.ratio p.value   cc  lower upper
#>  12vs34       2.96 2.12 12   1.398  0.1873 3.24  -3.89  9.80
#>  1vs2       -14.08 2.99 12  -4.707  0.0005 3.24 -23.77 -4.40
#>  3vs4        -2.50 2.99 12  -0.836  0.4197 3.24 -12.18  7.18
#> 
#> Results are averaged over the levels of: spacing 
#> PsyR CI method: Scheffe has been applied 
#> Family-wise correction assumes current contrasts are: between subject contrasts 
#> PsyR used an alpha rate of: 0.05 
#> PsyR assumed between subject factor(s) are: group 
#> PsyR used df between of: 3 
#> PsyR used df error of: 12
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
psyci(model=mod, contrast_tables = list(con_win), method="ph", family_list=list("w"), 
      within_factors = list("spacing"))
#> $w
#>  contrast estimate    SE df t.ratio p.value   cc lower   upper
#>  20vs40      -1.75 0.832 12  -2.103  0.0573 2.95 -4.20  0.7033
#>  20vs60      -3.00 1.010 12  -2.958  0.0120 2.95 -5.99 -0.0104
#>  Quad        -0.25 0.459 12  -0.545  0.5956 2.95 -1.60  1.1017
#> 
#> Results are averaged over the levels of: group 
#> PsyR CI method: post-hoc within has been applied 
#> Family-wise correction assumes current contrasts are: within subject contrasts 
#> PsyR used an alpha rate of: 0.05 
#> PsyR assumed within subject factor(s) are: spacing 
#> PsyR used df within of: 2 
#> PsyR used df error of: 12
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
# generate 95% CIs for the between x within subjects contrasts
psyci(model=mod, contrast_tables = list(con_int), method="ph", family_list=list("bw"),
      within_factors = list("spacing"), between_factors=list("group"))
#> $bw
#>  group_custom spacing_custom estimate    SE df t.ratio p.value   cc  lower
#>  12vs34       20vs40           -8.000 1.660 12  -4.806  0.0004 4.19 -14.98
#>  1vs2         20vs40            5.000 2.350 12   2.124  0.0551 4.19  -4.87
#>  3vs4         20vs40           -3.000 2.350 12  -1.274  0.2267 4.19 -12.87
#>  12vs34       20vs60          -19.250 2.030 12  -9.490  <.0001 4.19 -27.75
#>  1vs2         20vs60            8.750 2.870 12   3.050  0.0101 4.19  -3.28
#>  3vs4         20vs60           -6.750 2.870 12  -2.353  0.0365 4.19 -18.78
#>  12vs34       Quad              1.625 0.917 12   1.772  0.1018 4.19  -2.22
#>  1vs2         Quad              0.625 1.300 12   0.482  0.6386 4.19  -4.81
#>  3vs4         Quad              0.375 1.300 12   0.289  0.7774 4.19  -5.06
#>   upper
#>   -1.02
#>   14.87
#>    6.87
#>  -10.75
#>   20.78
#>    5.28
#>    5.47
#>    6.06
#>    5.81
#> 
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied 
#> Family-wise correction assumes current contrasts are: between x within subject contrasts 
#> PsyR used an alpha rate of: 0.05 
#> PsyR assumed between subject factor(s) are: group 
#> PsyR used df between of: 3 
#> PsyR assumed within subject factor(s) are: spacing 
#> PsyR used df within of: 2 
#> PsyR used df error of: 12
```

Or even better, you can do all three at once:

``` r

# enter contrast tables in one long list

psyci(model=mod, contrast_tables = list(btwn_con, con_win, con_int), 
      method="ph", family_list=list("b", "w", "bw"), 
      within_factors = list("spacing"), 
      between_factors=list("group"))
#> $b
#>  contrast estimate   SE df t.ratio p.value   cc  lower upper
#>  12vs34       2.96 2.12 12   1.398  0.1873 3.24  -3.89  9.80
#>  1vs2       -14.08 2.99 12  -4.707  0.0005 3.24 -23.77 -4.40
#>  3vs4        -2.50 2.99 12  -0.836  0.4197 3.24 -12.18  7.18
#> 
#> Results are averaged over the levels of: spacing 
#> PsyR CI method: Scheffe has been applied 
#> Family-wise correction assumes current contrasts are: between subject contrasts 
#> PsyR used an alpha rate of: 0.05 
#> PsyR assumed between subject factor(s) are: group 
#> PsyR used df between of: 3 
#> PsyR used df error of: 12 
#> 
#> $w
#>  contrast estimate    SE df t.ratio p.value   cc lower   upper
#>  20vs40      -1.75 0.832 12  -2.103  0.0573 2.95 -4.20  0.7033
#>  20vs60      -3.00 1.010 12  -2.958  0.0120 2.95 -5.99 -0.0104
#>  Quad        -0.25 0.459 12  -0.545  0.5956 2.95 -1.60  1.1017
#> 
#> Results are averaged over the levels of: group 
#> PsyR CI method: post-hoc within has been applied 
#> Family-wise correction assumes current contrasts are: within subject contrasts 
#> PsyR used an alpha rate of: 0.05 
#> PsyR assumed within subject factor(s) are: spacing 
#> PsyR used df within of: 2 
#> PsyR used df error of: 12 
#> 
#> $bw
#>  group_custom spacing_custom estimate    SE df t.ratio p.value   cc  lower
#>  12vs34       20vs40           -8.000 1.660 12  -4.806  0.0004 4.19 -14.98
#>  1vs2         20vs40            5.000 2.350 12   2.124  0.0551 4.19  -4.87
#>  3vs4         20vs40           -3.000 2.350 12  -1.274  0.2267 4.19 -12.87
#>  12vs34       20vs60          -19.250 2.030 12  -9.490  <.0001 4.19 -27.75
#>  1vs2         20vs60            8.750 2.870 12   3.050  0.0101 4.19  -3.28
#>  3vs4         20vs60           -6.750 2.870 12  -2.353  0.0365 4.19 -18.78
#>  12vs34       Quad              1.625 0.917 12   1.772  0.1018 4.19  -2.22
#>  1vs2         Quad              0.625 1.300 12   0.482  0.6386 4.19  -4.81
#>  3vs4         Quad              0.375 1.300 12   0.289  0.7774 4.19  -5.06
#>   upper
#>   -1.02
#>   14.87
#>    6.87
#>  -10.75
#>   20.78
#>    5.28
#>    5.47
#>    6.06
#>    5.81
#> 
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied 
#> Family-wise correction assumes current contrasts are: between x within subject contrasts 
#> PsyR used an alpha rate of: 0.05 
#> PsyR assumed between subject factor(s) are: group 
#> PsyR used df between of: 3 
#> PsyR assumed within subject factor(s) are: spacing 
#> PsyR used df within of: 2 
#> PsyR used df error of: 12
```
