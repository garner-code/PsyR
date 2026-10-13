
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

# define and mean-difference scale some between-group contrasts
emm_btwn <- emmeans(mod, "group")
raw_con_b <- list(
  "12vs34" = c(1, 1, -1, -1), # groups 1 & 2 vs groups 3 & 4
  "1vs2" = c(1, -1, 0, 0), # and so on...
  "3vs4" = c(0, 0, 1, -1)
)
con_b <- rescale_contrasts(raw_con_b, mode = "mean_difference")
btwn_con <- contrast(emm_btwn, con_b)

# feed the contrast table into psyci(); its family and participating factors
# are inferred from the model and the original emmGrid
psyci(model = mod, contrast_tables = list(btwn_con), method = "ph")
#> Raw contrast effects and confidence intervals (dependent-variable units)
#>
#> $b
#>  contrast estimate   SE df t.ratio p.value   cc  lower upper psyr_p_value
#>  12vs34       2.96 2.12 12   1.398  0.1873 3.24  -3.89  9.80      0.59698
#>  1vs2       -14.08 2.99 12  -4.707  0.0005 3.24 -23.77 -4.40      0.00461
#>  3vs4        -2.50 2.99 12  -0.836  0.4197 3.24 -12.18  7.18      0.87185
#>
#> Results are averaged over the levels of: spacing
#> PsyR CI method: Scheffe has been applied
#> PsyR p-value method: Scheffe has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: between subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used df between of: 3
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality with other tables was not applicable because only one table was supplied.
#>
#> Scaled contrast effects and confidence intervals (sample SD units)
#> Scaling divisor (pooled sample SD): 4.625563
#>
#> $b
#>  contrast estimate    SE df t.ratio p.value   cc  lower  upper psyr_p_value
#>  12vs34       0.64 0.457 12   1.398  0.1873 3.24 -0.841  2.120      0.59698
#>  1vs2        -3.04 0.647 12  -4.707  0.0005 3.24 -5.138 -0.952      0.00461
#>  3vs4        -0.54 0.647 12  -0.836  0.4197 3.24 -2.634  1.553      0.87185
#>
#> Results are averaged over the levels of: spacing
#> PsyR CI method: Scheffe has been applied
#> PsyR p-value method: Scheffe has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: between subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used df between of: 3
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality with other tables was not applicable because only one table was supplied.
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
```

`psyci()` prints each result twice. The first table reports the contrast
effects and confidence intervals in the dependent variable’s original
units. The second reports the same values in sample standard-deviation
units, using the pooled within-cell scaling method from the original Psy
program. The raw tables remain available through ordinary list indexing.
For example:

``` r
result <- psyci(mod, list(btwn_con), method = "ph")
raw_table <- result[[1]]
scaled_table <- attr(result, "scaled_tables")[[1]]
attr(result, "scale_info")$divisor
```

Here is an example, that carries on from the last, for obtaining within
subject confidence intervals:

``` r

# define and mean-difference scale some within-subject contrasts
emm_win <- emmeans(mod, "spacing")
raw_con_w <- list(
  "20vs40" = c(1, -1, 0),
  "20vs60" = c(1, 0, -1),
  "Quad" = c(1, -2, 1)
)
con_w <- rescale_contrasts(raw_con_w, mode = "mean_difference")
con_win <- contrast(emm_win, con_w)

# generate 95% CIs for the within subjects contrasts
psyci(model=mod, contrast_tables = list(con_win), method="ph")
#> Raw contrast effects and confidence intervals (dependent-variable units)
#>
#> $w
#>  contrast estimate    SE df t.ratio p.value   cc lower   upper psyr_p_value
#>  20vs40      -1.75 0.832 12  -2.103  0.0573 2.95 -4.20  0.7033       0.1782
#>  20vs60      -3.00 1.010 12  -2.958  0.0120 2.95 -5.99 -0.0104       0.0492
#>  Quad        -0.25 0.459 12  -0.545  0.5956 2.95 -1.60  1.1017       0.8741
#>
#> Results are averaged over the levels of: group
#> PsyR CI method: post-hoc within has been applied
#> PsyR p-value method: post-hoc within has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used df within of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality with other tables was not applicable because only one table was supplied.
#>
#> Scaled contrast effects and confidence intervals (sample SD units)
#> Scaling divisor (pooled sample SD): 4.625563
#>
#> $w
#>  contrast estimate     SE df t.ratio p.value   cc  lower    upper psyr_p_value
#>  20vs40     -0.378 0.1800 12  -2.103  0.0573 2.95 -0.909  0.15205       0.1782
#>  20vs60     -0.649 0.2190 12  -2.958  0.0120 2.95 -1.295 -0.00225       0.0492
#>  Quad       -0.054 0.0991 12  -0.545  0.5956 2.95 -0.346  0.23818       0.8741
#>
#> Results are averaged over the levels of: group
#> PsyR CI method: post-hoc within has been applied
#> PsyR p-value method: post-hoc within has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used df within of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality with other tables was not applicable because only one table was supplied.
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
```

For pairwise correlated within-subject means, PsyR also provides the
asymmetric noncentral-*t* standardized-effect interval proposed by
Algina and Keselman (2003). It is opt-in, so existing calls retain the
original Psy/Bird output. Use it with individual or Bonferroni
intervals:

``` r
result_ak <- psyci(
  mod,
  list(con_win),
  method = "bf",
  standardized_ci = "algina_keselman"
)
ak_table <- attr(result_ak, "scaled_tables")[[1]]
```

The new `algina_keselman_ci()` function also exposes the calculation
directly when an estimated contrast, its standard error, error degrees
of freedom, and the appropriate sample-SD scaling divisor are already
available.

And finally, here is an example of how to generate CIs for between x
within contrasts, using the same post-hoc method:

``` r

# get emms for each cell from the between x within design
emm_int <- emmeans(mod, c("group", "spacing"))

# emmeans forms products of the supplied between and within coefficients. It
# does not normalize those coefficients. Because con_b and con_w were each
# scaled above as order-0 mean differences, their products are correctly scaled
# order-1 interactions (differences between mean differences).

con_int <- contrast(emm_int, interaction=list(con_b, con_w))
# generate 95% CIs for the between x within subjects contrasts
psyci(model=mod, contrast_tables = list(con_int), method="ph")
#> Raw contrast effects and confidence intervals (dependent-variable units)
#>
#> $bw
#>  group_custom spacing_custom estimate    SE df t.ratio p.value   cc  lower
#>  12vs34       20vs40           -8.000 1.660 12  -4.806  0.0004 4.19 -14.98
#>  1vs2         20vs40            5.000 2.350 12   2.124  0.0551 4.19  -4.87
#>  3vs4         20vs40           -3.000 2.350 12  -1.274  0.2267 4.19 -12.87
#>  12vs34       20vs60          -19.250 2.030 12  -9.490 <0.0001 4.19 -27.75
#>  1vs2         20vs60            8.750 2.870 12   3.050  0.0101 4.19  -3.28
#>  3vs4         20vs60           -6.750 2.870 12  -2.353  0.0365 4.19 -18.78
#>  12vs34       Quad              1.625 0.917 12   1.772  0.1018 4.19  -2.22
#>  1vs2         Quad              0.625 1.300 12   0.482  0.6386 4.19  -4.81
#>  3vs4         Quad              0.375 1.300 12   0.289  0.7774 4.19  -5.06
#>   upper psyr_p_value
#>   -1.02     2.16e-02
#>   14.87     5.88e-01
#>    6.87     9.30e-01
#>  -10.75     8.16e-05
#>   20.78     2.24e-01
#>    5.28     4.81e-01
#>    5.47     7.55e-01
#>    6.06     1.00e+00
#>    5.81     1.00e+00
#>
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied
#> PsyR p-value method: post-hoc between x within (Roy's GCR) has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: between x within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used df between of: 3
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used df within of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality with other tables was not applicable because only one table was supplied.
#>
#> Scaled contrast effects and confidence intervals (sample SD units)
#> Scaling divisor (pooled sample SD): 4.625563
#>
#> $bw
#>  group_custom spacing_custom estimate    SE df t.ratio p.value   cc  lower
#>  12vs34       20vs40          -1.7295 0.360 12  -4.806  0.0004 4.19 -3.238
#>  1vs2         20vs40           1.0809 0.509 12   2.124  0.0551 4.19 -1.052
#>  3vs4         20vs40          -0.6486 0.509 12  -1.274  0.2267 4.19 -2.782
#>  12vs34       20vs60          -4.1617 0.439 12  -9.490 <0.0001 4.19 -6.000
#>  1vs2         20vs60           1.8917 0.620 12   3.050  0.0101 4.19 -0.708
#>  3vs4         20vs60          -1.4593 0.620 12  -2.353  0.0365 4.19 -4.059
#>  12vs34       Quad             0.3513 0.198 12   1.772  0.1018 4.19 -0.480
#>  1vs2         Quad             0.1351 0.280 12   0.482  0.6386 4.19 -1.040
#>  3vs4         Quad             0.0811 0.280 12   0.289  0.7774 4.19 -1.094
#>   upper psyr_p_value
#>  -0.221     2.16e-02
#>   3.214     5.88e-01
#>   1.485     9.30e-01
#>  -2.323     8.16e-05
#>   4.491     2.24e-01
#>   1.140     4.81e-01
#>   1.182     7.55e-01
#>   1.311     1.00e+00
#>   1.257     1.00e+00
#>
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied
#> PsyR p-value method: post-hoc between x within (Roy's GCR) has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: between x within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used df between of: 3
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used df within of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality with other tables was not applicable because only one table was supplied.
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
```

Or even better, you can do all three at once:

``` r

# enter contrast tables in one long list

psyci(model=mod, contrast_tables = list(btwn_con, con_win, con_int),
      method="ph")
#> Raw contrast effects and confidence intervals (dependent-variable units)
#>
#> $b
#>  contrast estimate   SE df t.ratio p.value   cc  lower upper psyr_p_value
#>  12vs34       2.96 2.12 12   1.398  0.1873 3.24  -3.89  9.80      0.59698
#>  1vs2       -14.08 2.99 12  -4.707  0.0005 3.24 -23.77 -4.40      0.00461
#>  3vs4        -2.50 2.99 12  -0.836  0.4197 3.24 -12.18  7.18      0.87185
#>
#> Results are averaged over the levels of: spacing
#> PsyR CI method: Scheffe has been applied
#> PsyR p-value method: Scheffe has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: between subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used df between of: 3
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[2]], contrast_tables[[3]].
#>
#> $w
#>  contrast estimate    SE df t.ratio p.value   cc lower   upper psyr_p_value
#>  20vs40      -1.75 0.832 12  -2.103  0.0573 2.95 -4.20  0.7033       0.1782
#>  20vs60      -3.00 1.010 12  -2.958  0.0120 2.95 -5.99 -0.0104       0.0492
#>  Quad        -0.25 0.459 12  -0.545  0.5956 2.95 -1.60  1.1017       0.8741
#>
#> Results are averaged over the levels of: group
#> PsyR CI method: post-hoc within has been applied
#> PsyR p-value method: post-hoc within has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used df within of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[1]], contrast_tables[[3]].
#>
#> $bw
#>  group_custom spacing_custom estimate    SE df t.ratio p.value   cc  lower
#>  12vs34       20vs40           -8.000 1.660 12  -4.806  0.0004 4.19 -14.98
#>  1vs2         20vs40            5.000 2.350 12   2.124  0.0551 4.19  -4.87
#>  3vs4         20vs40           -3.000 2.350 12  -1.274  0.2267 4.19 -12.87
#>  12vs34       20vs60          -19.250 2.030 12  -9.490 <0.0001 4.19 -27.75
#>  1vs2         20vs60            8.750 2.870 12   3.050  0.0101 4.19  -3.28
#>  3vs4         20vs60           -6.750 2.870 12  -2.353  0.0365 4.19 -18.78
#>  12vs34       Quad              1.625 0.917 12   1.772  0.1018 4.19  -2.22
#>  1vs2         Quad              0.625 1.300 12   0.482  0.6386 4.19  -4.81
#>  3vs4         Quad              0.375 1.300 12   0.289  0.7774 4.19  -5.06
#>   upper psyr_p_value
#>   -1.02     2.16e-02
#>   14.87     5.88e-01
#>    6.87     9.30e-01
#>  -10.75     8.16e-05
#>   20.78     2.24e-01
#>    5.28     4.81e-01
#>    5.47     7.55e-01
#>    6.06     1.00e+00
#>    5.81     1.00e+00
#>
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied
#> PsyR p-value method: post-hoc between x within (Roy's GCR) has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: between x within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used df between of: 3
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used df within of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[1]], contrast_tables[[2]].
#>
#> Scaled contrast effects and confidence intervals (sample SD units)
#> Scaling divisor (pooled sample SD): 4.625563
#>
#> $b
#>  contrast estimate    SE df t.ratio p.value   cc  lower  upper psyr_p_value
#>  12vs34       0.64 0.457 12   1.398  0.1873 3.24 -0.841  2.120      0.59698
#>  1vs2        -3.04 0.647 12  -4.707  0.0005 3.24 -5.138 -0.952      0.00461
#>  3vs4        -0.54 0.647 12  -0.836  0.4197 3.24 -2.634  1.553      0.87185
#>
#> Results are averaged over the levels of: spacing
#> PsyR CI method: Scheffe has been applied
#> PsyR p-value method: Scheffe has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: between subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used df between of: 3
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[2]], contrast_tables[[3]].
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
#>
#> $w
#>  contrast estimate     SE df t.ratio p.value   cc  lower    upper psyr_p_value
#>  20vs40     -0.378 0.1800 12  -2.103  0.0573 2.95 -0.909  0.15205       0.1782
#>  20vs60     -0.649 0.2190 12  -2.958  0.0120 2.95 -1.295 -0.00225       0.0492
#>  Quad       -0.054 0.0991 12  -0.545  0.5956 2.95 -0.346  0.23818       0.8741
#>
#> Results are averaged over the levels of: group
#> PsyR CI method: post-hoc within has been applied
#> PsyR p-value method: post-hoc within has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used df within of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[1]], contrast_tables[[3]].
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
#>
#> $bw
#>  group_custom spacing_custom estimate    SE df t.ratio p.value   cc  lower
#>  12vs34       20vs40          -1.7295 0.360 12  -4.806  0.0004 4.19 -3.238
#>  1vs2         20vs40           1.0809 0.509 12   2.124  0.0551 4.19 -1.052
#>  3vs4         20vs40          -0.6486 0.509 12  -1.274  0.2267 4.19 -2.782
#>  12vs34       20vs60          -4.1617 0.439 12  -9.490 <0.0001 4.19 -6.000
#>  1vs2         20vs60           1.8917 0.620 12   3.050  0.0101 4.19 -0.708
#>  3vs4         20vs60          -1.4593 0.620 12  -2.353  0.0365 4.19 -4.059
#>  12vs34       Quad             0.3513 0.198 12   1.772  0.1018 4.19 -0.480
#>  1vs2         Quad             0.1351 0.280 12   0.482  0.6386 4.19 -1.040
#>  3vs4         Quad             0.0811 0.280 12   0.289  0.7774 4.19 -1.094
#>   upper psyr_p_value
#>  -0.221     2.16e-02
#>   3.214     5.88e-01
#>   1.485     9.30e-01
#>  -2.323     8.16e-05
#>   4.491     2.24e-01
#>   1.140     4.81e-01
#>   1.182     7.55e-01
#>   1.311     1.00e+00
#>   1.257     1.00e+00
#>
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied
#> PsyR p-value method: post-hoc between x within (Roy's GCR) has been applied to `psyr_p_value`; `p.value` retains the original emmeans result and may use a table-specific adjustment.
#> Family-wise correction assumes current contrasts are: between x within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used df between of: 3
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used df within of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[1]], contrast_tables[[2]].
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
```
