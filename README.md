
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PsyR

<!-- badges: start -->

<!-- badges: end -->

PsyR utilises procedures described by Bird (2002) to compute
simultaneous confidence intervals for planned and post-hoc contrasts
from ANOVA models. Simultaneous confidence intervals give a researcher
greater insight into the magnitude of an effect captured by a contrast.
This can allow for more nuanced interpretation of results with respect
to directionality and effect size.

## Installation

Install the latest version of PsyR by running the below code.

``` r
install.packages("PsyR")
```

Alternatively, you can install the development version of PsyR from
[GitHub](eventual_link_UNSW_git). Only do so if you want to access early
features. Warning: this version may contain bugs!

``` r
# install.packages("devtools")
devtools::install_github("garner-code/PsyR")
```

## Cheat Sheet

Click [here](link%20to%20be%20inserted) to view and download a cheat
sheet for the package.

## Classic Workflow

An average workflow will go as follows:

- Load data, tidy into longform.

- Use `afex::aov_ez()` to generate an ANOVA model.

- Use `emmeans::emmeans()` to obtain estimated marginal means.

- Use `base::list()` to define contrasts.

- Use `emmeans::contrast()` to generate statistical output pertaining to
  each contrast.

- Use `psyci()` to generate simultaneous confidence intervals.

## Example

Below are some basic examples which shows you how to obtain simultaneous
confidence intervals. All of the following examples will use the
post-hoc method. See (Bird 2002) for details. Also, check out the
“Articles” tab for detailed examples of analysis pertaining to various
designs.

``` r
## load relevant packages
library(emmeans)
library(afex)
library(PsyR)
library(tidyverse)

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
  "1vs2" = c(1, -1, 0, 0),    # group 1 vs group 2
  "3vs4" = c(0, 0, 1, -1)     # group 3 vs group 4
)

#mean difference scale contrasts, then produce and store statistical output
con_b <- rescale_contrasts(raw_con_b, mode = "mean_difference")
btwn_con <- contrast(emm_btwn, con_b)

# feed the contrast table into psyci(); its family and participating factors
# are inferred from the model and the original emmGrid
result <- psyci(
  model = mod, 
  contrast_tables = list(btwn_con), 
  method = "ph"
)
```

| contrast | estimate | SE | df_num | df_error | F_value | p_value_adjusted | cc | lower | upper |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12vs34 | 2.958 | 2.116 | 1 | 12 | 1.955 | 0.597 | 3.236 | -3.888 | 9.805 |
| 1vs2 | -14.083 | 2.992 | 1 | 12 | 22.155 | 0.005 | 3.236 | -23.765 | -4.401 |
| 3vs4 | -2.500 | 2.992 | 1 | 12 | 0.698 | 0.872 | 3.236 | -12.182 | 7.182 |

`psyci()` prints each result twice. The first table reports the contrast
effects and confidence intervals in the dependent variable’s original
units. The second reports the same values in sample standard-deviation
units, using the pooled within-cell scaling method from the original Psy
program (Bird 2002). The raw tables remain available through ordinary
list indexing. For example:

``` r
result <- psyci(mod, list(btwn_con), method = "ph")
raw_table <- result[[1]]
scaled_table <- attr(result, "scaled_tables")[[1]]
```

| contrast | estimate | SE | df_num | df_error | F_value | p_value_adjusted | cc | lower | upper |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12vs34 | 0.640 | 0.457 | 1 | 12 | 1.955 | 0.597 | 3.236 | -0.841 | 2.120 |
| 1vs2 | -3.045 | 0.647 | 1 | 12 | 22.155 | 0.005 | 3.236 | -5.138 | -0.952 |
| 3vs4 | -0.540 | 0.647 | 1 | 12 | 0.698 | 0.872 | 3.236 | -2.634 | 1.553 |

It is possible to print Psy’s legacy SS/df/MS/F layout for the supplied
contrasts. The legacy table is printed once, before the raw and
standardized confidence interval tables. It includes only supplied
contrasts and their applicable error rows; it does not create absent
between-by-within cross-products.

To print Psy’s legacy layout, set `legay_anova` as TRUE within the
`psyci()` call.

``` r
result <- psyci(
  mod,
  list(btwn_con),
  method = "ph",
  legacy_anova = TRUE
)
legacy_table <- attr(result, "legacy_anova_table")
```

| section | label  | source | row_type |       SS |  df |       MS | F_value |
|:--------|:-------|:-------|:---------|---------:|----:|---------:|--------:|
| Between | 12vs34 | B1     | contrast |  105.021 |   1 |  105.021 |   1.955 |
| Between | 1vs2   | B2     | contrast | 1190.042 |   1 | 1190.042 |  22.155 |
| Between | 3vs4   | B3     | contrast |   37.500 |   1 |   37.500 |   0.698 |
| Between |        | Error  | error    |  644.583 |  12 |   53.715 |      NA |

Here is an example, that carries on from the last, for obtaining
within-subjects confidence intervals:

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
win_result <- psyci(
  model = mod, 
  contrast_tables = list(con_win),
  method = "ph"
)
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
psyci(
  model = mod, 
  contrast_tables = list(con_int), 
  method = "ph"
  )
#> Raw contrast effects and confidence intervals (dependent-variable units)
#> 
#> $bw
#>   group_custom spacing_custom estimate       SE df_num df_error     F_value
#> 1       12vs34         20vs40   -8.000 1.664582      1       12 23.09774436
#> 2         1vs2         20vs40    5.000 2.354074      1       12  4.51127820
#> 3         3vs4         20vs40   -3.000 2.354074      1       12  1.62406015
#> 4       12vs34         20vs60  -19.250 2.028444      1       12 90.06075949
#> 5         1vs2         20vs60    8.750 2.868652      1       12  9.30379747
#> 6         3vs4         20vs60   -6.750 2.868652      1       12  5.53670886
#> 7       12vs34           Quad    1.625 0.917140      1       12  3.13931889
#> 8         1vs2           Quad    0.625 1.297032      1       12  0.23219814
#> 9         3vs4           Quad    0.375 1.297032      1       12  0.08359133
#>         theta      p_value p_value_adjusted p_adjustment       cc      lower
#> 1 0.658097686 4.291254e-04     2.157666e-02      Roy GCR 4.191986 -14.977905
#> 2 0.273224044 5.514322e-02     5.884826e-01      Roy GCR 4.191986  -4.868248
#> 3 0.119205298 2.266526e-01     9.297047e-01      Roy GCR 4.191986 -12.868248
#> 4 0.882422980 6.285730e-07     8.158899e-05      Roy GCR 4.191986 -27.753208
#> 5 0.436720143 1.008073e-02     2.240525e-01      Roy GCR 4.191986  -3.275351
#> 6 0.315721091 3.651164e-02     4.807223e-01      Roy GCR 4.191986 -18.775351
#> 7 0.207361963 1.017931e-01     7.546456e-01      Roy GCR 4.191986  -2.219638
#> 8 0.018982536 6.385676e-01     9.995511e-01      Roy GCR 4.191986  -4.812140
#> 9 0.006917756 7.774224e-01     9.999771e-01      Roy GCR 4.191986  -5.062140
#>        upper
#> 1  -1.022095
#> 2  14.868248
#> 3   6.868248
#> 4 -10.746792
#> 5  20.775351
#> 6   5.275351
#> 7   5.469638
#> 8   6.062140
#> 9   5.812140
#> 
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied
#> PsyR reports `F_value` as the ordinary F(1, df_error) contrast statistic and `p_value` as its unadjusted upper-tail probability.
#> PsyR simultaneous-inference procedure: Roy GCR. `p_value_adjusted` and the confidence interval use this procedure.
#> Psy probability-calculator parameters for Roy GCR: s = 2, m = 0, n = 4.5. They are stored in `psyr_probability_parameters`; the `theta` column contains the value for Psy's GCR field.
#> Family-wise correction assumes current contrasts are: between x within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used post-hoc between dimension of: 3
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used post-hoc within dimension of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality with other tables was not applicable because only one table was supplied.
#> 
#> Scaled contrast effects and confidence intervals (sample SD units)
#> Scaling divisor (pooled sample SD): 4.625563
#> 
#> $bw
#>   group_custom spacing_custom    estimate        SE df_num df_error     F_value
#> 1       12vs34         20vs40 -1.72951919 0.3598658      1       12 23.09774436
#> 2         1vs2         20vs40  1.08094949 0.5089271      1       12  4.51127820
#> 3         3vs4         20vs40 -0.64856969 0.5089271      1       12  1.62406015
#> 4       12vs34         20vs60 -4.16165554 0.4385290      1       12 90.06075949
#> 5         1vs2         20vs60  1.89166161 0.6201737      1       12  9.30379747
#> 6         3vs4         20vs60 -1.45928181 0.6201737      1       12  5.53670886
#> 7       12vs34           Quad  0.35130858 0.1982764      1       12  3.13931889
#> 8         1vs2           Quad  0.13511869 0.2804052      1       12  0.23219814
#> 9         3vs4           Quad  0.08107121 0.2804052      1       12  0.08359133
#>         theta      p_value p_value_adjusted p_adjustment       cc      lower
#> 1 0.658097686 4.291254e-04     2.157666e-02      Roy GCR 4.191986 -3.2380717
#> 2 0.273224044 5.514322e-02     5.884826e-01      Roy GCR 4.191986 -1.0524660
#> 3 0.119205298 2.266526e-01     9.297047e-01      Roy GCR 4.191986 -2.7819852
#> 4 0.882422980 6.285730e-07     8.158899e-05      Roy GCR 4.191986 -5.9999631
#> 5 0.436720143 1.008073e-02     2.240525e-01      Roy GCR 4.191986 -0.7080979
#> 6 0.315721091 3.651164e-02     4.807223e-01      Roy GCR 4.191986 -4.0590413
#> 7 0.207361963 1.017931e-01     7.546456e-01      Roy GCR 4.191986 -0.4798634
#> 8 0.018982536 6.385676e-01     9.995511e-01      Roy GCR 4.191986 -1.0403360
#> 9 0.006917756 7.774224e-01     9.999771e-01      Roy GCR 4.191986 -1.0943835
#>        upper
#> 1 -0.2209666
#> 2  3.2143650
#> 3  1.4848458
#> 4 -2.3233480
#> 5  4.4914211
#> 6  1.1404777
#> 7  1.1824806
#> 8  1.3105734
#> 9  1.2565259
#> 
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied
#> PsyR reports `F_value` as the ordinary F(1, df_error) contrast statistic and `p_value` as its unadjusted upper-tail probability.
#> PsyR simultaneous-inference procedure: Roy GCR. `p_value_adjusted` and the confidence interval use this procedure.
#> Psy probability-calculator parameters for Roy GCR: s = 2, m = 0, n = 4.5. They are stored in `psyr_probability_parameters`; the `theta` column contains the value for Psy's GCR field.
#> Family-wise correction assumes current contrasts are: between x within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used post-hoc between dimension of: 3
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used post-hoc within dimension of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality with other tables was not applicable because only one table was supplied.
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
```

The “contrast_tables” argument can also take in all stored contrasts
within a `base::list()`

``` r
# enter contrast tables in one long list
psyci(
  model = mod, 
  contrast_tables = list(btwn_con, con_win, con_int), 
  method = "ph"
  )
#> Raw contrast effects and confidence intervals (dependent-variable units)
#> 
#> $b
#>   contrast   estimate       SE df_num df_error    F_value      p_value
#> 1   12vs34   2.958333 2.115720      1       12  1.9551390 0.1873413489
#> 2     1vs2 -14.083333 2.992081      1       12 22.1546218 0.0005082834
#> 3     3vs4  -2.500000 2.992081      1       12  0.6981254 0.4197397783
#>   p_value_adjusted p_adjustment       cc      lower     upper
#> 1      0.596976752      Scheffe 3.235875  -3.887873  9.804540
#> 2      0.004609298      Scheffe 3.235875 -23.765331 -4.401336
#> 3      0.871845486      Scheffe 3.235875 -12.181998  7.181998
#> 
#> Results are averaged over the levels of: spacing
#> PsyR CI method: Scheffe has been applied
#> PsyR reports `F_value` as the ordinary F(1, df_error) contrast statistic and `p_value` as its unadjusted upper-tail probability.
#> PsyR simultaneous-inference procedure: Scheffe. `p_value_adjusted` and the confidence interval use this procedure.
#> Family-wise correction assumes current contrasts are: between subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used post-hoc between dimension of: 3
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[2]], contrast_tables[[3]].
#> 
#> $w
#>   contrast estimate       SE df_num df_error   F_value    p_value
#> 1   20vs40    -1.75 0.832291      1       12 4.4210526 0.05727679
#> 2   20vs60    -3.00 1.014222      1       12 8.7493671 0.01196653
#> 3     Quad    -0.25 0.458570      1       12 0.2972136 0.59562225
#>   p_value_adjusted        p_adjustment       cc     lower      upper
#> 1       0.17815228 Hotelling T-squared 2.947652 -4.203304  0.7033040
#> 2       0.04920057 Hotelling T-squared 2.947652 -5.989572 -0.0104275
#> 3       0.87409756 Hotelling T-squared 2.947652 -1.601705  1.1017046
#> 
#> Results are averaged over the levels of: group
#> PsyR CI method: post-hoc within (Hotelling's T-squared) has been applied
#> PsyR reports `F_value` as the ordinary F(1, df_error) contrast statistic and `p_value` as its unadjusted upper-tail probability.
#> PsyR simultaneous-inference procedure: Hotelling T-squared. `p_value_adjusted` and the confidence interval use this procedure.
#> Family-wise correction assumes current contrasts are: within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used post-hoc within dimension of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[1]], contrast_tables[[3]].
#> 
#> $bw
#>   group_custom spacing_custom estimate       SE df_num df_error     F_value
#> 1       12vs34         20vs40   -8.000 1.664582      1       12 23.09774436
#> 2         1vs2         20vs40    5.000 2.354074      1       12  4.51127820
#> 3         3vs4         20vs40   -3.000 2.354074      1       12  1.62406015
#> 4       12vs34         20vs60  -19.250 2.028444      1       12 90.06075949
#> 5         1vs2         20vs60    8.750 2.868652      1       12  9.30379747
#> 6         3vs4         20vs60   -6.750 2.868652      1       12  5.53670886
#> 7       12vs34           Quad    1.625 0.917140      1       12  3.13931889
#> 8         1vs2           Quad    0.625 1.297032      1       12  0.23219814
#> 9         3vs4           Quad    0.375 1.297032      1       12  0.08359133
#>         theta      p_value p_value_adjusted p_adjustment       cc      lower
#> 1 0.658097686 4.291254e-04     2.157666e-02      Roy GCR 4.191986 -14.977905
#> 2 0.273224044 5.514322e-02     5.884826e-01      Roy GCR 4.191986  -4.868248
#> 3 0.119205298 2.266526e-01     9.297047e-01      Roy GCR 4.191986 -12.868248
#> 4 0.882422980 6.285730e-07     8.158899e-05      Roy GCR 4.191986 -27.753208
#> 5 0.436720143 1.008073e-02     2.240525e-01      Roy GCR 4.191986  -3.275351
#> 6 0.315721091 3.651164e-02     4.807223e-01      Roy GCR 4.191986 -18.775351
#> 7 0.207361963 1.017931e-01     7.546456e-01      Roy GCR 4.191986  -2.219638
#> 8 0.018982536 6.385676e-01     9.995511e-01      Roy GCR 4.191986  -4.812140
#> 9 0.006917756 7.774224e-01     9.999771e-01      Roy GCR 4.191986  -5.062140
#>        upper
#> 1  -1.022095
#> 2  14.868248
#> 3   6.868248
#> 4 -10.746792
#> 5  20.775351
#> 6   5.275351
#> 7   5.469638
#> 8   6.062140
#> 9   5.812140
#> 
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied
#> PsyR reports `F_value` as the ordinary F(1, df_error) contrast statistic and `p_value` as its unadjusted upper-tail probability.
#> PsyR simultaneous-inference procedure: Roy GCR. `p_value_adjusted` and the confidence interval use this procedure.
#> Psy probability-calculator parameters for Roy GCR: s = 2, m = 0, n = 4.5. They are stored in `psyr_probability_parameters`; the `theta` column contains the value for Psy's GCR field.
#> Family-wise correction assumes current contrasts are: between x within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used post-hoc between dimension of: 3
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used post-hoc within dimension of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[1]], contrast_tables[[2]].
#> 
#> Scaled contrast effects and confidence intervals (sample SD units)
#> Scaling divisor (pooled sample SD): 4.625563
#> 
#> $b
#>   contrast   estimate        SE df_num df_error    F_value      p_value
#> 1   12vs34  0.6395618 0.4573974      1       12  1.9551390 0.1873413489
#> 2     1vs2 -3.0446744 0.6468576      1       12 22.1546218 0.0005082834
#> 3     3vs4 -0.5404747 0.6468576      1       12  0.6981254 0.4197397783
#>   p_value_adjusted p_adjustment       cc      lower      upper
#> 1      0.596976752      Scheffe 3.235875 -0.8405188  2.1196424
#> 2      0.004609298      Scheffe 3.235875 -5.1378245 -0.9515243
#> 3      0.871845486      Scheffe 3.235875 -2.6336248  1.5526753
#> 
#> Results are averaged over the levels of: spacing
#> PsyR CI method: Scheffe has been applied
#> PsyR reports `F_value` as the ordinary F(1, df_error) contrast statistic and `p_value` as its unadjusted upper-tail probability.
#> PsyR simultaneous-inference procedure: Scheffe. `p_value_adjusted` and the confidence interval use this procedure.
#> Family-wise correction assumes current contrasts are: between subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used post-hoc between dimension of: 3
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[2]], contrast_tables[[3]].
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
#> 
#> $w
#>   contrast    estimate        SE df_num df_error   F_value    p_value
#> 1   20vs40 -0.37833232 0.1799329      1       12 4.4210526 0.05727679
#> 2   20vs60 -0.64856969 0.2192645      1       12 8.7493671 0.01196653
#> 3     Quad -0.05404747 0.0991382      1       12 0.2972136 0.59562225
#>   p_value_adjusted        p_adjustment       cc      lower        upper
#> 1       0.17815228 Hotelling T-squared 2.947652 -0.9087119  0.152047212
#> 2       0.04920057 Hotelling T-squared 2.947652 -1.2948851 -0.002254321
#> 3       0.87409756 Hotelling T-squared 2.947652 -0.3462724  0.238177416
#> 
#> Results are averaged over the levels of: group
#> PsyR CI method: post-hoc within (Hotelling's T-squared) has been applied
#> PsyR reports `F_value` as the ordinary F(1, df_error) contrast statistic and `p_value` as its unadjusted upper-tail probability.
#> PsyR simultaneous-inference procedure: Hotelling T-squared. `p_value_adjusted` and the confidence interval use this procedure.
#> Family-wise correction assumes current contrasts are: within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used post-hoc within dimension of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[1]], contrast_tables[[3]].
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
#> 
#> $bw
#>   group_custom spacing_custom    estimate        SE df_num df_error     F_value
#> 1       12vs34         20vs40 -1.72951919 0.3598658      1       12 23.09774436
#> 2         1vs2         20vs40  1.08094949 0.5089271      1       12  4.51127820
#> 3         3vs4         20vs40 -0.64856969 0.5089271      1       12  1.62406015
#> 4       12vs34         20vs60 -4.16165554 0.4385290      1       12 90.06075949
#> 5         1vs2         20vs60  1.89166161 0.6201737      1       12  9.30379747
#> 6         3vs4         20vs60 -1.45928181 0.6201737      1       12  5.53670886
#> 7       12vs34           Quad  0.35130858 0.1982764      1       12  3.13931889
#> 8         1vs2           Quad  0.13511869 0.2804052      1       12  0.23219814
#> 9         3vs4           Quad  0.08107121 0.2804052      1       12  0.08359133
#>         theta      p_value p_value_adjusted p_adjustment       cc      lower
#> 1 0.658097686 4.291254e-04     2.157666e-02      Roy GCR 4.191986 -3.2380717
#> 2 0.273224044 5.514322e-02     5.884826e-01      Roy GCR 4.191986 -1.0524660
#> 3 0.119205298 2.266526e-01     9.297047e-01      Roy GCR 4.191986 -2.7819852
#> 4 0.882422980 6.285730e-07     8.158899e-05      Roy GCR 4.191986 -5.9999631
#> 5 0.436720143 1.008073e-02     2.240525e-01      Roy GCR 4.191986 -0.7080979
#> 6 0.315721091 3.651164e-02     4.807223e-01      Roy GCR 4.191986 -4.0590413
#> 7 0.207361963 1.017931e-01     7.546456e-01      Roy GCR 4.191986 -0.4798634
#> 8 0.018982536 6.385676e-01     9.995511e-01      Roy GCR 4.191986 -1.0403360
#> 9 0.006917756 7.774224e-01     9.999771e-01      Roy GCR 4.191986 -1.0943835
#>        upper
#> 1 -0.2209666
#> 2  3.2143650
#> 3  1.4848458
#> 4 -2.3233480
#> 5  4.4914211
#> 6  1.1404777
#> 7  1.1824806
#> 8  1.3105734
#> 9  1.2565259
#> 
#> PsyR CI method: post-hoc between x within (Roy's GCR) has been applied
#> PsyR reports `F_value` as the ordinary F(1, df_error) contrast statistic and `p_value` as its unadjusted upper-tail probability.
#> PsyR simultaneous-inference procedure: Roy GCR. `p_value_adjusted` and the confidence interval use this procedure.
#> Psy probability-calculator parameters for Roy GCR: s = 2, m = 0, n = 4.5. They are stored in `psyr_probability_parameters`; the `theta` column contains the value for Psy's GCR field.
#> Family-wise correction assumes current contrasts are: between x within subject contrasts
#> PsyR used an alpha rate of: 0.05
#> PsyR assumed between subject factor(s) are: group
#> PsyR used post-hoc between dimension of: 3
#> PsyR assumed within subject factor(s) are: spacing
#> PsyR used post-hoc within dimension of: 2
#> PsyR used df error of: 12
#> PsyR found no non-contrast rows to drop.
#> Orthogonality: this table is orthogonal to contrast_tables[[1]], contrast_tables[[2]].
#> PsyR scaled estimate, SE, and confidence limits by sample SD: 4.625563
```

## Getting Help

To report a bug, create an issue [here](eventual_link_UNSW_git).

------------------------------------------------------------------------

<img src="./man/figures/logo.png" alt="" width="50%" style="display: block; margin: auto;" />

------------------------------------------------------------------------

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bird2002" class="csl-entry">

Bird, Kevin D. 2002. “Confidence Intervals for Effect Sizes in Analysis
of Variance.” *Educational and Psychological Measurement* 62 (2):
197–226. <https://doi.org/10.1177/0013164402062002001>.

</div>

</div>
