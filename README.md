
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
pak::pak("garner-code/PsyR")
```

## Cheat Sheet

[Download the PsyR cheat sheet](inst/doc/psyr-cheatsheet.pdf).

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

| contrast | is_polynomial | estimate | SE | df_num | df_error | F_value | p_value_adjusted | cc | lower | upper |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12vs34 | FALSE | 2.958 | 2.116 | 1 | 12 | 1.955 | 0.597 | 3.236 | -3.888 | 9.805 |
| 1vs2 | FALSE | -14.083 | 2.992 | 1 | 12 | 22.155 | 0.005 | 3.236 | -23.765 | -4.401 |
| 3vs4 | FALSE | -2.500 | 2.992 | 1 | 12 | 0.698 | 0.872 | 3.236 | -12.182 | 7.182 |

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

| contrast | is_polynomial | estimate | SE | df_num | df_error | F_value | p_value_adjusted | cc | lower | upper |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12vs34 | FALSE | 0.640 | 0.457 | 1 | 12 | 1.955 | 0.597 | 3.236 | -0.841 | 2.120 |
| 1vs2 | FALSE | -3.045 | 0.647 | 1 | 12 | 22.155 | 0.005 | 3.236 | -5.138 | -0.952 |
| 3vs4 | FALSE | -0.540 | 0.647 | 1 | 12 | 0.698 | 0.872 | 3.236 | -2.634 | 1.553 |

It is possible to print Psy’s legacy layout for the supplied contrasts.
This will contain values for each contrast’s sum of squares, degrees of
freedom, mean squares error, and F statistic.

To print Psy’s legacy layout, set `legacy_summary` as TRUE within the
`psyci()` call.

``` r
result <- psyci(
  mod,
  list(btwn_con),
  method = "ph",
  legacy_summary = TRUE
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
