library(PsyR)
library(tidyverse)
library(emmeans)
library(afex)

data(depression)
depression <- depression %>%
    mutate(
        Group = fct_relevel(Group, "Ctrl", "Tmt1", "Tmt2"),
        Time = fct_relevel(Time, "Pre", "Post", "FU"),
        Subject = as_factor(Subject)
    )

afex_options(emmeans_model = "multivariate")
mod <- aov_ez(
    "Subject",
    "Happiness",
    depression,
    within = "Time",
    between = "Group"
)

between_contrasts <- list(
    "Ts_v_Ctrl" = c(-1, 0.5, 0.5),
    "T1_v_T2" = c(0, 1, -1)
)
between_contrasts <- rescale_contrasts(between_contrasts)
within_contrasts <- list(
    "Post_vs_Pre" = c(-1, 1, 0),
    "FU_vs_Post" = c(0, -1, 1)
)
within_contrasts <- rescale_contrasts(within_contrasts)

btwn_fx <- contrast(
    emmeans(mod, "Group"),
    between_contrasts
)
win_fx <- contrast(
    emmeans(mod, "Time"),
    within_contrasts
)
int_fx <- contrast(
    emmeans(mod, c("Group", "Time")),
    interaction = list(between_contrasts, within_contrasts)
)

contrasts_w_cis <- psyci(
    model = mod,
    contrast_tables = c(btwn_fx, win_fx, int_fx),
    method = "bf"
)
