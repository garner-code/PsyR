library(devtools)
load_all()
#library(PsyR)
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
emms <- emmeans(mod, c("Group", "Time"))

between_contrasts <- list(
    "Ts_v_Ctrl" = c(-2, 1, 1),
    "T1_v_T2" = c(0, 1, -1)
)
between_contrasts <- rescale_contrasts(between_contrasts)
between_selectors <- list(
    "T1*" = c(0, 1, 0),
    "T2*" = c(0, 0, 1),
    "Ctrl*" = c(1, 0, 0)
)
between_contrasts <- c(between_contrasts, between_selectors)

within_contrasts <- list(
    "Post_vs_Pre" = c(-1, 1, 0),
    "FU_vs_Post" = c(0, -1, 1)
)
within_contrasts <- rescale_contrasts(within_contrasts)
within_selectors <- list(
    "Post*" = c(0, 1, 0),
    "FU*" = c(0, 0, 1)
)
within_contrasts <- c(within_contrasts, within_selectors)

interaction_fx <- contrast(
    emms,
    interaction = list(between_contrasts, within_contrasts),
    adjust = "bonferroni"
)
between_fx <- contrast(
    emmeans(mod, "Group"),
    between_contrasts[1:2],
    adjust = "bonferroni"
)
within_fx <- contrast(
    emmeans(mod, "Time"),
    within_contrasts[1:2],
    adjust = "bonferroni"
)

all_fx <- list(interaction_fx, between_fx, within_fx)
this_alpha <- 3 * 0.05
contrasts_w_cis <- psyci(
    mod,
    contrast_tables = all_fx,
    method = "bf",
    alpha = this_alpha,
    family_group = rep("all_factorial", length(all_fx))
)
