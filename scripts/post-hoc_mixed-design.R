library(afex)
library(emmeans)
library(PsyR)
library(tidyverse)

data(social_anxiety)
social_anxiety <- social_anxiety %>%
    mutate(
        Session = fct_relevel(Session, "Pre", "Post", "FU1", "FU2"),
        Group = factor(Group),
        Group = fct_recode(
            Group,
            "New" = "1",
            "Standard" = "2",
            "MinContact" = "3"
        )
    )

afex_options(emmeans_model = "multivariate")
mod <- aov_ez(
    "Subject",
    "Score",
    social_anxiety,
    between = "Group",
    within = "Session"
)

btwn_contrasts <- list(
    "Ts - C" = c(0.5, 0.5, -1),
    "NT - ST" = c(1, -1, 0)
)
btwn_contrasts <- rescale_contrasts(btwn_contrasts)
win_contrasts <- list(
    "Pre - Rest" = c(1, -1 / 3, -1 / 3, -1 / 3),
    "Post - FUs" = c(0, 1, -0.5, -0.5),
    "FU1 - FU2" = c(0, 0, 1, -1)
)
win_contrasts <- rescale_contrasts(win_contrasts)

btwn_fx <- contrast(emmeans(mod, "Group"), btwn_contrasts)
win_fx <- contrast(emmeans(mod, "Session"), win_contrasts)
int_fx <- contrast(
    emmeans(mod, c("Group", "Session")),
    interaction = list(btwn_contrasts, win_contrasts)
)

all_fx <- c(win_fx, btwn_fx, int_fx)
alphas <- 0.05
planned_results <- psyci(
    model = mod,
    contrast_tables = all_fx,
    method = "bf",
    alpha = alphas
)
post_hoc_results <- psyci(
    model = mod,
    contrast_tables = all_fx,
    method = "ph",
    alpha = alphas
)
