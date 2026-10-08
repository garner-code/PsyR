library(afex)
library(emmeans)
library(tidyverse)
library(devtools)
load_all()
#library(PsyR)

data(experience)
experience <- experience %>%
    mutate(
        Experience = factor(
            Experience,
            levels = c("10yrs", "3yrs", "<1yr")
        ),
        Therapy = factor(
            Therapy,
            levels = c(
                "Behavioural",
                "Cognitive",
                "Psychotherapy",
                "Group"
            )
        )
    )

mod <- aov_ez(
    id = "SubjID",
    dv = "Outcome",
    data = experience,
    between = c("Experience", "Therapy")
)
emms <- emmeans(mod, c("Experience", "Therapy"))

experience_contrasts <- list(
    "lots_vs_little" = c(0, 1, -1),
    "lots_vs_some" = c(1, -1, 0)
)
experience_contrasts <- rescale_contrasts(experience_contrasts)
therapy_contrasts <- list(
    "evidence_vs_woowoo" = c(0.5, 0.5, -0.5, -0.5),
    "behav_vs_cog" = c(1, -1, 0, 0),
    "woo1_vs_woo2" = c(0, 0, 1, -1)
)
therapy_contrasts <- rescale_contrasts(therapy_contrasts)

exp_fx <- contrast(
    emmeans(mod, "Experience"),
    experience_contrasts,
    adjust = "bonferroni"
)
Ab_fx <- contrast(
    emmeans(mod, ~ Experience | Therapy),
    method = experience_contrasts,
    adjust = "bonferroni"
)
interaction_fx <- contrast(
    emms,
    interaction = list(experience_contrasts, therapy_contrasts),
    adjust = "bonferroni"
)

all_fx <- c(exp_fx, Ab_fx, interaction_fx)
planned_results <- psyci(
    model = mod,
    contrast_tables = all_fx,
    method = "bf",
    alpha = rep(list(2 * 0.05), times = 3),
    independent = FALSE
)

nu1 <- length(levels(experience$Therapy)) *
    (length(levels(experience$Experience)) - 1)
post_hoc_results <- psyci(
    model = mod,
    contrast_tables = all_fx,
    method = "ph",
    alpha = rep(list(1 - (1 - 0.05)^2), times = 3),
    independent = FALSE,
    nu1 = nu1
)

thrpy_fx <- contrast(
    emmeans(mod, "Therapy"),
    therapy_contrasts
)
standard_fx <- list(exp_fx, thrpy_fx, interaction_fx)
standard_results <- psyci(
    model = mod,
    contrast_tables = standard_fx,
    method = "bf",
    alpha = 0.05
)
