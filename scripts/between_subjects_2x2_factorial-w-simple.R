library(PsyR)
library(tidyverse)
library(emmeans)
library(afex)
library(forcats)

data(drive_sleepy)
drive_sleepy <- drive_sleepy %>%
    mutate(
        sleep_deprivation = factor(sleep_deprivation),
        sleep_deprivation = fct_recode(
            sleep_deprivation,
            "None" = "0",
            "High" = "12"
        ),
        sleep_deprivation = fct_relevel(
            sleep_deprivation,
            "High",
            "None"
        ),
        nvh = factor(nvh),
        nvh = fct_relevel(nvh, "high", "low")
    )

mod <- aov_ez(
    "subID",
    "performance",
    drive_sleepy,
    between = c("sleep_deprivation", "nvh")
)
emms <- emmeans(mod, c("sleep_deprivation", "nvh"))

btwn_cont <- list(
    "Sleep_Dep" = c(0.5, -0.5, 0.5, -0.5),
    "NVH" = c(0.5, 0.5, -0.5, -0.5)
)
btwn_cont <- rescale_contrasts(btwn_cont)
mfx <- contrast(emms, btwn_cont)

interaction_contrast <- list(
    "int" = c(0.5, -0.5, -0.5, 0.5)
)
interaction_contrast <- rescale_contrasts(
    interaction_contrast,
    mode = "interaction",
    interaction_order = 1
)
int <- contrast(emms, interaction_contrast)

simp_fx_contrast <- list("NVH" = c(1, -1))
simp_fx_contrast <- rescale_contrasts(simp_fx_contrast)
nvh_simp_fx <- contrast(
    emmeans(mod, ~ nvh | sleep_deprivation),
    method = simp_fx_contrast
)
sdep_simp_fx <- contrast(
    emmeans(mod, ~ sleep_deprivation | nvh),
    method = simp_fx_contrast
)

all_fx <- c(mfx, int, nvh_simp_fx, sdep_simp_fx)
smr_params <- list(p = 2, q = 2, n_sim = 100000, seed = 42)
contrasts_w_cis <- psyci(
    model = mod,
    contrast_tables = all_fx,
    method = "smr",
    smr_params = smr_params,
    alpha = 0.05
)
