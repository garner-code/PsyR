# Basic Bonferroni example using the spacing data from the mixed-design vignette.
# Run from the PsyR repository root to load the current development version.
devtools::load_all()
library(afex)
library(emmeans)

data("spacing", package = "PsyR")
spacing$group <- factor(spacing$group, levels = 1:4)
spacing$spacing <- factor(
    spacing$spacing,
    levels = c("TWENTY", "FORTY", "SIXTY")
)

mod <- aov_ez(
    id = "subj",
    dv = "yield",
    data = spacing,
    between = "group",
    within = "spacing"
)

# Coefficients follow the equally spaced levels: twenty, forty, sixty.
within_contrasts <- list(
    linear = c(-1, 0, 1),
    quadratic = c(1, -2, 1),
    forty_vs_twenty = c(-1, 1, 0)
)
win_fx <- contrast(
    emmeans(mod, "spacing", model = "multivariate"),
    method = within_contrasts,
    adjust = "none"
)

# All three contrasts count toward the Bonferroni adjustment. Polynomial rows
# retain estimates and tests, but their confidence limits are suppressed.
# The ordinary forty-versus-twenty comparison retains its confidence interval.
contrasts_w_cis <- psyci(
    model = mod,
    contrast_tables = win_fx,
    method = "bf"
)
print(contrasts_w_cis)

# Inspect contrasts_w_cis[[1]]$is_polynomial to see the detected trend rows.
