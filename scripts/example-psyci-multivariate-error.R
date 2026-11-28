# Demonstrate PsyR's error for univariate within-subject contrast tables.
# Run from the PsyR repository root. The final call deliberately stops.
devtools::load_all()
library(afex)
library(emmeans)

data("spacing", package = "PsyR")
spacing$group <- factor(spacing$group, levels = 1:4)
spacing$spacing <- factor(
    spacing$spacing,
    levels = c("TWENTY", "FORTY", "SIXTY")
)

# include_aov = TRUE retains the univariate fit needed for this demonstration.
# Without it, afex may substitute the multivariate fit automatically.
mod <- aov_ez(
    id = "subj",
    dv = "yield",
    data = spacing,
    between = "group",
    within = "spacing",
    include_aov = TRUE
)

within_contrasts <- list("Twenty - Forty" = c(1, -1, 0))

# Deliberately select the incompatible basis without changing global options.
win_fx <- contrast(
    emmeans(mod, "spacing", model = "univariate"),
    method = within_contrasts,
    adjust = "none"
)

# Expected error identifies 'within' and explains how to rebuild the table.
psyci(
    model = mod,
    contrast_tables = list(within = win_fx),
    method = "bf"
)

# To fix the error, run these lines separately after inspecting the message:
# win_fx <- contrast(
#     emmeans(mod, "spacing", model = "multivariate"),
#     method = within_contrasts,
#     adjust = "none"
# )
# psyci(mod, list(within = win_fx), method = "bf")
