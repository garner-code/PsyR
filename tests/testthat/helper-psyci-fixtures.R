# Shared fixtures for end-to-end tests of psyci().

make_between_fixture <- function() {
  data("fatigue", package = "PsyR", envir = environment())

  required_columns <- c("id", "Group", "Errors")
  missing_columns <- setdiff(required_columns, names(fatigue))
  if (length(missing_columns) > 0) {
    stop(
      "The fatigue dataset is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  fatigue$Group <- factor(
    fatigue$Group,
    levels = c(1, 2, 3, 4),
    labels = c("Alert", "Mild", "Moderate", "Extreme")
  )

  model <- afex::aov_ez(
    id = "id",
    dv = "Errors",
    data = fatigue,
    between = "Group"
  )

  emmeans <- emmeans::emmeans(model, "Group")
  between_contrasts <- list(
    "Fatigue - Alert" = c(-1, 1 / 3, 1 / 3, 1 / 3),
    "Extreme - Mod, Mild" = c(0, -1 / 2, -1 / 2, 1),
    "Low - High" = c(1 / 2, 1 / 2, -1 / 2, -1 / 2),
    "Alert - Mild" = c(1, -1, 0, 0)
  )

  contrast_table <- emmeans::contrast(
    emmeans,
    method = between_contrasts,
    adjust = "scheffe"
  )

  list(
    data = fatigue,
    model = model,
    emmeans = emmeans,
    tables = list(contrast_table),
    families = list("b"),
    between_factors = list("Group"),
    within_factors = NA,
    contrast_names = names(between_contrasts)
  )
}

make_within_fixture <- function() {
  data("priming", package = "PsyR", envir = environment())

  required_columns <- c("subj", "Condition", "RT")
  missing_columns <- setdiff(required_columns, names(priming))
  if (length(missing_columns) > 0) {
    stop(
      "The priming dataset is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  priming$Condition <- factor(
    priming$Condition,
    levels = c("STEREO", "ATYPICAL", "NEUTRAL")
  )

  afex::afex_options(emmeans_model = "multivariate")
  model <- afex::aov_ez(
    id = "subj",
    dv = "RT",
    data = priming,
    within = "Condition"
  )

  emmeans <- emmeans::emmeans(model, "Condition")
  within_contrasts <- list(
    "Stereo - Neutral" = c(1, 0, -1),
    "Atypical - Neutral" = c(0, 1, -1)
  )
  contrast_table <- emmeans::contrast(
    emmeans,
    method = within_contrasts,
    adjust = "none"
  )

  list(
    data = priming,
    model = model,
    emmeans = emmeans,
    tables = list(contrast_table),
    families = list("w"),
    between_factors = NA,
    within_factors = list("Condition"),
    contrast_names = names(within_contrasts)
  )
}

make_mixed_fixture <- function() {
  data("spacing", package = "PsyR", envir = environment())

  required_columns <- c("subj", "group", "spacing", "yield")
  missing_columns <- setdiff(required_columns, names(spacing))
  if (length(missing_columns) > 0) {
    stop(
      "The spacing dataset is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  spacing$group <- factor(spacing$group)
  spacing$spacing <- factor(
    spacing$spacing,
    levels = c("TWENTY", "FORTY", "SIXTY")
  )

  afex::afex_options(emmeans_model = "multivariate")
  model <- afex::aov_ez(
    id = "subj",
    dv = "yield",
    data = spacing,
    between = "group",
    within = "spacing"
  )

  between_emmeans <- emmeans::emmeans(model, "group")
  within_emmeans <- emmeans::emmeans(model, "spacing")
  interaction_emmeans <- emmeans::emmeans(model, c("group", "spacing"))

  between_contrasts <- list(
    "Group 1 - Group 2" = c(1, -1, 0, 0),
    "Group 3 - Group 4" = c(0, 0, 1, -1)
  )
  within_contrasts <- list(
    "Twenty - Forty" = c(1, -1, 0),
    "Sixty - Forty" = c(0, -1, 1)
  )

  between_table <- emmeans::contrast(
    between_emmeans,
    method = between_contrasts,
    adjust = "none"
  )
  within_table <- emmeans::contrast(
    within_emmeans,
    method = within_contrasts,
    adjust = "none"
  )
  interaction_table <- emmeans::contrast(
    interaction_emmeans,
    interaction = list(between_contrasts, within_contrasts),
    adjust = "none"
  )

  list(
    data = spacing,
    model = model,
    emmeans = list(
      between = between_emmeans,
      within = within_emmeans,
      interaction = interaction_emmeans
    ),
    tables = list(between_table, within_table, interaction_table),
    families = list("b", "w", "bw"),
    between_factors = list("group"),
    within_factors = list("spacing"),
    contrast_names = list(
      between = names(between_contrasts),
      within = names(within_contrasts),
      interaction = names(interaction_table)
    )
  )
}

expect_psyci_table <- function(result, n_rows) {
  expect_true(is.data.frame(result))
  expect_true(all(c(
    "estimate", "SE", "df", "cc", "lower", "upper"
  ) %in% names(result)))
  expect_equal(nrow(result), n_rows)
  expect_true(all(is.finite(result$estimate)))
  expect_true(all(is.finite(result$SE)))
  expect_true(all(is.finite(result$df)))
  expect_true(all(is.finite(result$cc)))
  expect_true(all(is.finite(result$lower)))
  expect_true(all(is.finite(result$upper)))
  expect_true(all(result$lower <= result$upper))
  expect_equal(
    result$lower,
    result$estimate - result$cc * result$SE,
    tolerance = 1e-12
  )
  expect_equal(
    result$upper,
    result$estimate + result$cc * result$SE,
    tolerance = 1e-12
  )
  invisible(result)
}
