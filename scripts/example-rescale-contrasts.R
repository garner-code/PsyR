# Runnable example of rescale_contrasts() with emmeans.
#
# Run this script from the PsyR repository root. It demonstrates both supported
# scaling modes and verifies two equivalent ways to specify an interaction.

if (!file.exists("DESCRIPTION")) {
    stop(
        "Run this script from the PsyR repository root.",
        call. = FALSE
    )
}

required_packages <- c("devtools", "afex", "emmeans")
missing_packages <- required_packages[!vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
)]
if (length(missing_packages) > 0L) {
    stop(
        sprintf(
            "Install the required package(s) before running this example: %s.",
            paste(missing_packages, collapse = ", ")
        ),
        call. = FALSE
    )
}

devtools::load_all(".", quiet = TRUE)

local({
    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    afex::afex_options(emmeans_model = "multivariate")

    data("spacing", package = "PsyR", envir = environment())
    spacing$group <- factor(spacing$group, levels = 1:4)
    spacing$spacing <- factor(
        spacing$spacing,
        levels = c("TWENTY", "FORTY", "SIXTY")
    )

    model <- afex::aov_ez(
        id = "subj",
        dv = "yield",
        data = spacing,
        between = "group",
        within = "spacing"
    )

    # The coefficient order follows the factor levels set above. These raw
    # vectors express the intended comparisons but are not mean scaled.
    raw_between <- list(
        "Groups 1 and 2 versus groups 3 and 4" = c(1, 1, -1, -1)
    )
    raw_within <- list(
        "Twenty versus the mean of forty and sixty" = c(2, -1, -1)
    )

    between_contrasts <- rescale_contrasts(
        raw_between,
        mode = "mean_difference"
    )
    within_contrasts <- rescale_contrasts(
        raw_within,
        mode = "mean_difference"
    )

    cat("\nMean-difference-scaled component contrasts:\n")
    print(between_contrasts)
    print(within_contrasts)

    between_emmeans <- emmeans::emmeans(model, "group")
    within_emmeans <- emmeans::emmeans(model, "spacing")
    cell_emmeans <- emmeans::emmeans(model, c("group", "spacing"))

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

    cat("\nBetween-subjects contrast table:\n")
    print(between_table)
    cat("\nWithin-subjects contrast table:\n")
    print(within_table)

    # emmeans multiplies the two supplied component vectors. Each component is
    # an order-0 mean difference, so their product is already a correctly
    # scaled order-1 interaction. Do not apply interaction scaling to either
    # component before passing it to interaction = list(...).
    crossed_interaction_table <- emmeans::contrast(
        cell_emmeans,
        interaction = list(between_contrasts, within_contrasts),
        adjust = "none"
    )

    cat("\nInteraction formed by crossing the scaled components:\n")
    print(crossed_interaction_table)

    # Interaction mode is used when the complete product coefficient vector
    # has already been formed. outer() produces coefficients in the same order
    # as this emmeans cell grid: group changes fastest, then spacing.
    raw_product <- list(
        "Group by spacing interaction" = as.vector(outer(
            raw_between[[1L]],
            raw_within[[1L]]
        ))
    )
    direct_interaction <- rescale_contrasts(
        raw_product,
        mode = "interaction",
        interaction_order = 1L
    )

    # The directly scaled product must equal the product of the two separately
    # mean-scaled component vectors.
    expected_product <- as.vector(outer(
        between_contrasts[[1L]],
        within_contrasts[[1L]]
    ))
    stopifnot(isTRUE(all.equal(
        direct_interaction[[1L]],
        expected_product,
        check.attributes = FALSE
    )))

    direct_interaction_table <- emmeans::contrast(
        cell_emmeans,
        method = direct_interaction,
        adjust = "none"
    )

    cat("\nThe same interaction supplied as one preformed product vector:\n")
    print(direct_interaction_table)

    crossed_summary <- as.data.frame(summary(crossed_interaction_table))
    direct_summary <- as.data.frame(summary(direct_interaction_table))
    comparison_columns <- c("estimate", "SE", "df", "t.ratio", "p.value")
    stopifnot(isTRUE(all.equal(
        crossed_summary[comparison_columns],
        direct_summary[comparison_columns],
        check.attributes = FALSE
    )))

    cat("\nBoth interaction specifications give the same results.\n")
})
