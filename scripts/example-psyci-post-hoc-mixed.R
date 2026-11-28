# Post-hoc PsyR analysis for a mixed between-by-within design.
#
# This is a standalone version of the analysis in
# vignettes/post-hoc_mixed-design.Rmd. Run it from the PsyR repository root.
# The script keeps `post_hoc_results` available after source() returns.

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

post_hoc_results <- local({
    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    afex::afex_options(emmeans_model = "multivariate")

    data("social_anxiety", package = "PsyR", envir = environment())
    social_anxiety$Session <- factor(
        social_anxiety$Session,
        levels = c("Pre", "Post", "FU1", "FU2")
    )
    social_anxiety$Group <- factor(
        social_anxiety$Group,
        levels = c(1, 2, 3),
        labels = c("New", "Standard", "MinContact")
    )

    model <- afex::aov_ez(
        id = "Subject",
        dv = "Score",
        data = social_anxiety,
        between = "Group",
        within = "Session"
    )

    # Define the intended comparisons in convenient integer form, then apply
    # Psy's mean-difference scaling before passing them to emmeans.
    raw_between_contrasts <- list(
        "Ts - C" = c(1, 1, -2),
        "NT - ST" = c(1, -1, 0)
    )
    between_contrasts <- rescale_contrasts(
        raw_between_contrasts,
        mode = "mean_difference"
    )

    raw_within_contrasts <- list(
        "Pre - Rest" = c(3, -1, -1, -1),
        "Post - FUs" = c(0, 2, -1, -1),
        "FU1 - FU2" = c(0, 0, 1, -1)
    )
    within_contrasts <- rescale_contrasts(
        raw_within_contrasts,
        mode = "mean_difference"
    )

    between_emmeans <- emmeans::emmeans(model, "Group")
    within_emmeans <- emmeans::emmeans(model, "Session")
    cell_emmeans <- emmeans::emmeans(model, c("Group", "Session"))

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

    # Each component above is an order-0 mean difference. emmeans forms the
    # correctly scaled order-1 product, so it must not be rescaled again.
    interaction_table <- emmeans::contrast(
        cell_emmeans,
        interaction = list(between_contrasts, within_contrasts),
        adjust = "none"
    )

    contrast_tables <- list(
        within = within_table,
        between = between_table,
        between_by_within = interaction_table
    )
    alpha <- 0.05
    results <- psyci(
        model = model,
        contrast_tables = contrast_tables,
        method = "ph",
        alpha = alpha
    )

    expected_methods <- c("Hotelling T-squared", "Scheffe", "Roy GCR")

    for (index in seq_along(results)) {
        table <- results[[index]]
        available_p_values <- !is.na(table$p_adj)
        comparable <- stats::complete.cases(
            table$p_adj,
            table$lower,
            table$upper
        )

        stopifnot(
            all(c(
                "df1", "df2", "F_value", "p_value",
                "p_adj", "method"
            ) %in% names(table)),
            !any(c("df", "t.ratio", "p.value") %in% names(table)),
            all(table$df1 == 1L),
            isTRUE(all.equal(
                table$p_value,
                stats::pf(
                    table$F_value, 1, table$df2,
                    lower.tail = FALSE
                )
            )),
            all(is.finite(table$p_adj[available_p_values])),
            all(table$p_adj[available_p_values] >= 0),
            all(table$p_adj[available_p_values] <= 1),
            all(table$method == expected_methods[[index]]),
            any(grepl(
                "PsyR simultaneous-inference procedure:",
                attr(table, "mesg"),
                fixed = TRUE
            )),
            all(
                (table$p_adj[comparable] <= alpha) ==
                    (table$lower[comparable] > 0 |
                        table$upper[comparable] < 0)
            )
        )

        cat(
            "\n",
            names(results)[[index]],
            " family — ",
            unique(table$method),
            "\n",
            sep = ""
        )
        print(as.data.frame(table))
        cat("Messages:\n", paste(attr(table, "mesg"), collapse = "\n"), "\n")
    }

    gcr_parameters <- attr(
        results[[3L]],
        "psyr_probability_parameters"
    )
    stopifnot(
        identical(names(gcr_parameters), c("s", "m", "n", "theta")),
        isTRUE(all.equal(gcr_parameters$s, 2)),
        isTRUE(all.equal(gcr_parameters$m, 0)),
        isTRUE(all.equal(gcr_parameters$n, 41.5)),
        !"theta" %in% names(results[[3L]]),
        isTRUE(all.equal(
            gcr_parameters$theta,
            results[[3L]]$F_value /
                (results[[3L]]$df2 + results[[3L]]$F_value)
        ))
    )

    cat("\nAll mixed post-hoc p-value and CI checks passed.\n")
    results
})
