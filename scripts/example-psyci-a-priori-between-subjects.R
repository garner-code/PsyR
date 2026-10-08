# A-priori Bonferroni PsyR analysis for an all-between-subjects J x K design.
#
# This is a standalone version of the planned analysis in
# vignettes/J-K-Design_Analysis-allowing-for-inferences-on-simple-effect-contrasts.Rmd.
# Run it from the PsyR repository root. The script keeps
# `a_priori_between_results` available after source() returns.

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

a_priori_between_results <- local({
    data("experience", package = "PsyR", envir = environment())
    experience$Experience <- factor(
        experience$Experience,
        levels = c("10yrs", "3yrs", "<1yr")
    )
    experience$Therapy <- factor(
        experience$Therapy,
        levels = c(
            "Behavioural",
            "Cognitive",
            "Psychotherapy",
            "Group"
        )
    )

    model <- afex::aov_ez(
        id = "SubjID",
        dv = "Outcome",
        data = experience,
        between = c("Experience", "Therapy")
    )

    raw_experience_contrasts <- list(
        "lots_vs_little" = c(0, 1, -1),
        "lots_vs_some" = c(1, -1, 0)
    )
    experience_contrasts <- rescale_contrasts(
        raw_experience_contrasts,
        mode = "mean_difference"
    )

    raw_therapy_contrasts <- list(
        "evidence_vs_woowoo" = c(1, 1, -1, -1),
        "behav_vs_cog" = c(1, -1, 0, 0),
        "woo1_vs_woo2" = c(0, 0, 1, -1)
    )
    therapy_contrasts <- rescale_contrasts(
        raw_therapy_contrasts,
        mode = "mean_difference"
    )

    cell_emmeans <- emmeans::emmeans(
        model,
        c("Experience", "Therapy")
    )
    experience_emmeans <- emmeans::emmeans(model, "Experience")
    experience_main_table <- emmeans::contrast(
        experience_emmeans,
        method = experience_contrasts,
        adjust = "bonferroni"
    )

    experience_by_therapy <- emmeans::emmeans(
        model,
        ~ Experience | Therapy
    )
    experience_simple_table <- emmeans::contrast(
        experience_by_therapy,
        method = experience_contrasts,
        adjust = "bonferroni"
    )

    # The component lists are already mean scaled. Their emmeans product is a
    # correctly scaled order-1 interaction and must not be rescaled again.
    interaction_table <- emmeans::contrast(
        cell_emmeans,
        interaction = list(experience_contrasts, therapy_contrasts),
        adjust = "bonferroni"
    )

    contrast_tables <- list(
        experience_main_effects = experience_main_table,
        experience_simple_effects = experience_simple_table,
        interactions = interaction_table
    )
    original_p_values <- lapply(
        contrast_tables,
        function(table) as.data.frame(summary(table))$p.value
    )

    # The vignette treats these three tables as one non-independent family of
    # 16 planned contrasts and uses alpha = 2 * .05 for that family.
    alpha <- 0.10
    results <- psyci(
        model = model,
        contrast_tables = contrast_tables,
        method = "bf",
        alpha = rep(list(alpha), length(contrast_tables)),
        independent = FALSE
    )

    expected_method <- "Bonferroni"
    for (index in seq_along(results)) {
        table <- results[[index]]
        available_p_values <- !is.na(table$psyr_p_value)
        comparable <- stats::complete.cases(
            table$psyr_p_value,
            table$lower,
            table$upper
        )

        stopifnot(
            "p.value" %in% names(table),
            "psyr_p_value" %in% names(table),
            isTRUE(all.equal(
                table$p.value,
                original_p_values[[index]],
                check.attributes = FALSE
            )),
            all(is.finite(table$psyr_p_value[available_p_values])),
            all(table$psyr_p_value[available_p_values] >= 0),
            all(table$psyr_p_value[available_p_values] <= 1),
            identical(
                attr(table, "psyr_p_value_method"),
                expected_method
            ),
            any(startsWith(
                attr(table, "mesg"),
                "PsyR p-value method:"
            )),
            all(
                (table$psyr_p_value[comparable] <= alpha) ==
                    (table$lower[comparable] > 0 |
                        table$upper[comparable] < 0)
            )
        )

        cat(
            "\nTable ",
            index,
            " — ",
            attr(table, "psyr_p_value_method"),
            "\n",
            sep = ""
        )
        print(as.data.frame(table))
        cat("Messages:\n", paste(attr(table, "mesg"), collapse = "\n"), "\n")
    }

    cat("\nAll a-priori between-subjects p-value and CI checks passed.\n")
    results
})
