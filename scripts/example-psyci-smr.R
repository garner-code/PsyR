# Studentized Maximum Root (SMR) PsyR analysis for a 2 x 2 design.
#
# This is a standalone version of the analysis in
# vignettes/between_subjects_2x2_factorial-w-simple.Rmd. Run it from the PsyR
# repository root. The script keeps `smr_results` available after source()
# returns.

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

smr_results <- local({
    data("drive_sleepy", package = "PsyR", envir = environment())
    drive_sleepy$sleep_deprivation <- factor(
        drive_sleepy$sleep_deprivation,
        levels = c(12, 0),
        labels = c("High", "None")
    )
    drive_sleepy$nvh <- factor(
        drive_sleepy$nvh,
        levels = c("high", "low")
    )

    model <- afex::aov_ez(
        id = "subID",
        dv = "performance",
        data = drive_sleepy,
        between = c("sleep_deprivation", "nvh")
    )
    cell_emmeans <- emmeans::emmeans(
        model,
        c("sleep_deprivation", "nvh")
    )

    raw_main_effects <- list(
        "Sleep_Dep" = c(1, -1, 1, -1),
        "NVH" = c(1, 1, -1, -1)
    )
    main_effects <- rescale_contrasts(
        raw_main_effects,
        mode = "mean_difference"
    )
    main_table <- emmeans::contrast(
        cell_emmeans,
        method = main_effects,
        adjust = "none"
    )

    # This is a preformed two-factor product vector. Psy interaction scaling
    # changes the vignette's c(.5, -.5, -.5, .5) to c(1, -1, -1, 1).
    raw_interaction <- list(
        "Sleep deprivation x NVH" = c(0.5, -0.5, -0.5, 0.5)
    )
    interaction_contrast <- rescale_contrasts(
        raw_interaction,
        mode = "interaction",
        interaction_order = 1L
    )
    stopifnot(isTRUE(all.equal(
        interaction_contrast[[1L]],
        c(1, -1, -1, 1),
        check.attributes = FALSE
    )))
    interaction_table <- emmeans::contrast(
        cell_emmeans,
        method = interaction_contrast,
        adjust = "none"
    )

    raw_nvh_simple <- list("NVH" = c(1, -1))
    nvh_simple_contrast <- rescale_contrasts(
        raw_nvh_simple,
        mode = "mean_difference"
    )
    nvh_by_sleep_deprivation <- emmeans::emmeans(
        model,
        ~ nvh | sleep_deprivation
    )
    nvh_simple_table <- emmeans::contrast(
        nvh_by_sleep_deprivation,
        method = nvh_simple_contrast,
        adjust = "none"
    )

    raw_sleep_simple <- list("Sleep deprivation" = c(1, -1))
    sleep_simple_contrast <- rescale_contrasts(
        raw_sleep_simple,
        mode = "mean_difference"
    )
    sleep_deprivation_by_nvh <- emmeans::emmeans(
        model,
        ~ sleep_deprivation | nvh
    )
    sleep_simple_table <- emmeans::contrast(
        sleep_deprivation_by_nvh,
        method = sleep_simple_contrast,
        adjust = "none"
    )

    contrast_tables <- list(
        main_effects = main_table,
        interaction = interaction_table,
        nvh_simple_effects = nvh_simple_table,
        sleep_deprivation_simple_effects = sleep_simple_table
    )
    original_p_values <- lapply(
        contrast_tables,
        function(table) as.data.frame(summary(table))$p.value
    )

    alpha <- 0.05
    smr_seed <- 42L
    smr_parameters <- list(
        p = 2L,
        q = 2L,
        n_sim = 100000L
    )

    # psyci() fits one seeded SMR reference distribution here and reuses it for
    # the confidence intervals and p-values in every contrast table.
    results <- psyci(
        model = model,
        contrast_tables = contrast_tables,
        method = "smr",
        smr_params = smr_parameters,
        seed = smr_seed,
        alpha = alpha
    )

    expected_method <- "Studentized Maximum Root (SMR)"
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

    cat(
        "\nAll SMR p-value and CI checks passed using ",
        smr_parameters$n_sim,
        " simulations and seed ",
        smr_seed,
        " from one shared SMR reference distribution",
        ".\n",
        sep = ""
    )
    results
})
