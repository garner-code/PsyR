# Run the vignette-derived PsyR analyses and compare their numeric results with
# an independent port of the original Psy Pascal calculations.

source("scripts/legacy-psy-reference.R")
devtools::load_all(quiet = TRUE)

run_analysis_script <- function(file) {
    environment <- new.env(parent = globalenv())
    sys.source(file, envir = environment)
    environment
}

flatten_psyr <- function(result) {
    if (!is.list(result) || inherits(result, "data.frame")) result <- list(result)
    tables <- lapply(seq_along(result), function(i) {
        table <- as.data.frame(result[[i]])
        table$.table <- i
        table
    })
    do.call(rbind, tables)
}

compare_result <- function(analysis, psyr, legacy, note = "") {
    psyr <- flatten_psyr(psyr)
    if (nrow(psyr) != nrow(legacy)) {
        stop(sprintf("%s: PsyR has %d rows but the legacy port has %d.",
            analysis, nrow(psyr), nrow(legacy)))
    }
    columns <- c("estimate", "SE", "cc", "lower", "upper")
    differences <- vapply(columns, function(column) {
        max(abs(psyr[[column]] - legacy[[column]]), na.rm = TRUE)
    }, numeric(1))
    comparison_details[[analysis]] <<- data.frame(
        row = seq_len(nrow(legacy)),
        contrast = legacy$contrast,
        family = legacy$family,
        legacy_estimate = legacy$estimate,
        psyr_estimate = psyr$estimate,
        estimate_difference = psyr$estimate - legacy$estimate,
        legacy_SE = legacy$SE,
        psyr_SE = psyr$SE,
        SE_difference = psyr$SE - legacy$SE,
        legacy_df = legacy$df,
        psyr_df = psyr$df,
        legacy_cc = legacy$cc,
        psyr_cc = psyr$cc,
        cc_difference = psyr$cc - legacy$cc,
        legacy_lower = legacy$lower,
        psyr_lower = psyr$lower,
        lower_difference = psyr$lower - legacy$lower,
        legacy_upper = legacy$upper,
        psyr_upper = psyr$upper,
        upper_difference = psyr$upper - legacy$upper,
        stringsAsFactors = FALSE
    )
    data.frame(
        analysis = analysis,
        rows = nrow(legacy),
        core_matches = differences[["estimate"]] < 1e-10 &&
            differences[["SE"]] < 1e-10,
        max_estimate_difference = differences[["estimate"]],
        max_se_difference = differences[["SE"]],
        max_cc_difference = differences[["cc"]],
        max_limit_difference = max(differences[c("lower", "upper")]),
        note = note,
        stringsAsFactors = FALSE
    )
}

make_item <- function(label, family, a = NULL, w = NULL) {
    list(label = label, family = family, a = a, w = w)
}

results <- list()
comparison_details <- list()

# Analysis of a 3-between x 3-within design -------------------------------
analysis <- run_analysis_script("scripts/Analysis-of-3-between-x-3-within-design.R")
design <- legacy_psy_design(
    analysis$depression, "Subject", "Happiness", "Group", "Time"
)
a <- list(
    Ts_v_Ctrl = c(-1, 0.5, 0.5), T1_v_T2 = c(0, 1, -1),
    `T1*` = c(0, 1, 0), `T2*` = c(0, 0, 1), `Ctrl*` = c(1, 0, 0)
)
w <- list(
    Post_vs_Pre = c(-1, 1, 0), FU_vs_Post = c(0, -1, 1),
    `Post*` = c(0, 1, 0), `FU*` = c(0, 0, 1)
)
items <- c(
    legacy_psy_cross(a, w[1:2], FALSE, FALSE, TRUE),
    legacy_psy_cross(a[1:2], w[3:4], FALSE, FALSE, TRUE),
    legacy_psy_cross(a[1:2], list(), TRUE, FALSE, FALSE),
    legacy_psy_cross(list(), w[1:2], FALSE, TRUE, FALSE)
)
legacy <- legacy_psy_analyse(
    design, items, "bf", alpha = 0.15, n_k = 18L
)
results[[length(results) + 1L]] <- compare_result(
    "3-between x 3-within", analysis$contrasts_w_cis, legacy,
    paste0("One 18-contrast all-factorial family; selector-by-selector cell ",
        "means are dropped before correction.")
)

# Between-subjects 2 x 2 simple effects, SMR ------------------------------
analysis <- run_analysis_script("scripts/between_subjects_2x2_factorial-w-simple.R")
design <- legacy_psy_design(
    analysis$drive_sleepy, "subID", "performance",
    c("sleep_deprivation", "nvh")
)
items <- list(
    make_item("Sleep_Dep", "b", c(0.5, -0.5, 0.5, -0.5)),
    make_item("NVH", "b", c(0.5, 0.5, -0.5, -0.5)),
    make_item("int", "b", c(1, -1, -1, 1)),
    make_item("NVH | High", "b", c(1, 0, -1, 0)),
    make_item("NVH | None", "b", c(0, 1, 0, -1)),
    make_item("Sleep | high", "b", c(1, -1, 0, 0)),
    make_item("Sleep | low", "b", c(0, 0, 1, -1))
)
legacy <- legacy_psy_analyse(
    design, items, "smr", alpha = 0.05, smr = list(p = 2L, q = 2L)
)
results[[length(results) + 1L]] <- compare_result(
    "between-subjects 2 x 2 SMR", analysis$contrasts_w_cis, legacy,
    paste0("Psy uses deterministic quadrature; PsyR uses the vignette's ",
        "seeded 100,000-draw simulation, so its critical constant is approximate.")
)

# J x K design -------------------------------------------------------------
analysis <- run_analysis_script(
    "scripts/J-K-Design_Analysis-allowing-for-inferences-on-simple-effect-contrasts.R"
)
design <- legacy_psy_design(
    analysis$experience, "SubjID", "Outcome", c("Experience", "Therapy")
)
experience_c <- list(
    lots_vs_little = c(0, 1, -1), lots_vs_some = c(1, -1, 0)
)
therapy_c <- list(
    evidence_vs_woowoo = c(0.5, 0.5, -0.5, -0.5),
    behav_vs_cog = c(1, -1, 0, 0), woo1_vs_woo2 = c(0, 0, 1, -1)
)
exp_main <- lapply(names(experience_c), function(name) {
    make_item(name, "b", rep(experience_c[[name]] / 4, 4))
})
exp_simple <- list()
for (therapy in seq_len(4)) for (name in names(experience_c)) {
    coefficient <- numeric(12)
    coefficient[(3 * therapy - 2):(3 * therapy)] <- experience_c[[name]]
    exp_simple[[length(exp_simple) + 1L]] <- make_item(
        paste(name, "therapy", therapy), "b", coefficient
    )
}
interactions <- list()
for (therapy_name in names(therapy_c)) for (exp_name in names(experience_c)) {
    interactions[[length(interactions) + 1L]] <- make_item(
        paste(exp_name, therapy_name), "b",
        as.vector(outer(experience_c[[exp_name]], therapy_c[[therapy_name]]))
    )
}
planned_items <- c(exp_main, exp_simple, interactions)
legacy <- legacy_psy_analyse(
    design, planned_items, "bf", alpha = 0.10, n_k = 16L
)
results[[length(results) + 1L]] <- compare_result(
    "J x K pooled planned", analysis$planned_results, legacy
)
legacy <- legacy_psy_analyse(
    design, planned_items, "ph", alpha = 1 - 0.95^2
)
results[[length(results) + 1L]] <- compare_result(
    "J x K pooled post-hoc", analysis$post_hoc_results, legacy,
    paste0("Expected critical-constant difference: original Psy uses all 12 cell ",
        "groups (11 df); the PsyR workflow explicitly supplies nu1 = 8.")
)
therapy_main <- lapply(names(therapy_c), function(name) {
    make_item(name, "b", as.vector(outer(rep(1 / 3, 3), therapy_c[[name]])))
})
standard_items <- c(exp_main, therapy_main, interactions)
legacy <- legacy_psy_analyse(
    design, standard_items, "bf", alpha = 0.05,
    n_k = c(rep(2L, 2), rep(3L, 3), rep(6L, 6))
)
results[[length(results) + 1L]] <- compare_result(
    "J x K standard planned", analysis$standard_results, legacy
)

# Planned between, within and interaction contrasts -----------------------
analysis <- run_analysis_script("scripts/planned-contrasts_between-within.R")
design <- legacy_psy_design(
    analysis$depression, "Subject", "Happiness", "Group", "Time"
)
a <- list(Ts_v_Ctrl = c(-1, 0.5, 0.5), T1_v_T2 = c(0, 1, -1))
w <- list(Post_vs_Pre = c(-1, 1, 0), FU_vs_Post = c(0, -1, 1))
items <- legacy_psy_cross(a, w)
legacy <- legacy_psy_analyse(
    design, items, "bf", alpha = 0.05,
    n_k = c(rep(2L, 2), rep(2L, 2), rep(4L, 4))
)
results[[length(results) + 1L]] <- compare_result(
    "planned between-within", analysis$contrasts_w_cis, legacy
)

# Planned factorial mixed-design simple effects ---------------------------
analysis <- run_analysis_script(
    "scripts/planned-factorial_mixed-designs_simple-effects.R"
)
design <- legacy_psy_design(
    analysis$depression, "Subject", "Happiness", "Group", "Time"
)
a <- list(Ts_v_Ctrl = c(-1, 0.5, 0.5), T1_v_T2 = c(0, 1, -1))
w <- list(Post_vs_Pre = c(-1, 1, 0), FU_vs_Post = c(0, -1, 1))
items <- c(
    legacy_psy_cross(list(), w, FALSE, TRUE, FALSE),
    legacy_psy_cross(a, list(), TRUE, FALSE, FALSE),
    legacy_psy_cross(a, w, FALSE, FALSE, TRUE)
)
for (group in seq_len(3)) for (name in names(w)) {
    selector <- numeric(3)
    selector[group] <- 1
    items[[length(items) + 1L]] <- make_item(
        paste(name, "group", group), "bw", selector, w[[name]]
    )
}
within_selectors <- list(`Post*` = c(0, 1, 0), `FU*` = c(0, 0, 1))
items <- c(items, legacy_psy_cross(a, within_selectors, FALSE, FALSE, TRUE))
legacy <- legacy_psy_analyse(
    design, items, "bf", alpha = 0.15, n_k = 18L
)
results[[length(results) + 1L]] <- compare_result(
    "planned mixed simple effects", analysis$contrasts_w_cis, legacy,
    paste0("Psy's user-supplied-constant route is used for the selected pooled ",
        "family; its standard interface would generate unwanted cross-products.")
)

# Post-hoc mixed design ----------------------------------------------------
analysis <- run_analysis_script("scripts/post-hoc_mixed-design.R")
design <- legacy_psy_design(
    analysis$social_anxiety, "Subject", "Score", "Group", "Session"
)
a <- list(`Ts - C` = c(0.5, 0.5, -1), `NT - ST` = c(1, -1, 0))
w <- list(
    `Pre - Rest` = c(1, -1 / 3, -1 / 3, -1 / 3),
    `Post - FUs` = c(0, 1, -0.5, -0.5), `FU1 - FU2` = c(0, 0, 1, -1)
)
items <- c(
    legacy_psy_cross(list(), w, FALSE, TRUE, FALSE),
    legacy_psy_cross(a, list(), TRUE, FALSE, FALSE),
    legacy_psy_cross(a, w, FALSE, FALSE, TRUE)
)
legacy <- legacy_psy_analyse(
    design, items, "bf", alpha = 0.05,
    n_k = c(rep(3L, 3), rep(2L, 2), rep(6L, 6))
)
results[[length(results) + 1L]] <- compare_result(
    "mixed design planned", analysis$planned_results, legacy
)
legacy <- legacy_psy_analyse(design, items, "ph", alpha = 0.05)
results[[length(results) + 1L]] <- compare_result(
    "mixed design post-hoc", analysis$post_hoc_results, legacy
)

comparison_summary <- do.call(rbind, results)
print(comparison_summary, row.names = FALSE, digits = 8)

report_file <- "LEGACY_PSY_COMPARISON_RESULTS.txt"
report <- c(
    "PsyR versus original Psy source calculations",
    "==============================================",
    "",
    paste0(
        "This report compares PsyR with an independent R translation of the ",
        "calculation routines in PSY-source_no-touchy. Differences are shown ",
        "as PsyR minus legacy Psy."
    ),
    "",
    "SUMMARY",
    "-------",
    capture.output(print(
        comparison_summary[setdiff(names(comparison_summary), "note")],
        row.names = FALSE,
        digits = 12
    )),
    "",
    "INTERPRETATION NOTES",
    "--------------------",
    unlist(lapply(seq_len(nrow(comparison_summary)), function(i) {
        note <- comparison_summary$note[[i]]
        if (!nzchar(note)) return(NULL)
        paste0("- ", comparison_summary$analysis[[i]], ": ", note)
    })),
    ""
)
for (analysis_name in names(comparison_details)) {
    detail <- comparison_details[[analysis_name]]
    core <- detail[c(
        "row", "contrast", "family", "legacy_estimate", "psyr_estimate",
        "estimate_difference", "legacy_SE", "psyr_SE", "SE_difference",
        "legacy_df", "psyr_df"
    )]
    critical <- detail[c(
        "row", "contrast", "legacy_cc", "psyr_cc", "cc_difference"
    )]
    limits <- detail[c(
        "row", "contrast", "legacy_lower", "psyr_lower", "lower_difference",
        "legacy_upper", "psyr_upper", "upper_difference"
    )]
    report <- c(
        report,
        toupper(analysis_name),
        paste(rep("-", nchar(analysis_name)), collapse = ""),
        "",
        "Estimates, standard errors, and error degrees of freedom",
        capture.output(print(
            core,
            row.names = FALSE,
            digits = 12,
            width = 220
        )),
        "",
        "Critical constants",
        capture.output(print(
            critical,
            row.names = FALSE,
            digits = 12,
            width = 220
        )),
        "",
        "Confidence limits",
        capture.output(print(
            limits,
            row.names = FALSE,
            digits = 12,
            width = 220
        )),
        ""
    )
}
writeLines(report, report_file, useBytes = TRUE)
message("Detailed comparison written to ", normalizePath(report_file))

invisible(comparison_summary)
