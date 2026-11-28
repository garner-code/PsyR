# Formatting psyci() output for a vignette.
#
# Follows vignettes/planned-contrasts_between-within.Rmd, using the same
# depression data, factor order, planned contrasts, and Bonferroni procedure.
# Run from the repository root:
# source("scripts/example-psyci-vignette-tables.R")
#
# Objects remain available afterwards. In particular:
#   output                  Original numerical results (full precision).
#   printed_output          Exact console output, as a character vector.
#   raw_display             One display table containing all three families.
#   standardized_display    The corresponding table in sample SD units.
#   raw_flextable           Styled table for a vignette.
#   standardized_flextable  Styled standardized table for a vignette.
#   table_preview_path      HTML preview containing both styled tables.

# 1. Load the development package and the presentation dependencies ----------
if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the PsyR repository root.", call. = FALSE)
}

required_packages <- c("devtools", "afex", "emmeans", "dplyr", "knitr", "flextable")
missing_packages <- required_packages[!vapply(
    required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
    stop(
        "Install the required packages: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}
devtools::load_all(".", quiet = TRUE)

# 2. Reproduce the vignette's model and planned contrasts -------------------
# local() restores the user's afex setting when the analysis finishes.
output <- local({
    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(afex::afex_options(emmeans_model = old_emmeans_model))
    afex::afex_options(emmeans_model = "multivariate")

    data("depression", package = "PsyR", envir = environment())
    depression$Group <- factor(depression$Group, c("Ctrl", "Tmt1", "Tmt2"))
    depression$Time <- factor(depression$Time, c("Pre", "Post", "FU"))
    depression$Subject <- factor(depression$Subject)

    model <- afex::aov_ez(
        id = "Subject",
        dv = "Happiness",
        data = depression,
        between = "Group",
        within = "Time",
        anova_table = list(correction = "none")
    )

    between_contrasts <- list(
        "Ts - Ctrl" = c(-1, 0.5, 0.5),
        "T1 - T2" = c(0, 1, -1)
    )
    within_contrasts <- list(
        "Post - Pre" = c(-1, 1, 0),
        "FU - Post" = c(0, -1, 1)
    )

    between_table <- emmeans::contrast(
        emmeans::emmeans(model, "Group"), method = between_contrasts
    )
    within_table <- emmeans::contrast(
        emmeans::emmeans(model, "Time"), method = within_contrasts
    )
    interaction_table <- emmeans::contrast(
        emmeans::emmeans(model, c("Group", "Time")),
        interaction = list(between_contrasts, within_contrasts)
    )

    PsyR::psyci(
        model = model,
        contrast_tables = list(between_table, within_table, interaction_table),
        method = "bf"
    )
})

# 3. Capture the actual console output -------------------------------------
# Produces headings, raw and standardized tables, and the default notes.
# capture.output() returns text; print(output) itself invisibly returns output,
# so assigning printed_output <- print(output) would still give you the list.
printed_output <- capture.output(print(output))
cat(printed_output, sep = "\n")

# In an R Markdown chunk, the preceding cat() call gives console-style output.
# Use capture.output(print(output, verbose = TRUE)) to retain detailed notes.

# 4. Reuse the printed table's formatting as a data frame -------------------
# This small example helper uses the SAME internal formatter as print.psyr_ci:
# three decimal places, blank missing values, combined interaction labels,
# and hidden p_value/polynomial columns. It also combines contrast families.
#
# The ::: access is intentional: PsyR currently has no public display-table
# accessor. This example may need updating if that internal helper changes.
# All display columns are character strings; use output for calculations.
combine_display_tables <- function(tables) {
    displays <- lapply(tables, PsyR:::.psyci_display_table)
    display <- dplyr::bind_rows(displays, .id = "family")
    labels <- c(b = "between", w = "within", bw = "mixed")
    matched <- match(display$family, names(labels))
    known <- !is.na(matched)
    display$family[known] <- unname(labels[matched[known]])
    # Different families can have different columns. The console printer
    # fills those absent cells with blanks too.
    display[] <- lapply(display, function(column) {
        column[is.na(column)] <- ""
        column
    })
    display
}

raw_display <- combine_display_tables(output)
standardized_tables <- attr(output, "scaled_tables", exact = TRUE)
if (is.null(standardized_tables)) {
    stop("Standardized results are unavailable for this model.", call. = FALSE)
}
standardized_display <- combine_display_tables(standardized_tables)

# 5. A simple knitr table --------------------------------------------------
# Produces a single eight-row table: two between, two within, four mixed.
# Here print() shows Markdown table text in the console. In an R Markdown
# chunk, use knitr::kable(...) as the final expression (without print()) to
# render the table, or use results='asis' if you explicitly print it.
raw_kable <- knitr::kable(
    raw_display,
    format = "pipe",
    row.names = FALSE,
    caption = "Raw contrast effects and confidence intervals (happiness points)"
)
cat("\n\n--- Combined table as Markdown ---\n\n")
print(raw_kable)

# 6. A styled flextable, matching the vignette's presentation approach -------
# Keep every displayed column. Set alignment and headings without changing
# the estimates, intervals, or correction method. Notes printed by psyci()
# are separate from its table, so add the relevant interpretation as a footer.
style_display_table <- function(display, caption) {
    table <- flextable::flextable(display)
    table <- flextable::set_caption(table, caption)
    table <- flextable::set_header_labels(
        table, family = "Family", contrast = "Contrast",
        estimate = "Estimate", SE = "SE", df1 = "df1", df2 = "df2",
        F_value = "F", F_crit = "Critical F", p_adj = "Adjusted p",
        method = "Procedure", sig = "Sig.", lower = "Lower CI", upper = "Upper CI"
    )
    table <- flextable::theme_booktabs(table)
    table <- flextable::fontsize(table, size = 9, part = "all")
    table <- flextable::align(table, align = "center", part = "all")
    table <- flextable::align(
        table, j = c("family", "contrast"), align = "left", part = "all"
    )
    table <- flextable::add_footer_lines(
        table,
        "95% confidence intervals; Bonferroni adjustment within each contrast family."
    )
    flextable::autofit(table)
}

raw_flextable <- style_display_table(
    raw_display, "Raw contrast effects and confidence intervals (happiness points)"
)
standardized_flextable <- style_display_table(
    standardized_display,
    "Standardized contrast effects and confidence intervals (sample SD units)"
)

# In a vignette, finish one chunk with raw_flextable, and another chunk with
# standardized_flextable. Each will render as a formatted table. To show only
# the raw table, include only raw_flextable. Add echo=FALSE to hide its code.
# Detailed per-family notes remain in attr(output[[i]], "mesg").

# 7. Preview both formatted tables without knitting a whole vignette --------
# The preview goes to R's temporary directory, not the repository.
table_preview_path <- tempfile("psyci-vignette-tables-", fileext = ".html")
flextable::save_as_html(
    "Raw effects" = raw_flextable,
    "Standardized effects" = standardized_flextable,
    path = table_preview_path,
    title = "PsyR vignette table examples"
)
cat("\nHTML preview saved to:\n", table_preview_path, "\n", sep = "")
# Run this line to open it in your browser:
# utils::browseURL(table_preview_path)
