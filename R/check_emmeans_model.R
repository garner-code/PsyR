#' Require multivariate emmeans tables for within and mixed families
#'
#' Checks the model basis retained by original `emmGrid` contrast objects.
#' Within-subject (`"w"`) and between-by-within (`"bw"`) families require
#' the multivariate basis. Both declared and detected families are checked;
#' a conflicting declaration cannot bypass this requirement.
#'
#' @param model An ANOVA model inheriting from `afex_aov`.
#' @param contrast_tables An original `emmGrid` contrast table or a non-empty
#'   list of original tables from [emmeans::contrast()].
#' @param family_list Optional list or character vector of family codes
#'   (`"b"`, `"w"`, `"bw"`). A single code is recycled over multiple tables.
#'   When `NULL`, families are inferred using [check_contrast_families()].
#' @return Invisibly, a data frame with `table`, `family`, `detected_basis`,
#'   `multivariate_required`, and `valid` columns. Stops if a required basis
#'   is univariate or unknown, or if an inferred family cannot be resolved.
#'   Between-subject tables do not require a particular basis. This function
#'   checks the basis requirement, not agreement between family declarations
#'   and detected families; use [check_contrast_families()] for that check.
#' @details The current global `afex_options("emmeans_model")` setting cannot
#'   establish which basis an existing table uses. Changing that option does
#'   not update existing tables. Recreate the marginal means with
#'   `emmeans::emmeans(model, ..., model = "multivariate")`, then recreate
#'   contrasts with the original definitions, averaging weights, conditioning,
#'   and other analysis settings. The ANOVA normally does not need refitting.
#'   This check does not change options, refit models, or repair tables.
#'   [psyci()] enforces the same requirement even with `family_check = "warn"`
#'   or `family_check = "none"`.
#' @export
#' @examples
#' data(spacing)
#' model <- afex::aov_ez(
#'     "subj", "yield", spacing, within = "spacing", between = "group"
#' )
#' marginal_means <- emmeans::emmeans(
#'     model, "spacing", model = "multivariate"
#' )
#' contrasts <- emmeans::contrast(
#'     marginal_means, list("Twenty - Forty" = c(1, -1, 0))
#' )
#' diagnostics <- check_emmeans_model(model, contrasts)
check_emmeans_model <- function(model, contrast_tables, family_list = NULL) {
    if (inherits(contrast_tables, "emmGrid")) {
        contrast_tables <- list(contrast_tables)
    }
    if (!is.list(contrast_tables) || length(contrast_tables) == 0L ||
        !all(vapply(contrast_tables, inherits, logical(1), "emmGrid"))) {
        stop(
            "`contrast_tables` must be an original emmGrid or a non-empty list of original emmGrid objects.",
            call. = FALSE
        )
    }
    report <- check_contrast_families(
        model, contrast_tables, family_list = family_list, mode = "none"
    )
    .check_emmeans_model(contrast_tables, report)
}

# Share the resolved family report with psyci() rather than classifying twice.
.check_emmeans_model <- function(contrast_tables, family_report) {
    required <- vapply(seq_along(contrast_tables), function(index) {
        classification <- family_report$tables[[index]]
        families <- c(
            family_report$families[[index]],
            classification$detected_family,
            classification$contrasts$detected_family
        )
        any(families %in% c("w", "bw"))
    }, logical(1))
    basis <- vapply(contrast_tables, .get_emm_basis_type, character(1))
    valid <- !required | (!is.na(basis) & basis == "multivariate")
    diagnostics <- data.frame(
        table = family_report$summary$table,
        family = family_report$families,
        detected_basis = unname(basis),
        multivariate_required = required,
        valid = unname(valid),
        stringsAsFactors = FALSE
    )
    if (any(!valid)) {
        affected <- paste0(
            diagnostics$table[!valid], " (basis: ",
            ifelse(is.na(basis[!valid]), "unknown", basis[!valid]), ")"
        )
        stop(
            paste0(
                "PsyR requires multivariate emmeans contrast tables for within-subject ",
                "and between-by-within families. Affected tables: ",
                paste(affected, collapse = "; "), ". ",
                "Recreate the marginal means with emmeans::emmeans(model, ..., ",
                "model = \"multivariate\"), then recreate the contrasts using the ",
                "original definitions and analysis settings. Changing afex_options() ",
                "alone does not update existing tables; the ANOVA normally does not need refitting."
            ),
            call. = FALSE
        )
    }
    if (anyNA(family_report$families)) {
        stop(
            paste0(
                "PsyR could not resolve the family for ",
                paste(diagnostics$table[is.na(family_report$families)], collapse = ", "),
                ". Supply original, verifiable emmGrid contrasts or an explicit family_list."
            ),
            call. = FALSE
        )
    }
    invisible(diagnostics)
}
