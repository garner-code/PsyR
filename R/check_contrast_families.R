#' Infer or check families for emmeans contrast tables
#'
#' Infers whether each supplied [emmeans::contrast()] table belongs to a
#' between-subject (`"b"`), within-subject (`"w"`), or between-by-within
#' (`"bw"`) family. When `family_list` is supplied, the inferred result is used
#' to check the declaration. The check uses factor roles stored in an
#' `afex_aov` model and the linear functions retained by the original
#' `emmGrid` object. When `emmeans` has marginalized or conditioned away a
#' model factor, the checker recovers coefficients on the complete model-factor
#' grid before classifying the table.
#'
#' This function checks the subject-domain classification of each table. It
#' deliberately does not assess whether contrast coefficients are zero-sum,
#' orthogonal, product/rank-one, planned independently of the data, or members
#' of a coherent multiplicity family.
#'
#' For custom contrasts on a full cell-means grid, same-direction weights are
#' treated as marginal averaging only when their normalized joint profile is
#' common across the other subject-domain strata. A one-level selector is
#' treated as conditioning. A partial-support averaging pattern with more than
#' one but fewer than all levels is reported as unverifiable when its role
#' could change the family classification.
#'
#' @param model An ANOVA model inheriting from `afex_aov`.
#' @param contrast_tables An `emmGrid`, `summary_emm`, or list of such contrast
#'   tables produced by [emmeans::contrast()]. Complete checking requires the
#'   original `emmGrid`; a `summary_emm` has lost its coefficient grid.
#' @param family_list Optional list or character vector containing one family
#'   code per contrast table. Valid codes are `"b"`, `"w"`, and `"bw"`. A
#'   single code is recycled over multiple tables. When `NULL` (the default),
#'   the family is inferred from the model and the retained contrast
#'   coefficients. An explicit value is checked against the inferred family.
#' @param between_factors Optional between-subject factor names to cross-check
#'   against the model metadata and the factors detected in each table.
#' @param within_factors Optional within-subject factor names to cross-check
#'   against the model metadata and the factors detected in each table.
#' @param mode How mismatches or unverifiable tables are handled: `"error"`,
#'   `"warn"`, or `"none"`. With `"none"`, classifications are returned but
#'   no family conditions are signalled. Strict checking requires the original
#'   `emmGrid`; a `summary_emm` is necessarily unverifiable because its
#'   coefficient grid is no longer available.
#' @param tol Numeric tolerance used only when detecting whether coefficient
#'   patterns vary over a factor. It is not used to judge contrast validity.
#'
#' @return An object of class `psyr_family_check`. It is a list containing an
#'   overall validity flag, the resolved `families`, whether they were
#'   `inferred`, model factor roles, normalized fallback factor arguments,
#'   detailed table classifications, a table-level summary, and any issues.
#' @export
#'
#' @examples
#' data(fatigue)
#' fatigue$Group <- factor(fatigue$Group)
#' model <- afex::aov_ez("id", "Errors", fatigue, between = "Group")
#' marginal_means <- emmeans::emmeans(model, "Group")
#' contrasts <- emmeans::contrast(
#'     marginal_means,
#'     list("Alert versus others" = c(-1, 1 / 3, 1 / 3, 1 / 3))
#' )
#' check_contrast_families(
#'     model,
#'     contrasts
#' )
check_contrast_families <- function(
    model,
    contrast_tables,
    family_list = NULL,
    between_factors = NA,
    within_factors = NA,
    mode = c("error", "warn", "none"),
    tol = sqrt(.Machine$double.eps)
) {
    mode <- match.arg(mode)
    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol < 0) {
        stop("`tol` must be one finite, non-negative number.", call. = FALSE)
    }

    if (inherits(contrast_tables, "emmGrid") ||
        inherits(contrast_tables, "summary_emm")) {
        contrast_tables <- list(contrast_tables)
    }
    if (!is.list(contrast_tables) || length(contrast_tables) == 0L) {
        stop(
            "`contrast_tables` must be an emmGrid, summary_emm, or non-empty list.",
            call. = FALSE
        )
    }

    valid_table <- vapply(
        contrast_tables,
        function(x) inherits(x, "emmGrid") || inherits(x, "summary_emm"),
        logical(1)
    )
    if (!all(valid_table)) {
        stop(
            "Every element of `contrast_tables` must be an emmGrid or summary_emm.",
            call. = FALSE
        )
    }

    infer_families <- is.null(family_list)
    if (infer_families) {
        families <- rep(NA_character_, length(contrast_tables))
    } else if (is.list(family_list)) {
        families <- unlist(family_list, use.names = FALSE)
    } else if (is.character(family_list)) {
        families <- family_list
    } else {
        stop(
            "`family_list` must be NULL, a list, or a character vector.",
            call. = FALSE
        )
    }
    families <- as.character(families)
    if (length(families) == 1L && length(contrast_tables) > 1L) {
        families <- rep(families, length(contrast_tables))
    }
    if (length(families) != length(contrast_tables)) {
        stop(
            "The length of `family_list` must match `contrast_tables`.",
            call. = FALSE
        )
    }
    invalid_families <- setdiff(
        unique(families[!is.na(families)]),
        c("b", "w", "bw")
    )
    if (!infer_families &&
        (length(invalid_families) > 0L || anyNA(families))) {
        stop(
            "Every family code must be one of 'b', 'w', or 'bw'.",
            call. = FALSE
        )
    }

    normalize_factors <- function(x) {
        x <- unlist(x, use.names = FALSE)
        if (length(x) == 0L || all(is.na(x))) {
            return(character())
        }
        unique(as.character(x[!is.na(x)]))
    }

    roles <- .get_afex_factor_roles(model)
    if (mode == "none") {
        factor_arguments <- list(
            between = normalize_factors(between_factors),
            within = normalize_factors(within_factors)
        )
    } else {
        if (!roles$verifiable) {
            stop(
                paste(
                    "PsyR could not establish between/within factor roles:",
                    paste(roles$issues, collapse = " ")
                ),
                call. = FALSE
            )
        }
        factor_arguments <- .validate_factor_role_arguments(
            roles,
            between_factors = between_factors,
            within_factors = within_factors
        )
    }

    table_names <- names(contrast_tables)
    if (is.null(table_names)) {
        table_names <- rep("", length(contrast_tables))
    }
    table_ids <- vapply(seq_along(contrast_tables), function(index) {
        if (nzchar(table_names[[index]])) {
            sprintf("contrast_tables[['%s']]", table_names[[index]])
        } else {
            sprintf("contrast_tables[[%d]]", index)
        }
    }, character(1))

    table_bases <- vapply(
        contrast_tables,
        .get_emm_basis_type,
        character(1)
    )
    required_bases <- unique(table_bases[!is.na(table_bases)])
    reference_grids <- stats::setNames(
        lapply(required_bases, function(basis) {
            .build_family_reference_grid(model, roles, basis = basis)
        }),
        required_bases
    )
    classifications <- lapply(seq_along(contrast_tables), function(index) {
        basis <- table_bases[[index]]
        reference_grid <- if (!is.na(basis) && basis %in% names(reference_grids)) {
            reference_grids[[basis]]
        } else {
            NULL
        }
        .classify_contrast_table(
            contrast_table = contrast_tables[[index]],
            roles = roles,
            reference_grid = reference_grid,
            tol = tol
        )
    })
    detected_families <- vapply(
        classifications,
        function(x) x$detected_family,
        character(1)
    )
    resolved_families <- if (infer_families) {
        ifelse(
            detected_families %in% c("b", "w", "bw"),
            detected_families,
            NA_character_
        )
    } else {
        families
    }

    validations <- vector("list", length(classifications))
    for (index in seq_along(classifications)) {
        validations[[index]] <- if (infer_families) {
            .validate_inferred_family(
                classification = classifications[[index]],
                table_id = table_ids[[index]],
                mode = mode
            )
        } else {
            .validate_family_assignment(
                classification = classifications[[index]],
                declared_family = families[[index]],
                table_id = table_ids[[index]],
                factor_arguments = factor_arguments,
                mode = mode
            )
        }
        if (!is.null(validations[[index]]$message)) {
            classifications[[index]]$issues <- unique(c(
                classifications[[index]]$issues,
                validations[[index]]$message
            ))
        }
    }

    status <- vapply(validations, function(x) x$status, character(1))
    summary_table <- data.frame(
        table = table_ids,
        declared_family = if (infer_families) {
            rep(NA_character_, length(families))
        } else {
            families
        },
        detected_family = detected_families,
        resolved_family = resolved_families,
        scope = vapply(
            classifications,
            function(x) x$scope,
            character(1)
        ),
        status = status,
        stringsAsFactors = FALSE
    )
    issues <- c(
        unlist(lapply(reference_grids, function(x) x$issues), use.names = FALSE),
        unlist(
            lapply(classifications, function(x) x$issues),
            use.names = FALSE
        ),
        unlist(
            lapply(validations, function(x) x$message),
            use.names = FALSE
        )
    )

    structure(
        list(
            valid = all(status == "pass"),
            roles = roles,
            factor_arguments = factor_arguments,
            families = resolved_families,
            inferred = infer_families,
            tables = classifications,
            summary = summary_table,
            issues = unique(issues)
        ),
        class = c("psyr_family_check", "list")
    )
}
