# Internal helpers for classifying and validating contrast families.

.classify_contrast_table <- function(
    contrast_table,
    roles,
    reference_grid = NULL,
    tol = sqrt(.Machine$double.eps)
) {
    extracted <- .extract_emm_family_data(
        contrast_table,
        reference_grid = reference_grid,
        tol = tol
    )
    if (extracted$provenance == "unavailable") {
        contrast_report <- extracted$contrasts
        if (nrow(contrast_report) > 0L) {
            contrast_report$detected_family <- "unknown"
            contrast_report$scope <- "unknown"
            contrast_report$between_factors <- ""
            contrast_report$within_factors <- ""
        }

        return(list(
            detected_family = "unknown",
            scope = "unknown",
            between_factors = character(),
            within_factors = character(),
            contrasts = contrast_report,
            extraction = extracted,
            participation = list(),
            issues = extracted$issues
        ))
    }

    participation <- .detect_factor_participation(
        coefficients = extracted$coefficients,
        grid = extracted$grid,
        roles = roles,
        by_factors = extracted$by_factors,
        averaged_factors = extracted$averaged_factors,
        tol = tol
    )
    classifications <- lapply(
        participation,
        .classify_contrast_family,
        roles = roles
    )

    contrast_report <- data.frame(
        index = extracted$contrasts$index,
        label = extracted$contrasts$label,
        detected_family = vapply(
            classifications,
            function(x) x$family,
            character(1)
        ),
        scope = vapply(
            classifications,
            function(x) x$scope,
            character(1)
        ),
        between_factors = vapply(
            classifications,
            function(x) paste(x$between_factors, collapse = ", "),
            character(1)
        ),
        within_factors = vapply(
            classifications,
            function(x) paste(x$within_factors, collapse = ", "),
            character(1)
        ),
        stringsAsFactors = FALSE
    )

    detected <- unique(contrast_report$detected_family)
    known_detected <- setdiff(detected, "unknown")
    table_family <- if (length(known_detected) > 1L) {
        "mixed"
    } else if ("unknown" %in% detected) {
        "unknown"
    } else if (length(known_detected) == 1L) {
        known_detected
    } else {
        "unknown"
    }

    scopes <- unique(contrast_report$scope)
    table_scope <- if (length(scopes) == 1L) scopes else "mixed"
    issues <- extracted$issues
    if (table_family == "mixed") {
        issues <- c(
            issues,
            "The table contains contrasts from more than one detected family."
        )
    }
    if (table_family == "unknown") {
        issues <- c(
            issues,
            "The family could not be verified from the available emmeans metadata."
        )
    }
    missing_context <- unique(unlist(
        lapply(participation, function(x) x$missing_context),
        use.names = FALSE
    ))
    if (length(missing_context) > 0L) {
        issues <- c(
            issues,
            sprintf(
                "Model factors absent from the emmeans coefficient grid: %s.",
                paste(missing_context, collapse = ", ")
            )
        )
    }
    ambiguous_context <- unique(unlist(
        lapply(participation, function(x) x$ambiguous_context),
        use.names = FALSE
    ))
    if (length(ambiguous_context) > 0L) {
        issues <- c(
            issues,
            sprintf(
                paste0(
                    "Participation could not be resolved for model factors: ",
                    "%s. This can occur with incomplete grids or custom ",
                    "weights containing structural zeros."
                ),
                paste(ambiguous_context, collapse = ", ")
            )
        )
    }
    incomplete_level_details <- unique(unlist(lapply(
        participation,
        function(x) {
            if (length(x$incomplete_levels) == 0L) {
                return(character())
            }
            vapply(names(x$incomplete_levels), function(factor_name) {
                sprintf(
                    "%s missing {%s}",
                    factor_name,
                    paste(x$incomplete_levels[[factor_name]], collapse = ", ")
                )
            }, character(1))
        }
    ), use.names = FALSE))
    if (length(incomplete_level_details) > 0L) {
        issues <- c(
            issues,
            sprintf(
                "The coefficient grid omits model levels: %s.",
                paste(incomplete_level_details, collapse = "; ")
            )
        )
    }
    if (any(vapply(
        participation,
        function(x) isTRUE(x$joint_grid_incomplete),
        logical(1)
    ))) {
        issues <- c(
            issues,
            paste(
                "The coefficient grid omits combinations of model-factor",
                "levels, so family participation is ambiguous."
            )
        )
    }
    unknown_grid_factors <- unique(unlist(
        lapply(participation, function(x) x$unknown_grid_factors),
        use.names = FALSE
    ))
    if (length(unknown_grid_factors) > 0L) {
        issues <- c(
            issues,
            sprintf(
                "Coefficient-grid factors absent from the afex role metadata: %s.",
                paste(unknown_grid_factors, collapse = ", ")
            )
        )
    }

    list(
        detected_family = table_family,
        scope = table_scope,
        between_factors = unique(unlist(
            lapply(classifications, function(x) x$between_factors),
            use.names = FALSE
        )),
        within_factors = unique(unlist(
            lapply(classifications, function(x) x$within_factors),
            use.names = FALSE
        )),
        contrasts = contrast_report,
        extraction = extracted,
        participation = participation,
        issues = unique(issues)
    )
}

.validate_family_assignment <- function(
    classification,
    declared_family,
    table_id,
    factor_arguments = list(between = character(), within = character()),
    mode = c("error", "warn", "none")
) {
    mode <- match.arg(mode)
    declared_family <- as.character(declared_family)[1L]
    allowed <- c("b", "w", "bw")
    if (length(declared_family) != 1L || is.na(declared_family) ||
        !declared_family %in% allowed) {
        stop(
            sprintf(
                "Declared family for %s must be one of 'b', 'w', or 'bw'.",
                table_id
            ),
            call. = FALSE
        )
    }

    detected <- classification$detected_family
    known_families <- if (nrow(classification$contrasts) > 0L) {
        setdiff(unique(classification$contrasts$detected_family), "unknown")
    } else {
        character()
    }
    definite_mismatch <- length(known_families) > 0L &&
        !declared_family %in% known_families

    factor_mismatch <- FALSE
    factor_mismatch_details <- character()
    missing_between <- setdiff(
        classification$between_factors,
        factor_arguments$between
    )
    if (length(factor_arguments$between) > 0L &&
        length(missing_between) > 0L) {
        factor_mismatch <- TRUE
        factor_mismatch_details <- c(
            factor_mismatch_details,
            sprintf(
                paste0(
                    "between_factors supplied {%s}, detected {%s}; ",
                    "missing {%s}"
                ),
                paste(factor_arguments$between, collapse = ", "),
                paste(classification$between_factors, collapse = ", "),
                paste(missing_between, collapse = ", ")
            )
        )
    }
    missing_within <- setdiff(
        classification$within_factors,
        factor_arguments$within
    )
    if (length(factor_arguments$within) > 0L &&
        length(missing_within) > 0L) {
        factor_mismatch <- TRUE
        factor_mismatch_details <- c(
            factor_mismatch_details,
            sprintf(
                paste0(
                    "within_factors supplied {%s}, detected {%s}; ",
                    "missing {%s}"
                ),
                paste(factor_arguments$within, collapse = ", "),
                paste(classification$within_factors, collapse = ", "),
                paste(missing_within, collapse = ", ")
            )
        )
    }

    status <- if (factor_mismatch) {
        "factor_mismatch"
    } else if (identical(detected, declared_family)) {
        "pass"
    } else if (identical(detected, "mixed")) {
        "mixed"
    } else if (identical(detected, "unknown") && !definite_mismatch) {
        "unverifiable"
    } else {
        "mismatch"
    }

    if (status == "pass") {
        return(list(status = status, message = NULL))
    }

    contrast_details <- classification$contrasts
    detail_text <- if (nrow(contrast_details) > 0L) {
        paste(
            sprintf(
                "'%s'=%s",
                contrast_details$label,
                contrast_details$detected_family
            ),
            collapse = "; "
        )
    } else {
        "no per-contrast classification was available"
    }

    message <- switch(
        status,
        mismatch = sprintf(
            paste0(
                "PsyR family check failed for %s. Declared family: '%s'; ",
                "table-level detected family: '%s'; known per-contrast ",
                "families: {%s}. Contrasts: %s."
            ),
            table_id,
            declared_family,
            detected,
            if (length(known_families) > 0L) {
                paste(known_families, collapse = ", ")
            } else {
                "none"
            },
            detail_text
        ),
        mixed = sprintf(
            paste0(
                "PsyR family check found mixed families in %s, which has one ",
                "family_list entry ('%s'). Split the contrasts into separate ",
                "emmeans tables. Contrasts: %s."
            ),
            table_id,
            declared_family,
            detail_text
        ),
        unverifiable = sprintf(
            paste0(
                "PsyR could not verify family '%s' for %s. The available ",
                "emmeans object does not retain complete coefficient-grid and ",
                "model-factor context. %s Use mode = 'warn' or mode = 'none' ",
                "only when this limitation is acceptable."
            ),
            declared_family,
            table_id,
            if (length(classification$issues) > 0L) {
                paste0("Details: ", paste(classification$issues, collapse = " "))
            } else {
                ""
            }
        ),
        factor_mismatch = sprintf(
            paste0(
                "PsyR factor check failed for %s: %s. The factors used to ",
                "describe the family must include the factors detected in ",
                "the contrast table."
            ),
            table_id,
            paste(factor_mismatch_details, collapse = "; ")
        )
    )

    if (mode == "none") {
        return(list(status = status, message = message))
    }

    condition_class <- if (mode == "error") {
        c("psyr_family_error", "error", "condition")
    } else {
        c("psyr_family_warning", "warning", "condition")
    }
    condition <- structure(
        list(
            message = message,
            call = NULL,
            table_id = table_id,
            declared_family = declared_family,
            detected_family = detected,
            classification = classification
        ),
        class = condition_class
    )

    if (mode == "error") {
        stop(condition)
    }
    warning(condition)
    list(status = status, message = message)
}

.validate_inferred_family <- function(
    classification,
    table_id,
    mode = c("error", "warn", "none")
) {
    mode <- match.arg(mode)
    detected <- classification$detected_family
    if (detected %in% c("b", "w", "bw")) {
        return(list(status = "pass", message = NULL))
    }

    contrast_details <- classification$contrasts
    detail_text <- if (nrow(contrast_details) > 0L) {
        paste(
            sprintf(
                "'%s'=%s",
                contrast_details$label,
                contrast_details$detected_family
            ),
            collapse = "; "
        )
    } else {
        "no per-contrast classification was available"
    }
    status <- if (identical(detected, "mixed")) "mixed" else "unverifiable"
    message <- if (status == "mixed") {
        sprintf(
            paste0(
                "PsyR inferred more than one family in %s. Split the ",
                "contrasts into separate emmeans tables before calling ",
                "psyci(). Contrasts: %s."
            ),
            table_id,
            detail_text
        )
    } else {
        sprintf(
            paste0(
                "PsyR could not infer a family for %s from the available ",
                "emmeans coefficient-grid and model-factor context. %s ",
                "Supply the original emmGrid, or provide family_list ",
                "explicitly when the family is known."
            ),
            table_id,
            if (length(classification$issues) > 0L) {
                paste0("Details: ", paste(classification$issues, collapse = " "))
            } else {
                ""
            }
        )
    }

    if (mode == "none") {
        return(list(status = status, message = message))
    }
    condition_class <- if (mode == "error") {
        c("psyr_family_error", "error", "condition")
    } else {
        c("psyr_family_warning", "warning", "condition")
    }
    condition <- structure(
        list(
            message = message,
            call = NULL,
            table_id = table_id,
            declared_family = NULL,
            detected_family = detected,
            classification = classification
        ),
        class = condition_class
    )
    if (mode == "error") {
        stop(condition)
    }
    warning(condition)
    list(status = status, message = message)
}

