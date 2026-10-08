# Internal helpers for checking declared contrast families.

.get_afex_factor_roles <- function(model) {
    if (!inherits(model, "afex_aov")) {
        stop("`model` must inherit from class 'afex_aov'.", call. = FALSE)
    }

    role_names <- function(x) {
        if (is.null(x)) {
            return(character())
        }

        object_names <- names(x)
        if (!is.null(object_names) && all(nzchar(object_names))) {
            return(object_names)
        }

        if (is.character(x)) {
            return(as.character(x))
        }

        character()
    }

    between_metadata <- attr(model, "between", exact = TRUE)
    within_metadata <- attr(model, "within", exact = TRUE)
    between <- role_names(between_metadata)
    within <- role_names(within_metadata)
    sources <- character()

    if (!is.null(between_metadata)) {
        sources <- c(sources, "afex_between_attribute")
    }
    if (!is.null(within_metadata)) {
        sources <- c(sources, "afex_within_attribute")
    }

    if (is.null(within_metadata) && !is.null(model$data$idata)) {
        within <- names(model$data$idata)
        sources <- c(sources, "afex_idata")
    }

    if (is.null(between_metadata) && !is.null(model$lm)) {
        term_labels <- attr(stats::terms(model$lm), "term.labels")
        long_names <- names(model$data$long)
        between <- intersect(term_labels, long_names)
        sources <- c(sources, "afex_lm_terms")
    }

    between <- unique(between[nzchar(between)])
    within <- unique(within[nzchar(within)])
    overlap <- intersect(between, within)
    issues <- character()

    if (length(overlap) > 0L) {
        issues <- c(
            issues,
            sprintf(
                "Factors assigned to both between- and within-subject roles: %s.",
                paste(overlap, collapse = ", ")
            )
        )
    }

    long_data <- model$data$long
    model_factors <- unique(c(between, within))
    missing_from_data <- setdiff(model_factors, names(long_data))
    if (length(missing_from_data) > 0L) {
        issues <- c(
            issues,
            sprintf(
                "Model factors missing from the afex long data: %s.",
                paste(missing_from_data, collapse = ", ")
            )
        )
    }

    factor_levels <- stats::setNames(vector("list", length(model_factors)), model_factors)
    for (factor_name in model_factors) {
        if (!is.null(long_data) && factor_name %in% names(long_data)) {
            values <- long_data[[factor_name]]
            factor_levels[[factor_name]] <- if (is.factor(values)) {
                levels(values)
            } else {
                unique(as.character(values))
            }
        } else if (factor_name %in% names(between_metadata)) {
            factor_levels[[factor_name]] <- as.character(between_metadata[[factor_name]])
        } else if (factor_name %in% names(within_metadata)) {
            factor_levels[[factor_name]] <- as.character(within_metadata[[factor_name]])
        }
    }

    list(
        between = between,
        within = within,
        levels = factor_levels,
        source = unique(sources),
        verifiable = length(model_factors) > 0L && length(issues) == 0L,
        issues = issues
    )
}

.validate_factor_role_arguments <- function(
    roles,
    between_factors = NA,
    within_factors = NA
) {
    normalize_factors <- function(x) {
        if (is.null(x)) {
            return(character())
        }

        x <- unlist(x, use.names = FALSE)
        if (length(x) == 0L || all(is.na(x))) {
            return(character())
        }

        unique(as.character(x[!is.na(x)]))
    }

    supplied_between <- normalize_factors(between_factors)
    supplied_within <- normalize_factors(within_factors)
    supplied <- unique(c(supplied_between, supplied_within))
    known <- unique(c(roles$between, roles$within))

    unknown <- setdiff(supplied, known)
    if (length(unknown) > 0L) {
        stop(
            sprintf(
                "Factor names not found in the afex model metadata: %s.",
                paste(unknown, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    between_as_within <- intersect(supplied_between, roles$within)
    if (length(between_as_within) > 0L) {
        stop(
            sprintf(
                "Factors supplied as between-subject are within-subject in the model: %s.",
                paste(between_as_within, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    within_as_between <- intersect(supplied_within, roles$between)
    if (length(within_as_between) > 0L) {
        stop(
            sprintf(
                "Factors supplied as within-subject are between-subject in the model: %s.",
                paste(within_as_between, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    overlap <- intersect(supplied_between, supplied_within)
    if (length(overlap) > 0L) {
        stop(
            sprintf(
                "Factors supplied as both between- and within-subject: %s.",
                paste(overlap, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    list(
        between = supplied_between,
        within = supplied_within
    )
}

.get_emm_basis_type <- function(emm_grid) {
    if (!inherits(emm_grid, "emmGrid")) {
        return(NA_character_)
    }

    model_info <- tryCatch(
        methods::slot(emm_grid, "model.info"),
        error = function(error) NULL
    )
    model_call <- model_info$call
    if (is.null(model_call) || length(model_call) == 0L) {
        return(NA_character_)
    }

    call_head <- paste(deparse(model_call[[1L]]), collapse = "")
    if (grepl("(^|::)aov$", call_head)) {
        return("univariate")
    }
    if (grepl("(^|::)lm$", call_head)) {
        return("multivariate")
    }
    NA_character_
}

.build_family_reference_grid <- function(
    model,
    roles,
    basis = c("multivariate", "univariate")
) {
    basis <- match.arg(basis)
    factor_names <- unique(c(roles$between, roles$within))
    if (length(factor_names) == 0L) {
        return(list(
            available = FALSE,
            grid = data.frame(),
            linfct = matrix(numeric(), nrow = 0L, ncol = 0L),
            issues = "No model factors were available for a full reference grid."
        ))
    }

    reference <- tryCatch(
        suppressMessages(emmeans::emmeans(
            model,
            specs = factor_names,
            model = basis
        )),
        error = function(error) error
    )
    if (inherits(reference, "error")) {
        return(list(
            available = FALSE,
            grid = data.frame(),
            linfct = matrix(numeric(), nrow = 0L, ncol = 0L),
            issues = sprintf(
                "Could not construct the full emmeans reference grid: %s",
                conditionMessage(reference)
            )
        ))
    }

    actual_basis <- .get_emm_basis_type(reference)
    if (!identical(actual_basis, basis)) {
        return(list(
            available = FALSE,
            grid = data.frame(),
            linfct = matrix(numeric(), nrow = 0L, ncol = 0L),
            basis = actual_basis,
            issues = sprintf(
                paste0(
                    "Could not construct a full emmeans reference grid in ",
                    "the required %s basis."
                ),
                basis
            )
        ))
    }

    reference_grid <- methods::slot(reference, "grid")
    missing_factors <- setdiff(factor_names, names(reference_grid))
    if (length(missing_factors) > 0L) {
        return(list(
            available = FALSE,
            grid = data.frame(),
            linfct = matrix(numeric(), nrow = 0L, ncol = 0L),
            issues = sprintf(
                "Full emmeans reference grid omitted model factors: %s.",
                paste(missing_factors, collapse = ", ")
            )
        ))
    }

    reference_grid <- reference_grid[, factor_names, drop = FALSE]
    expected_cells <- prod(vapply(
        roles$levels[factor_names],
        length,
        integer(1)
    ))
    level_complete <- all(vapply(factor_names, function(factor_name) {
        setequal(
            unique(as.character(reference_grid[[factor_name]])),
            roles$levels[[factor_name]]
        )
    }, logical(1)))
    if (!level_complete || nrow(unique(reference_grid)) != expected_cells) {
        return(list(
            available = FALSE,
            grid = reference_grid,
            linfct = matrix(numeric(), nrow = 0L, ncol = 0L),
            issues = paste(
                "The emmeans reference grid was not a complete Cartesian grid",
                "over the afex model factors."
            )
        ))
    }

    list(
        available = TRUE,
        grid = reference_grid,
        linfct = methods::slot(reference, "linfct"),
        basis = basis,
        issues = character()
    )
}

.reconstruct_family_coefficients <- function(
    contrast_table,
    reference_grid,
    tol = sqrt(.Machine$double.eps)
) {
    if (is.null(reference_grid) || !isTRUE(reference_grid$available)) {
        return(list(success = FALSE, coefficients = NULL, issue = NULL))
    }

    contrast_linfct <- tryCatch(
        methods::slot(contrast_table, "linfct"),
        error = function(error) NULL
    )
    reference_linfct <- reference_grid$linfct
    if (is.null(contrast_linfct) ||
        ncol(contrast_linfct) != ncol(reference_linfct)) {
        return(list(
            success = FALSE,
            coefficients = NULL,
            issue = paste(
                "The contrast and full reference grids used incompatible",
                "linear-function bases."
            )
        ))
    }

    coefficients <- tryCatch(
        qr.solve(
            t(reference_linfct),
            t(contrast_linfct),
            tol = sqrt(.Machine$double.eps)
        ),
        error = function(error) NULL
    )
    if (is.null(coefficients)) {
        return(list(
            success = FALSE,
            coefficients = NULL,
            issue = "The full-grid contrast coefficients could not be recovered."
        ))
    }
    if (is.null(dim(coefficients))) {
        coefficients <- matrix(coefficients, ncol = nrow(contrast_linfct))
    }

    reconstructed <- t(reference_linfct) %*% coefficients
    target <- t(contrast_linfct)
    target_scale <- max(abs(target), na.rm = TRUE)
    if (!is.finite(target_scale)) {
        target_scale <- 1
    }
    residual_tolerance <- sqrt(.Machine$double.eps) *
        max(target_scale, .Machine$double.xmin)
    residual <- max(abs(reconstructed - target), na.rm = TRUE)
    if (!is.finite(residual) || residual > residual_tolerance) {
        return(list(
            success = FALSE,
            coefficients = NULL,
            issue = paste(
                "Recovered full-grid coefficients did not reproduce the",
                "emmeans linear functions within numerical tolerance."
            )
        ))
    }

    for (index in seq_len(ncol(coefficients))) {
        coefficient_scale <- max(abs(coefficients[, index]), na.rm = TRUE)
        threshold <- tol * coefficient_scale
        coefficients[abs(coefficients[, index]) <= threshold, index] <- 0
    }
    colnames(coefficients) <- paste0("c.", seq_len(ncol(coefficients)))

    list(success = TRUE, coefficients = coefficients, issue = NULL)
}

.extract_emm_family_data <- function(
    contrast_table,
    reference_grid = NULL,
    tol = sqrt(.Machine$double.eps)
) {
    if (inherits(contrast_table, "summary_emm")) {
        summary_table <- as.data.frame(contrast_table)
        contrast_labels <- if ("contrast" %in% names(summary_table)) {
            as.character(summary_table$contrast)
        } else {
            rownames(summary_table)
        }

        return(list(
            grid = data.frame(),
            coefficients = matrix(numeric(), nrow = 0L, ncol = 0L),
            contrasts = data.frame(
                index = seq_along(contrast_labels),
                label = contrast_labels,
                stringsAsFactors = FALSE
            ),
            by_factors = character(),
            averaged_factors = character(),
            level_order = list(),
            provenance = "unavailable",
            issues = paste(
                "A summary_emm object no longer contains the coefficient grid",
                "needed to verify its family. Pass the original emmGrid."
            )
        ))
    }

    if (!inherits(contrast_table, "emmGrid")) {
        stop(
            "Each contrast table must inherit from 'emmGrid' or 'summary_emm'.",
            call. = FALSE
        )
    }

    coefficient_data <- tryCatch(
        stats::coef(contrast_table),
        error = function(error) NULL
    )

    if (is.null(coefficient_data) || !is.data.frame(coefficient_data)) {
        return(list(
            grid = data.frame(),
            coefficients = matrix(numeric(), nrow = 0L, ncol = 0L),
            contrasts = data.frame(),
            by_factors = character(),
            averaged_factors = character(),
            level_order = list(),
            provenance = "unavailable",
            issues = "emmeans did not provide a coefficient grid for this table."
        ))
    }

    coefficient_names <- grep(
        "^c\\.[0-9]+$",
        names(coefficient_data),
        value = TRUE
    )
    if (length(coefficient_names) == 0L) {
        return(list(
            grid = coefficient_data,
            coefficients = matrix(numeric(), nrow = nrow(coefficient_data), ncol = 0L),
            contrasts = data.frame(),
            by_factors = character(),
            averaged_factors = character(),
            level_order = lapply(coefficient_data, function(x) unique(as.character(x))),
            provenance = "unavailable",
            issues = "No emmeans coefficient columns were found."
        ))
    }

    suffix <- as.integer(sub("^c\\.", "", coefficient_names))
    coefficient_names <- coefficient_names[order(suffix)]
    grid_names <- setdiff(names(coefficient_data), coefficient_names)
    grid <- coefficient_data[, grid_names, drop = FALSE]
    coefficients <- as.matrix(
        coefficient_data[, coefficient_names, drop = FALSE]
    )
    storage.mode(coefficients) <- "double"

    summary_table <- as.data.frame(summary(contrast_table))
    misc <- tryCatch(
        methods::slot(contrast_table, "misc"),
        error = function(error) list()
    )
    by_factors <- if (!is.null(misc$by.vars)) {
        as.character(misc$by.vars)
    } else {
        character()
    }
    averaged_factors <- if (!is.null(misc$avgd.over)) {
        as.character(misc$avgd.over)
    } else {
        character()
    }

    factor_names <- if (!is.null(reference_grid) &&
        isTRUE(reference_grid$available)) {
        names(reference_grid$grid)
    } else {
        character()
    }
    original_grid_complete <- length(factor_names) > 0L &&
        all(factor_names %in% names(grid)) &&
        all(vapply(factor_names, function(factor_name) {
            setequal(
                unique(as.character(grid[[factor_name]])),
                unique(as.character(reference_grid$grid[[factor_name]]))
            )
        }, logical(1))) &&
        nrow(unique(grid[, factor_names, drop = FALSE])) ==
            nrow(reference_grid$grid)

    reconstruction_issue <- NULL
    reconstructed <- FALSE
    reconstruction_failed <- FALSE
    if (length(factor_names) > 0L && !original_grid_complete) {
        reconstruction <- .reconstruct_family_coefficients(
            contrast_table,
            reference_grid = reference_grid,
            tol = tol
        )
        if (isTRUE(reconstruction$success)) {
            grid <- reference_grid$grid
            coefficients <- reconstruction$coefficients
            averaged_factors <- character()
            reconstructed <- TRUE
        } else {
            reconstruction_issue <- reconstruction$issue
            reconstruction_failed <- TRUE
        }
    }

    contrast_labels <- if ("contrast" %in% names(summary_table)) {
        as.character(summary_table$contrast)
    } else {
        custom_columns <- grep("_custom$", names(summary_table), value = TRUE)
        if (length(custom_columns) > 0L) {
            apply(
                summary_table[, custom_columns, drop = FALSE],
                1L,
                function(values) paste(values, collapse = " x ")
            )
        } else {
            rownames(summary_table)
        }
    }
    if (is.null(contrast_labels)) {
        contrast_labels <- paste("Contrast", seq_len(nrow(summary_table)))
    }

    present_by <- intersect(by_factors, names(summary_table))
    if (length(present_by) > 0L && nrow(summary_table) > 0L) {
        by_labels <- apply(
            summary_table[, present_by, drop = FALSE],
            1L,
            function(values) {
                paste(paste0(present_by, "=", values), collapse = ", ")
            }
        )
        contrast_labels <- paste0(contrast_labels, " [", by_labels, "]")
    }

    issues <- character()
    if (!is.null(reconstruction_issue)) {
        issues <- c(issues, reconstruction_issue)
    }
    provenance <- if (reconstruction_failed) {
        "unavailable"
    } else if (reconstructed) {
        "reconstructed_full_grid"
    } else {
        "complete"
    }
    if (ncol(coefficients) != nrow(summary_table)) {
        if (provenance != "unavailable") {
            provenance <- "partial"
        }
        issues <- c(
            issues,
            sprintf(
                paste(
                    "The coefficient grid contains %d contrasts but the summary",
                    "contains %d rows."
                ),
                ncol(coefficients),
                nrow(summary_table)
            )
        )
    }

    n_contrasts <- ncol(coefficients)
    if (length(contrast_labels) < n_contrasts) {
        contrast_labels <- c(
            contrast_labels,
            paste("Contrast", seq.int(length(contrast_labels) + 1L, n_contrasts))
        )
    }

    list(
        grid = grid,
        coefficients = coefficients,
        contrasts = data.frame(
            index = seq_len(n_contrasts),
            label = contrast_labels[seq_len(n_contrasts)],
            stringsAsFactors = FALSE
        ),
        by_factors = by_factors,
        averaged_factors = averaged_factors,
        level_order = lapply(grid, function(x) unique(as.character(x))),
        provenance = provenance,
        issues = issues
    )
}

.detect_factor_participation <- function(
    coefficients,
    grid,
    roles,
    by_factors = character(),
    averaged_factors = character(),
    tol = sqrt(.Machine$double.eps),
    .collapse_domains = TRUE
) {
    if (!is.matrix(coefficients)) {
        coefficients <- as.matrix(coefficients)
    }

    role_names <- unique(c(roles$between, roles$within))
    role_lookup <- c(
        stats::setNames(rep("between", length(roles$between)), roles$between),
        stats::setNames(rep("within", length(roles$within)), roles$within)
    )
    unknown_grid_factors <- setdiff(names(grid), role_names)
    unknown_grid_factors <- unknown_grid_factors[
        !grepl("^\\.", unknown_grid_factors)
    ]
    role_levels_complete <- length(role_names) > 0L &&
        all(role_names %in% names(grid)) &&
        all(vapply(role_names, function(factor_name) {
            setequal(
                unique(as.character(grid[[factor_name]])),
                roles$levels[[factor_name]]
            )
        }, logical(1)))
    joint_grid_incomplete <- length(role_names) > 1L &&
        role_levels_complete &&
        nrow(unique(grid[, role_names, drop = FALSE])) <
            prod(vapply(roles$levels[role_names], length, integer(1)))

    classify_one <- function(coefficient) {
        finite_coefficients <- coefficient[is.finite(coefficient)]
        coefficient_scale <- if (length(finite_coefficients) > 0L) {
            max(abs(finite_coefficients))
        } else {
            0
        }
        threshold <- tol * coefficient_scale
        factor_status <- stats::setNames(rep("absent", length(role_names)), role_names)
        comparable <- stats::setNames(rep(FALSE, length(role_names)), role_names)
        strong_evidence <- stats::setNames(rep(FALSE, length(role_names)), role_names)
        profile_mismatch_evidence <- stats::setNames(
            rep(FALSE, length(role_names)),
            role_names
        )
        incomplete_levels <- list()

        for (factor_name in role_names) {
            if (!factor_name %in% names(grid)) {
                if (factor_name %in% averaged_factors) {
                    factor_status[[factor_name]] <- "averaged"
                }
                next
            }

            factor_values <- grid[[factor_name]]
            observed_levels <- unique(as.character(factor_values))
            expected_levels <- roles$levels[[factor_name]]
            omitted_levels <- setdiff(expected_levels, observed_levels)
            if (length(omitted_levels) > 0L) {
                factor_status[[factor_name]] <- "unknown"
                incomplete_levels[[factor_name]] <- omitted_levels
                next
            }

            if (factor_name %in% averaged_factors) {
                factor_status[[factor_name]] <- "averaged"
                next
            }

            if (length(observed_levels) <= 1L) {
                if (length(expected_levels) > 1L && any(abs(coefficient) > threshold)) {
                    factor_status[[factor_name]] <- "conditioning"
                    strong_evidence[[factor_name]] <- TRUE
                } else {
                    factor_status[[factor_name]] <- "inactive"
                }
                next
            }

            other_factors <- intersect(
                setdiff(role_names, factor_name),
                names(grid)
            )
            if (length(other_factors) == 0L) {
                groups <- list(seq_along(coefficient))
            } else {
                group_key <- interaction(
                    grid[, other_factors, drop = FALSE],
                    drop = TRUE,
                    lex.order = TRUE
                )
                groups <- split(seq_along(coefficient), group_key)
            }

            varies <- FALSE
            averaging_pattern <- TRUE
            conditioning_pattern <- FALSE
            sign_change_pattern <- FALSE
            profile_mismatch <- FALSE
            reference_profile <- NULL
            selector_pattern <- TRUE
            saw_nonzero <- FALSE
            for (indices in groups) {
                group_levels <- unique(as.character(factor_values[indices]))
                if (length(group_levels) <= 1L) {
                    next
                }

                comparable[[factor_name]] <- TRUE
                group_coefficients <- coefficient[indices]
                nonzero <- abs(group_coefficients) > threshold
                nonzero_coefficients <- group_coefficients[nonzero]
                if (length(nonzero_coefficients) > 0L) {
                    saw_nonzero <- TRUE
                    has_positive <- any(nonzero_coefficients > 0)
                    has_negative <- any(nonzero_coefficients < 0)
                    if (has_positive && has_negative) {
                        averaging_pattern <- FALSE
                        sign_change_pattern <- TRUE
                    }
                    if (sum(nonzero) < length(group_coefficients)) {
                        averaging_pattern <- FALSE
                        conditioning_pattern <- TRUE
                    }
                    if (sum(nonzero) != 1L) {
                        selector_pattern <- FALSE
                    }

                    profile <- vapply(
                        observed_levels,
                        function(level) {
                            sum(abs(group_coefficients[
                                as.character(factor_values[indices]) == level
                            ]))
                        },
                        numeric(1)
                    )
                    profile <- profile / sum(profile)
                    if (is.null(reference_profile)) {
                        reference_profile <- profile
                    } else if (any(abs(profile - reference_profile) > tol)) {
                        profile_mismatch <- TRUE
                        averaging_pattern <- FALSE
                    }
                }
                if (
                    any(!is.finite(group_coefficients)) ||
                    diff(range(group_coefficients, na.rm = TRUE)) > threshold
                ) {
                    varies <- TRUE
                }
            }

            if (varies) {
                factor_status[[factor_name]] <- if (factor_name %in% by_factors) {
                    "conditioning"
                } else if (sign_change_pattern || profile_mismatch) {
                    "varying"
                } else if (conditioning_pattern && selector_pattern) {
                    "conditioning"
                } else if (conditioning_pattern) {
                    "unknown"
                } else if (averaging_pattern && saw_nonzero) {
                    "averaged"
                } else {
                    "varying"
                }
            } else if (factor_name %in% by_factors && any(abs(coefficient) > threshold)) {
                factor_status[[factor_name]] <- "conditioning"
            } else if (comparable[[factor_name]]) {
                factor_status[[factor_name]] <- "inactive"
            } else {
                factor_status[[factor_name]] <- "unknown"
            }
            strong_evidence[[factor_name]] <-
                sign_change_pattern ||
                (conditioning_pattern && selector_pattern) ||
                (factor_name %in% by_factors &&
                    any(abs(coefficient) > threshold))
            profile_mismatch_evidence[[factor_name]] <- profile_mismatch
        }

        missing_context <- names(factor_status)[factor_status == "absent"]
        ambiguous_context <- names(factor_status)[factor_status == "unknown"]
        metadata_complete <- length(missing_context) == 0L &&
            length(ambiguous_context) == 0L &&
            length(unknown_grid_factors) == 0L &&
            !joint_grid_incomplete

        list(
            factors = data.frame(
                factor = role_names,
                role = unname(role_lookup[role_names]),
                status = unname(factor_status[role_names]),
                participates = unname(
                    factor_status[role_names] %in% c("varying", "conditioning")
                ),
                .strong = unname(strong_evidence[role_names]),
                .profile_mismatch = unname(
                    profile_mismatch_evidence[role_names]
                ),
                stringsAsFactors = FALSE
            ),
            missing_context = missing_context,
            ambiguous_context = ambiguous_context,
            incomplete_levels = incomplete_levels,
            unknown_grid_factors = unknown_grid_factors,
            joint_grid_incomplete = joint_grid_incomplete,
            metadata_complete = metadata_complete
        )
    }

    results <- lapply(seq_len(ncol(coefficients)), function(index) {
        classify_one(coefficients[, index])
    })

    combine_factor_values <- function(factor_names) {
        if (length(factor_names) == 0L) {
            return(rep(".all", nrow(grid)))
        }
        if (length(factor_names) == 1L) {
            return(as.character(grid[[factor_names]]))
        }
        as.character(interaction(
            grid[, factor_names, drop = FALSE],
            drop = TRUE,
            lex.order = TRUE,
            sep = "\r"
        ))
    }

    can_refine_joint_profiles <- .collapse_domains &&
        role_levels_complete &&
        !joint_grid_incomplete
    if (can_refine_joint_profiles) {
        for (index in seq_along(results)) {
            for (role in c("between", "within")) {
                role_rows <- which(results[[index]]$factors$role == role)
                if (length(role_rows) < 2L) {
                    next
                }

                role_report <- results[[index]]$factors[role_rows, , drop = FALSE]
                weak_factors <- role_report$factor[
                    role_report$status == "varying" &
                        role_report$.profile_mismatch &
                        !role_report$.strong
                ]
                strong_factors <- role_report$factor[role_report$.strong]
                unresolved <- role_report$status %in% c("unknown", "absent")
                if (length(weak_factors) == 0L ||
                    length(strong_factors) == 0L || any(unresolved)) {
                    next
                }

                nuisance_factors <- role_report$factor[
                    !role_report$.strong &
                        role_report$status %in% c(
                            "varying",
                            "averaged",
                            "inactive"
                        )
                ]
                profile_nuisance_factors <- unique(c(
                    nuisance_factors,
                    results[[index]]$factors$factor[
                        results[[index]]$factors$status %in% c(
                            "averaged",
                            "inactive"
                        )
                    ]
                ))
                active_strata <- results[[index]]$factors$factor[
                    !results[[index]]$factors$factor %in%
                        profile_nuisance_factors
                ]
                profile_grid <- data.frame(
                    .active_stratum = combine_factor_values(active_strata),
                    .nuisance_profile = combine_factor_values(
                        profile_nuisance_factors
                    ),
                    stringsAsFactors = FALSE
                )
                profile_roles <- list(
                    between = ".active_stratum",
                    within = ".nuisance_profile",
                    levels = list(
                        .active_stratum = unique(profile_grid$.active_stratum),
                        .nuisance_profile = unique(profile_grid$.nuisance_profile)
                    )
                )
                profile_result <- .detect_factor_participation(
                    coefficients = coefficients[, index, drop = FALSE],
                    grid = profile_grid,
                    roles = profile_roles,
                    by_factors = character(),
                    averaged_factors = character(),
                    tol = tol,
                    .collapse_domains = FALSE
                )[[1L]]
                nuisance_row <- profile_result$factors[
                    profile_result$factors$factor == ".nuisance_profile",
                    ,
                    drop = FALSE
                ]
                active_row <- profile_result$factors[
                    profile_result$factors$factor == ".active_stratum",
                    ,
                    drop = FALSE
                ]
                if (nrow(nuisance_row) == 1L &&
                    identical(nuisance_row$status[[1L]], "averaged") &&
                    nrow(active_row) == 1L &&
                    isTRUE(active_row$participates[[1L]])) {
                    nuisance_rows <- results[[index]]$factors$factor %in%
                        nuisance_factors
                    results[[index]]$factors$status[nuisance_rows] <- "averaged"
                    results[[index]]$factors$participates[nuisance_rows] <- FALSE
                    results[[index]]$ambiguous_context <- setdiff(
                        results[[index]]$ambiguous_context,
                        nuisance_factors
                    )
                }
            }
        }
    }

    can_collapse_domains <- .collapse_domains &&
        length(roles$between) > 0L &&
        length(roles$within) > 0L &&
        all(role_names %in% names(grid)) &&
        all(vapply(role_names, function(factor_name) {
            setequal(
                unique(as.character(grid[[factor_name]])),
                roles$levels[[factor_name]]
            )
        }, logical(1))) &&
        nrow(unique(grid[, role_names, drop = FALSE])) ==
            prod(vapply(roles$levels[role_names], length, integer(1)))

    if (can_collapse_domains) {
        domain_grid <- data.frame(
            .between_domain = combine_factor_values(roles$between),
            .within_domain = combine_factor_values(roles$within),
            stringsAsFactors = FALSE
        )
        domain_roles <- list(
            between = ".between_domain",
            within = ".within_domain",
            levels = list(
                .between_domain = unique(domain_grid$.between_domain),
                .within_domain = unique(domain_grid$.within_domain)
            )
        )
        domain_by <- c(
            if (any(roles$between %in% by_factors)) ".between_domain",
            if (any(roles$within %in% by_factors)) ".within_domain"
        )
        domain_results <- .detect_factor_participation(
            coefficients = coefficients,
            grid = domain_grid,
            roles = domain_roles,
            by_factors = domain_by,
            averaged_factors = character(),
            tol = tol,
            .collapse_domains = FALSE
        )

        for (index in seq_along(results)) {
            results[[index]]$domains <- domain_results[[index]]$factors
            for (role in c("between", "within")) {
                domain_name <- paste0(".", role, "_domain")
                domain_row <- domain_results[[index]]$factors[
                    domain_results[[index]]$factors$factor == domain_name,
                    ,
                    drop = FALSE
                ]
                role_factors <- roles[[role]]
                if (nrow(domain_row) == 0L) {
                    next
                }

                if (!domain_row$participates[[1L]] &&
                    domain_row$status[[1L]] %in% c("averaged", "inactive")) {
                    factor_rows <- results[[index]]$factors$factor %in% role_factors
                    results[[index]]$factors$status[factor_rows] <-
                        domain_row$status[[1L]]
                    results[[index]]$factors$participates[factor_rows] <- FALSE
                    results[[index]]$ambiguous_context <- setdiff(
                        results[[index]]$ambiguous_context,
                        role_factors
                    )
                } else if (domain_row$participates[[1L]] &&
                    !any(results[[index]]$factors$participates[
                        results[[index]]$factors$role == role
                    ])) {
                    factor_rows <- results[[index]]$factors$factor %in% role_factors
                    results[[index]]$factors$status[factor_rows] <- "varying"
                    results[[index]]$factors$participates[factor_rows] <- TRUE
                }
            }
            results[[index]]$metadata_complete <-
                length(results[[index]]$missing_context) == 0L &&
                length(results[[index]]$ambiguous_context) == 0L &&
                length(results[[index]]$unknown_grid_factors) == 0L &&
                !results[[index]]$joint_grid_incomplete &&
                domain_results[[index]]$metadata_complete
        }
    } else {
        for (index in seq_along(results)) {
            results[[index]]$domains <- data.frame()
        }
    }

    results
}

.classify_contrast_family <- function(participation, roles) {
    factor_report <- participation$factors
    active <- factor_report[factor_report$participates, , drop = FALSE]
    has_between <- any(active$role == "between")
    has_within <- any(active$role == "within")
    unresolved_between <- any(
        factor_report$role == "between" &
            factor_report$status %in% c("unknown", "absent")
    )
    unresolved_within <- any(
        factor_report$role == "within" &
            factor_report$status %in% c("unknown", "absent")
    )
    has_unknown_grid_factors <-
        length(participation$unknown_grid_factors) > 0L

    detected_family <- if (has_unknown_grid_factors ||
        isTRUE(participation$joint_grid_incomplete)) {
        "unknown"
    } else if (has_between && has_within) {
        "bw"
    } else if (has_between && !unresolved_within) {
        "b"
    } else if (has_within && !unresolved_between) {
        "w"
    } else {
        "unknown"
    }

    scope <- if (detected_family == "unknown") {
        "unknown"
    } else if (unresolved_between || unresolved_within) {
        "unknown"
    } else if (any(active$status == "conditioning")) {
        "conditioned"
    } else if (any(factor_report$status %in% c("averaged", "inactive"))) {
        "marginal"
    } else {
        "joint"
    }

    list(
        family = detected_family,
        scope = scope,
        between_factors = active$factor[active$role == "between"],
        within_factors = active$factor[active$role == "within"],
        evidence = stats::setNames(factor_report$status, factor_report$factor)
    )
}

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
