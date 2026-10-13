# Internal helpers for detecting factor participation.

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

