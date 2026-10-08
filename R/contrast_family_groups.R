.validate_family_group <- function(family_group, n_tables, method) {
    if (is.null(family_group)) {
        return(NULL)
    }
    if (method != "bf") {
        stop(
            "`family_group` is currently supported only when `method = \"bf\"`.",
            call. = FALSE
        )
    }
    if (length(family_group) != n_tables) {
        stop(
            "`family_group` must contain one identifier per contrast table.",
            call. = FALSE
        )
    }
    if (is.list(family_group)) {
        family_group <- unlist(family_group, use.names = FALSE)
    }
    if (!is.atomic(family_group) || anyNA(family_group)) {
        stop("`family_group` identifiers must be non-missing atomic values.", call. = FALSE)
    }
    family_group <- as.character(family_group)
    if (any(!nzchar(trimws(family_group)))) {
        stop("`family_group` identifiers cannot be empty.", call. = FALSE)
    }
    family_group
}

.contrast_grid_key <- function(grid, factor_names) {
    if (!all(factor_names %in% names(grid))) {
        return(NULL)
    }
    do.call(
        paste,
        c(lapply(grid[factor_names], as.character), sep = "\r")
    )
}

.align_contrast_coefficients <- function(extraction, roles) {
    factor_names <- unique(c(roles$between, roles$within))
    grid_key <- .contrast_grid_key(extraction$grid, factor_names)
    if (is.null(grid_key) || anyDuplicated(grid_key)) {
        return(NULL)
    }
    level_grid <- expand.grid(
        roles$levels[factor_names],
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    canonical_key <- .contrast_grid_key(level_grid, factor_names)
    order <- match(canonical_key, grid_key)
    if (anyNA(order)) {
        return(NULL)
    }
    extraction$coefficients[order, , drop = FALSE]
}

.assess_contrast_table_orthogonality <- function(
    extractions,
    roles,
    table_ids,
    tol
) {
    aligned <- lapply(extractions, .align_contrast_coefficients, roles = roles)
    table_reports <- lapply(seq_along(extractions), function(index) {
        list(
            table = table_ids[[index]],
            orthogonal_to = character(),
            not_orthogonal_to = character(),
            unavailable_for = character()
        )
    })
    pairs <- list()
    if (length(extractions) < 2L) {
        return(list(tables = table_reports, pairs = pairs))
    }

    for (left in seq_len(length(extractions) - 1L)) {
        for (right in seq.int(left + 1L, length(extractions))) {
            if (is.null(aligned[[left]]) || is.null(aligned[[right]])) {
                status <- "unavailable"
                nonorthogonal_pairs <- data.frame()
            } else {
                products <- crossprod(aligned[[left]], aligned[[right]])
                left_norm <- sqrt(colSums(aligned[[left]]^2))
                right_norm <- sqrt(colSums(aligned[[right]]^2))
                scale <- outer(left_norm, right_norm)
                relative <- abs(products) / scale
                relative[!is.finite(relative)] <- Inf
                locations <- which(relative > tol, arr.ind = TRUE)
                status <- if (nrow(locations) > 0L) {
                    "not_orthogonal"
                } else {
                    "orthogonal"
                }
                nonorthogonal_pairs <- if (nrow(locations) == 0L) {
                    data.frame()
                } else {
                    data.frame(
                        left_index = locations[, 1L],
                        left_label = extractions[[left]]$contrasts$label[locations[, 1L]],
                        right_index = locations[, 2L],
                        right_label = extractions[[right]]$contrasts$label[locations[, 2L]],
                        inner_product = products[locations],
                        relative_inner_product = relative[locations],
                        stringsAsFactors = FALSE
                    )
                }
            }
            pairs[[length(pairs) + 1L]] <- list(
                left = left,
                right = right,
                left_table = table_ids[[left]],
                right_table = table_ids[[right]],
                status = status,
                nonorthogonal_contrasts = nonorthogonal_pairs
            )
            field <- switch(
                status,
                orthogonal = "orthogonal_to",
                not_orthogonal = "not_orthogonal_to",
                unavailable = "unavailable_for"
            )
            table_reports[[left]][[field]] <- c(
                table_reports[[left]][[field]], table_ids[[right]]
            )
            table_reports[[right]][[field]] <- c(
                table_reports[[right]][[field]], table_ids[[left]]
            )
        }
    }
    list(tables = table_reports, pairs = pairs)
}

.prepare_psyci_contrast_tables <- function(
    contrast_tables,
    family_report,
    tol = sqrt(.Machine$double.eps)
) {
    if (!all(vapply(contrast_tables, inherits, logical(1), "emmGrid"))) {
        stop(
            paste0(
                "PsyR requires the original emmGrid objects to identify genuine ",
                "contrasts and assess orthogonality. A summary_emm no longer ",
                "contains the required coefficient information."
            ),
            call. = FALSE
        )
    }
    table_ids <- family_report$summary$table
    prepared <- vector("list", length(contrast_tables))
    dropped <- vector("list", length(contrast_tables))
    retained_extractions <- vector("list", length(contrast_tables))

    for (index in seq_along(contrast_tables)) {
        extraction <- family_report$tables[[index]]$extraction
        if (extraction$provenance == "unavailable" ||
            ncol(extraction$coefficients) == 0L) {
            stop(
                sprintf(
                    "PsyR could not recover coefficients for %s. Pass its original emmGrid.",
                    table_ids[[index]]
                ),
                call. = FALSE
            )
        }
        coefficient_scale <- colSums(abs(extraction$coefficients))
        coefficient_sum <- colSums(extraction$coefficients)
        genuine <- coefficient_scale > 0 &
            abs(coefficient_sum) <= tol * coefficient_scale
        if (!any(genuine)) {
            stop(
                sprintf(
                    "%s contains no genuine contrasts after non-contrast rows are removed.",
                    table_ids[[index]]
                ),
                call. = FALSE
            )
        }
        removed <- which(!genuine)
        dropped[[index]] <- data.frame(
            index = removed,
            label = extraction$contrasts$label[removed],
            coefficient_sum = coefficient_sum[removed],
            stringsAsFactors = FALSE
        )
        prepared[[index]] <- if (all(genuine)) {
            contrast_tables[[index]]
        } else {
            contrast_tables[[index]][genuine]
        }
        retained <- extraction
        retained$coefficients <- extraction$coefficients[, genuine, drop = FALSE]
        retained$contrasts <- extraction$contrasts[genuine, , drop = FALSE]
        retained_extractions[[index]] <- retained
    }

    list(
        tables = prepared,
        dropped = dropped,
        extractions = retained_extractions,
        counts = vapply(retained_extractions, function(x) ncol(x$coefficients), integer(1)),
        table_ids = table_ids,
        orthogonality = .assess_contrast_table_orthogonality(
            retained_extractions,
            roles = family_report$roles,
            table_ids = table_ids,
            tol = tol
        )
    )
}

.bonferroni_group_sizes <- function(counts, groups, alphas) {
    alpha_values <- unlist(alphas, use.names = FALSE)
    for (group in unique(groups)) {
        members <- which(groups == group)
        if (length(unique(alpha_values[members])) != 1L) {
            stop(
                sprintf(
                    "All contrast tables in family_group '%s' must use the same alpha.",
                    group
                ),
                call. = FALSE
            )
        }
    }
    totals <- vapply(unique(groups), function(group) {
        sum(counts[groups == group])
    }, integer(1))
    names(totals) <- unique(groups)
    unname(totals[groups])
}

.warn_separate_nonorthogonal_tables <- function(
    orthogonality,
    correction_groups,
    method,
    independent
) {
    separately_corrected <- if (method == "bf") {
        function(left, right) correction_groups[[left]] != correction_groups[[right]]
    } else if (method == "ph" && independent) {
        function(left, right) TRUE
    } else {
        function(left, right) FALSE
    }
    problematic <- Filter(function(pair) {
        pair$status == "not_orthogonal" &&
            separately_corrected(pair$left, pair$right)
    }, orthogonality$pairs)
    if (length(problematic) == 0L) {
        return(invisible(NULL))
    }
    labels <- vapply(problematic, function(pair) {
        paste(pair$left_table, "and", pair$right_table)
    }, character(1))
    warning(
        paste0(
            "Non-orthogonality was detected between separately corrected ",
            "contrast tables: ", paste(labels, collapse = "; "),
            ". If these tables form one planned inferential family, assign ",
            "them the same `family_group`."
        ),
        call. = FALSE
    )
    invisible(NULL)
}

.attach_psyci_contrast_diagnostics <- function(
    result,
    index,
    preparation,
    correction_group,
    family_size
) {
    dropped <- preparation$dropped[[index]]
    orthogonality <- preparation$orthogonality$tables[[index]]
    messages <- character()
    if (nrow(dropped) > 0L) {
        messages <- c(messages, sprintf(
            "PsyR dropped %d non-contrast row(s) before correction: %s.",
            nrow(dropped), paste(dropped$label, collapse = ", ")
        ))
    } else {
        messages <- c(messages, "PsyR found no non-contrast rows to drop.")
    }
    if (!is.na(correction_group) && !is.na(family_size)) {
        messages <- c(messages, sprintf(
            paste0(
                "PsyR Bonferroni confidence intervals and `psyr_p_value` use ",
                "family_group '%s', containing %d genuine contrast(s); ",
                "`p.value` retains the original emmeans result."
            ),
            correction_group, family_size
        ))
    }
    if (length(preparation$orthogonality$tables) == 1L) {
        messages <- c(messages,
            "Orthogonality with other tables was not applicable because only one table was supplied.")
    } else {
        if (length(orthogonality$not_orthogonal_to) > 0L) {
            messages <- c(messages, paste0(
                "Orthogonality: this table is not orthogonal to ",
                paste(orthogonality$not_orthogonal_to, collapse = ", "), "."
            ))
        }
        if (length(orthogonality$orthogonal_to) > 0L) {
            messages <- c(messages, paste0(
                "Orthogonality: this table is orthogonal to ",
                paste(orthogonality$orthogonal_to, collapse = ", "), "."
            ))
        }
        if (length(orthogonality$unavailable_for) > 0L) {
            messages <- c(messages, paste0(
                "Orthogonality could not be assessed against ",
                paste(orthogonality$unavailable_for, collapse = ", "), "."
            ))
        }
    }
    attr(result, "psyr_dropped_noncontrasts") <- dropped
    attr(result, "psyr_orthogonality") <- orthogonality
    attr(result, "psyr_family_group") <- correction_group
    attr(result, "psyr_family_size") <- family_size
    attr(result, "mesg") <- c(attr(result, "mesg"), messages)
    result
}
