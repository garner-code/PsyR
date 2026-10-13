# Internal helpers for recovering contrast coefficients.

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

