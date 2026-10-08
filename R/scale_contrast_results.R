.psyci_sample_sd <- function(model, tol = sqrt(.Machine$double.eps)) {
    long_data <- model$data$long
    response <- attr(model, "dv", exact = TRUE)
    subject <- attr(model, "id", exact = TRUE)
    between_factors <- names(attr(model, "between", exact = TRUE))
    within_factors <- names(attr(model, "within", exact = TRUE))

    required_columns <- unique(c(
        response,
        subject,
        between_factors,
        within_factors
    ))
    missing_columns <- setdiff(required_columns, names(long_data))
    if (!is.data.frame(long_data) || length(response) != 1L ||
        !nzchar(response) || length(subject) != 1L || !nzchar(subject) ||
        length(missing_columns) > 0L) {
        return(.unavailable_psyci_scale(
            "the afex model does not retain the response and design data needed for scaling"
        ))
    }

    values <- long_data[[response]]
    if (!is.numeric(values) || any(!is.finite(values))) {
        return(.unavailable_psyci_scale(
            "the model response contains missing or non-finite values"
        ))
    }

    subject_cells <- .psyci_cell_factor(long_data, subject)
    repeated_cells <- .psyci_cell_factor(long_data, within_factors)
    between_cells <- .psyci_cell_factor(long_data, between_factors)
    observation_cells <- interaction(
        subject_cells,
        repeated_cells,
        drop = TRUE,
        lex.order = TRUE
    )
    if (anyDuplicated(observation_cells)) {
        return(.unavailable_psyci_scale(
            "the afex model data contain more than one response per subject and repeated-measure cell"
        ))
    }

    repeated_rows <- split(seq_len(nrow(long_data)), repeated_cells)
    cell_variances <- vapply(repeated_rows, function(rows) {
        group <- droplevels(between_cells[rows])
        cell_means <- ave(values[rows], group, FUN = mean)
        residual_df <- length(rows) - nlevels(group)
        if (residual_df <= 0L) {
            return(NA_real_)
        }
        sum((values[rows] - cell_means)^2) / residual_df
    }, numeric(1))

    if (any(!is.finite(cell_variances)) || any(cell_variances < 0)) {
        return(.unavailable_psyci_scale(
            "the pooled within-cell variance has no residual degrees of freedom"
        ))
    }

    divisor <- sqrt(mean(cell_variances))
    if (!is.finite(divisor) || divisor <= tol) {
        return(.unavailable_psyci_scale(
            "the pooled within-cell sample SD is zero"
        ))
    }

    structure(
        list(
            available = TRUE,
            divisor = divisor,
            cell_sd = sqrt(cell_variances),
            response = response,
            between_factors = between_factors,
            within_factors = within_factors,
            reason = NULL
        ),
        class = "psyr_scale_info"
    )
}

.psyci_cell_factor <- function(data, columns) {
    if (length(columns) == 0L) {
        return(factor(rep("1", nrow(data))))
    }
    interaction(data[columns], drop = TRUE, lex.order = TRUE)
}

.unavailable_psyci_scale <- function(reason) {
    structure(
        list(
            available = FALSE,
            divisor = NA_real_,
            cell_sd = numeric(),
            response = NULL,
            between_factors = character(),
            within_factors = character(),
            reason = reason
        ),
        class = "psyr_scale_info"
    )
}

.scale_psyci_table <- function(table, divisor) {
    scaled <- table
    scaled_columns <- intersect(
        c("estimate", "SE", "lower", "upper"),
        names(scaled)
    )
    scaled[scaled_columns] <- lapply(
        scaled[scaled_columns],
        function(column) column / divisor
    )
    attr(scaled, "psyr_scale_divisor") <- divisor
    attr(scaled, "psyr_scale_units") <- "sample SD"
    attr(scaled, "mesg") <- c(
        attr(scaled, "mesg"),
        paste(
            "PsyR scaled estimate, SE, and confidence limits by sample SD:",
            format(divisor, digits = 7)
        )
    )
    scaled
}

.psyci_ak_scale_data <- function(model) {
    long_data <- model$data$long
    response <- attr(model, "dv", exact = TRUE)
    subject <- attr(model, "id", exact = TRUE)
    between_factors <- names(attr(model, "between", exact = TRUE))
    within_factors <- names(attr(model, "within", exact = TRUE))
    if (length(within_factors) == 0L) {
        stop(
            "Algina-Keselman scaling requires a within-subject factor.",
            call. = FALSE
        )
    }

    subject_cells <- .psyci_cell_factor(long_data, subject)
    repeated_cells <- .psyci_cell_factor(long_data, within_factors)
    between_cells <- .psyci_cell_factor(long_data, between_factors)
    observation_cells <- interaction(
        subject_cells,
        repeated_cells,
        drop = TRUE,
        lex.order = TRUE
    )
    if (anyDuplicated(observation_cells)) {
        stop(
            paste(
                "Algina-Keselman scaling requires one response per subject",
                "and repeated-measure cell."
            ),
            call. = FALSE
        )
    }

    repeated_rows <- split(seq_len(nrow(long_data)), repeated_cells)
    cell_variance <- vapply(repeated_rows, function(rows) {
        group <- droplevels(between_cells[rows])
        cell_means <- ave(long_data[[response]][rows], group, FUN = mean)
        residual_df <- length(rows) - nlevels(group)
        if (residual_df <= 0L) {
            return(NA_real_)
        }
        sum((long_data[[response]][rows] - cell_means)^2) / residual_df
    }, numeric(1))
    if (any(!is.finite(cell_variance)) || any(cell_variance <= 0)) {
        stop(
            paste(
                "Algina-Keselman scaling requires positive finite variance",
                "in every repeated-measure cell."
            ),
            call. = FALSE
        )
    }
    cell_grid <- do.call(rbind, lapply(repeated_rows, function(rows) {
        long_data[rows[[1L]], within_factors, drop = FALSE]
    }))
    rownames(cell_grid) <- NULL

    list(
        within_factors = within_factors,
        grid = cell_grid,
        variance = unname(cell_variance)
    )
}

.psyci_grid_key <- function(grid, factors) {
    do.call(
        paste,
        c(lapply(grid[factors], as.character), list(sep = "\r"))
    )
}

.psyci_ak_contrast_scales <- function(
    extraction,
    scale_data,
    tol = sqrt(.Machine$double.eps)
) {
    factors <- scale_data$within_factors
    if (!all(factors %in% names(extraction$grid))) {
        stop(
            paste(
                "PsyR could not map the contrast coefficient grid to the",
                "repeated-measure cells needed for Algina-Keselman scaling."
            ),
            call. = FALSE
        )
    }
    extraction_key <- .psyci_grid_key(extraction$grid, factors)
    variance_key <- .psyci_grid_key(scale_data$grid, factors)
    if (anyDuplicated(variance_key)) {
        stop(
            "Repeated-measure cells were not unique during Algina-Keselman scaling.",
            call. = FALSE
        )
    }

    vapply(seq_len(ncol(extraction$coefficients)), function(index) {
        coefficients <- extraction$coefficients[, index]
        coefficient_scale <- max(abs(coefficients))
        participating_rows <- abs(coefficients) > tol * coefficient_scale
        participating_keys <- unique(extraction_key[participating_rows])
        matched <- match(participating_keys, variance_key)
        if (length(matched) != 2L || anyNA(matched)) {
            stop(
                paste(
                    "Each Algina-Keselman contrast must involve exactly two",
                    "identifiable repeated-measure cells, as specified by",
                    "the published pairwise method."
                ),
                call. = FALSE
            )
        }
        sqrt(mean(scale_data$variance[matched]))
    }, numeric(1))
}

.scale_psyci_table_ak <- function(table, divisor, alpha, adjust) {
    interval <- algina_keselman_ci(
        estimate = table$estimate,
        se_contrast = table$SE,
        df = table$df,
        scale = divisor,
        alpha = alpha,
        adjust = adjust
    )
    scaled <- table
    scaled$estimate <- interval$effect_size
    scaled$SE <- interval$standard_error
    scaled$lower <- interval$lower
    scaled$upper <- interval$upper
    scaled$ncp <- interval$ncp
    scaled$ncp_lower <- interval$ncp_lower
    scaled$ncp_upper <- interval$ncp_upper
    scaled$scale <- interval$scale
    attr(scaled, "psyr_scale_divisor") <- divisor
    attr(scaled, "psyr_scale_units") <- "sample SD"
    attr(scaled, "psyr_standardized_ci_method") <- "Algina-Keselman"
    attr(scaled, "mesg") <- c(
        attr(scaled, "mesg"),
        paste0(
            "PsyR used Algina-Keselman noncentral-t confidence limits with ",
            "contrast-specific sample-SD scaling and Bonferroni adjust = ",
            adjust,
            "."
        )
    )
    scaled
}

.new_psyci_result <- function(
    tables,
    model,
    standardized_ci = "bird",
    extractions = NULL,
    alphas = NULL,
    adjustments = NULL
) {
    scale_info <- .psyci_sample_sd(model)
    scaled_tables <- if (!isTRUE(scale_info$available)) {
        NULL
    } else if (identical(standardized_ci, "algina_keselman")) {
        scale_data <- .psyci_ak_scale_data(model)
        divisors <- lapply(
            extractions,
            .psyci_ak_contrast_scales,
            scale_data = scale_data
        )
        Map(
            .scale_psyci_table_ak,
            tables,
            divisors,
            alphas,
            adjustments
        )
    } else {
        lapply(tables, .scale_psyci_table, divisor = scale_info$divisor)
    }
    if (!is.null(scaled_tables)) {
        names(scaled_tables) <- names(tables)
    }

    attr(tables, "scaled_tables") <- scaled_tables
    attr(tables, "scale_info") <- scale_info
    if (!identical(standardized_ci, "bird")) {
        attr(tables, "standardized_ci_method") <- standardized_ci
    }
    class(tables) <- unique(c("psyr_ci", class(tables)))
    tables
}

#' Print PsyR confidence-interval results
#'
#' Prints confidence intervals first in dependent-variable units and then in
#' sample standard-deviation units. By default, the scaled tables divide the
#' estimate, standard error, and confidence limits by the pooled sample SD used
#' by the original Psy program. When Algina-Keselman intervals were requested,
#' the standardized limits instead come from noncentral-t inversion and the
#' scaling divisor is contrast-specific. Inferential columns are unchanged.
#'
#' @param x A result returned by `psyci()`.
#' @param ... Additional arguments passed to the individual table print
#'   methods.
#'
#' @return `x`, invisibly.
#' @export
print.psyr_ci <- function(x, ...) {
    cat("Raw contrast effects and confidence intervals (dependent-variable units)\n")
    .print_psyci_tables(x, ...)

    scale_info <- attr(x, "scale_info", exact = TRUE)
    scaled_tables <- attr(x, "scaled_tables", exact = TRUE)
    standardized_ci <- attr(x, "standardized_ci_method", exact = TRUE)
    if (identical(standardized_ci, "algina_keselman")) {
        cat(paste0(
            "\nScaled contrast effects and Algina-Keselman confidence ",
            "intervals (sample SD units)\n"
        ))
    } else {
        cat("\nScaled contrast effects and confidence intervals (sample SD units)\n")
    }
    if (is.null(scaled_tables) || !isTRUE(scale_info$available)) {
        cat("Unavailable: ", scale_info$reason, ".\n", sep = "")
    } else if (identical(standardized_ci, "algina_keselman")) {
        cat("Scaling divisors are contrast-specific; see each table's scale column.\n")
        .print_psyci_tables(scaled_tables, ...)
    } else {
        cat(
            "Scaling divisor (pooled sample SD): ",
            format(scale_info$divisor, digits = 7),
            "\n",
            sep = ""
        )
        .print_psyci_tables(scaled_tables, ...)
    }

    invisible(x)
}

.print_psyci_tables <- function(tables, ...) {
    table_names <- names(tables)
    for (index in seq_along(tables)) {
        if (length(tables) > 1L || !is.null(table_names)) {
            label <- if (!is.null(table_names) && nzchar(table_names[[index]])) {
                table_names[[index]]
            } else {
                as.character(index)
            }
            cat("\n$", label, "\n", sep = "")
        }
        print(tables[[index]], ...)
    }
    invisible(tables)
}
