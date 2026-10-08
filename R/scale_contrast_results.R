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

.new_psyci_result <- function(tables, model) {
    scale_info <- .psyci_sample_sd(model)
    scaled_tables <- if (isTRUE(scale_info$available)) {
        lapply(tables, .scale_psyci_table, divisor = scale_info$divisor)
    } else {
        NULL
    }
    if (!is.null(scaled_tables)) {
        names(scaled_tables) <- names(tables)
    }

    attr(tables, "scaled_tables") <- scaled_tables
    attr(tables, "scale_info") <- scale_info
    class(tables) <- unique(c("psyr_ci", class(tables)))
    tables
}

#' Print PsyR confidence-interval results
#'
#' Prints confidence intervals first in dependent-variable units and then in
#' sample standard-deviation units. The scaled tables divide the estimate,
#' standard error, and confidence limits by the pooled sample SD used by the
#' original Psy program; inferential columns are unchanged.
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
    cat("\nScaled contrast effects and confidence intervals (sample SD units)\n")
    if (is.null(scaled_tables) || !isTRUE(scale_info$available)) {
        cat("Unavailable: ", scale_info$reason, ".\n", sep = "")
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
