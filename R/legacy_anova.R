.psyci_legacy_grid_key <- function(grid, factors) {
    if (length(factors) == 0L) {
        return(rep("1", nrow(grid)))
    }
    do.call(paste, c(lapply(grid[factors], as.character), list(sep = "\r")))
}

.psyci_canonical_direction <- function(x, tol) {
    scale <- sqrt(sum(x^2))
    if (!is.finite(scale) || scale <= 0) {
        stop("A recovered legacy contrast direction was empty.", call. = FALSE)
    }
    result <- x / scale
    first <- which(abs(result) > tol)[1L]
    if (!is.na(first) && result[[first]] < 0) {
        result <- -result
    }
    unname(result)
}

.psyci_direction_id <- function(direction, directions, tol) {
    if (length(directions) == 0L) {
        return(NA_integer_)
    }
    matches <- vapply(directions, function(candidate) {
        max(abs(direction - candidate)) <= tol
    }, logical(1))
    match(TRUE, matches, nomatch = NA_integer_)
}

.psyci_add_direction <- function(direction, directions, tol) {
    id <- .psyci_direction_id(direction, directions, tol)
    if (!is.na(id)) {
        return(list(id = id, directions = directions))
    }
    directions[[length(directions) + 1L]] <- direction
    list(id = length(directions), directions = directions)
}

.psyci_legacy_design <- function(model, roles, reference_grid) {
    long_data <- model$data$long
    response <- attr(model, "dv", exact = TRUE)
    subject_id <- attr(model, "id", exact = TRUE)
    between <- roles$between
    within <- roles$within
    required <- unique(c(subject_id, response, between, within))
    missing_columns <- setdiff(required, names(long_data))
    if (length(missing_columns) > 0L) {
        stop(
            sprintf(
                "The afex model data omitted columns required for the legacy ANOVA table: %s.",
                paste(missing_columns, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    if (length(between) == 0L) {
        group_grid <- data.frame(.group = "1", stringsAsFactors = FALSE)
        group_levels <- "1"
        row_group <- rep("1", nrow(long_data))
    } else {
        group_grid <- unique(reference_grid[, between, drop = FALSE])
        group_levels <- .psyci_legacy_grid_key(group_grid, between)
        row_group <- .psyci_legacy_grid_key(long_data, between)
    }
    if (length(within) == 0L) {
        repeat_grid <- data.frame(.repeat = "1", stringsAsFactors = FALSE)
        repeat_levels <- "1"
        row_repeat <- rep("1", nrow(long_data))
    } else {
        repeat_grid <- unique(reference_grid[, within, drop = FALSE])
        repeat_levels <- .psyci_legacy_grid_key(repeat_grid, within)
        row_repeat <- .psyci_legacy_grid_key(long_data, within)
    }

    subjects <- unique(as.character(long_data[[subject_id]]))
    wide <- matrix(
        NA_real_,
        nrow = length(subjects),
        ncol = length(repeat_levels),
        dimnames = list(subjects, repeat_levels)
    )
    groups <- integer(length(subjects))
    response_values <- long_data[[response]]
    for (index in seq_along(subjects)) {
        rows <- which(as.character(long_data[[subject_id]]) == subjects[[index]])
        subject_groups <- unique(match(row_group[rows], group_levels))
        if (length(subject_groups) != 1L || is.na(subject_groups)) {
            stop(
                "Each subject must belong to exactly one between-subject cell for legacy ANOVA calculations.",
                call. = FALSE
            )
        }
        repeated_cells <- match(row_repeat[rows], repeat_levels)
        if (anyNA(repeated_cells) || anyDuplicated(repeated_cells)) {
            stop(
                "Each subject must have one observation in every repeated-measure cell for legacy ANOVA calculations.",
                call. = FALSE
            )
        }
        if (length(repeated_cells) != length(repeat_levels) ||
            anyNA(response_values[rows]) ||
            any(!is.finite(response_values[rows]))) {
            stop(
                "Legacy ANOVA calculations require complete finite repeated-measure observations for every subject.",
                call. = FALSE
            )
        }
        groups[[index]] <- subject_groups
        wide[index, repeated_cells] <- response_values[rows]
    }
    if (anyNA(wide)) {
        stop(
            "Legacy ANOVA calculations require complete repeated-measure data.",
            call. = FALSE
        )
    }

    n_group <- tabulate(groups, nbins = nrow(group_grid))
    if (any(n_group == 0L)) {
        stop("Every between-subject cell must contain at least one subject.", call. = FALSE)
    }
    dfe <- nrow(wide) - nrow(group_grid)
    if (dfe <= 0L) {
        stop("The legacy ANOVA error degrees of freedom must be positive.", call. = FALSE)
    }

    list(
        data = wide,
        group = groups,
        n_group = n_group,
        dfe = dfe,
        between = between,
        within = within,
        group_grid = group_grid,
        repeat_grid = repeat_grid,
        group_levels = group_levels,
        repeat_levels = repeat_levels
    )
}

.psyci_legacy_coefficient_matrix <- function(extraction, design, column) {
    grid <- extraction$grid
    group_key <- if (length(design$between) == 0L) {
        rep("1", nrow(grid))
    } else {
        .psyci_legacy_grid_key(grid, design$between)
    }
    repeat_key <- if (length(design$within) == 0L) {
        rep("1", nrow(grid))
    } else {
        .psyci_legacy_grid_key(grid, design$within)
    }
    group_index <- match(group_key, design$group_levels)
    repeat_index <- match(repeat_key, design$repeat_levels)
    if (anyNA(group_index) || anyNA(repeat_index) ||
        anyDuplicated(paste(group_index, repeat_index, sep = "\r"))) {
        stop(
            "PsyR could not map a contrast coefficient grid uniquely to the legacy design cells.",
            call. = FALSE
        )
    }
    matrix_values <- matrix(
        0,
        nrow = length(design$group_levels),
        ncol = length(design$repeat_levels)
    )
    matrix_values[cbind(group_index, repeat_index)] <-
        extraction$coefficients[, column]
    matrix_values
}

.psyci_legacy_directions <- function(matrix_values, family, tol) {
    if (family == "b") {
        a <- rowMeans(matrix_values)
        w <- rep(1, ncol(matrix_values))
    } else if (family == "w") {
        a <- rep(1, nrow(matrix_values))
        w <- colMeans(matrix_values)
    } else {
        decomposition <- svd(matrix_values)
        leading <- decomposition$d[[1L]]
        remainder <- if (length(decomposition$d) > 1L) {
            max(decomposition$d[-1L])
        } else {
            0
        }
        if (!is.finite(leading) || leading <= 0 ||
            remainder > tol * leading) {
            stop(
                "A between-by-within contrast could not be factored into one between and one within direction.",
                call. = FALSE
            )
        }
        a <- decomposition$u[, 1L] * leading
        w <- decomposition$v[, 1L]
    }
    list(
        a = unname(a),
        w = unname(w),
        canonical_a = .psyci_canonical_direction(a, tol),
        canonical_w = .psyci_canonical_direction(w, tol)
    )
}

.psyci_legacy_row_statistics <- function(design, a, w) {
    psy_tolerance <- 1e-30
    w_norm <- sqrt(sum(w^2))
    if (!is.finite(w_norm) || w_norm <= 0) {
        stop("A recovered within direction had zero length.", call. = FALSE)
    }
    scores <- as.vector(design$data %*% (w / w_norm))
    fitted <- ave(scores, design$group, FUN = mean)
    error_ss <- sum((scores - fitted)^2)
    if (error_ss < psy_tolerance) {
        error_ss <- 0
    }
    error_ms <- error_ss / design$dfe
    group_means <- vapply(seq_along(design$n_group), function(group) {
        mean(scores[design$group == group])
    }, numeric(1))
    sample_value <- sum(a * group_means)
    denominator <- sum(a^2 / design$n_group)
    contrast_ss <- if (abs(denominator) < psy_tolerance) {
        0
    } else {
        sample_value^2 / denominator
    }
    if (abs(contrast_ss) < psy_tolerance) {
        contrast_ss <- 0
    }
    f_value <- if (error_ms < psy_tolerance) {
        Inf
    } else {
        contrast_ss / error_ms
    }
    list(
        SS = contrast_ss,
        df = 1,
        MS = contrast_ss,
        F_value = f_value,
        error_SS = error_ss,
        error_df = design$dfe,
        error_MS = error_ms
    )
}

.psyci_legacy_statistics <- function(
    model,
    contrast_tables,
    extractions,
    families,
    roles,
    tol
) {
    reference_grid <- extractions[[1L]]$grid
    design <- .psyci_legacy_design(model, roles, reference_grid)
    items <- list()
    f_values <- vector("list", length(extractions))
    between_directions <- list()
    within_directions <- list()

    for (table_index in seq_along(extractions)) {
        extraction <- extractions[[table_index]]
        family <- families[[table_index]]
        summary_table <- as.data.frame(summary(
            contrast_tables[[table_index]], infer = c(FALSE, FALSE)
        ))
        f_values[[table_index]] <- rep(NA_real_, ncol(extraction$coefficients))
        for (row_index in seq_len(ncol(extraction$coefficients))) {
            coefficient_matrix <- .psyci_legacy_coefficient_matrix(
                extraction,
                design,
                row_index
            )
            directions <- .psyci_legacy_directions(
                coefficient_matrix,
                family,
                tol
            )
            between_id <- NA_integer_
            within_id <- NA_integer_
            if (family %in% c("b", "bw")) {
                addition <- .psyci_add_direction(
                    directions$canonical_a,
                    between_directions,
                    tol
                )
                between_id <- addition$id
                between_directions <- addition$directions
            }
            if (family %in% c("w", "bw")) {
                addition <- .psyci_add_direction(
                    directions$canonical_w,
                    within_directions,
                    tol
                )
                within_id <- addition$id
                within_directions <- addition$directions
            }
            statistics <- .psyci_legacy_row_statistics(
                design,
                directions$a,
                directions$w
            )
            estimable <- row_index <= nrow(summary_table) &&
                !is.na(summary_table$estimate[[row_index]]) &&
                !is.na(summary_table$SE[[row_index]])
            if (estimable) {
                f_values[[table_index]][[row_index]] <- statistics$F_value
            }
            items[[length(items) + 1L]] <- list(
                table = table_index,
                row = row_index,
                order = length(items) + 1L,
                family = family,
                label = extraction$contrasts$label[[row_index]],
                between_id = between_id,
                within_id = within_id,
                statistics = statistics,
                estimable = estimable
            )
        }
    }

    list(
        dfe = design$dfe,
        f_values = f_values,
        table = .psyci_build_legacy_anova_table(items, design)
    )
}

.psyci_legacy_table_row <- function(
    section,
    label,
    source,
    row_type,
    SS,
    df,
    MS,
    F_value,
    table = NA_integer_,
    row = NA_integer_
) {
    data.frame(
        section = section,
        label = label,
        source = source,
        row_type = row_type,
        SS = SS,
        df = df,
        MS = MS,
        F_value = F_value,
        table = table,
        row = row,
        stringsAsFactors = FALSE
    )
}

.psyci_build_legacy_anova_table <- function(items, design) {
    rows <- list()
    between_items <- Filter(function(item) item$family == "b", items)
    if (length(between_items) > 0L) {
        for (item in between_items) {
            statistics <- item$statistics
            rows[[length(rows) + 1L]] <- .psyci_legacy_table_row(
                "Between",
                item$label,
                paste0("B", item$between_id),
                "contrast",
                statistics$SS,
                1,
                statistics$MS,
                if (item$estimable) statistics$F_value else NA_real_,
                item$table, item$row
            )
        }
        error <- between_items[[1L]]$statistics
        rows[[length(rows) + 1L]] <- .psyci_legacy_table_row(
            "Between", "", "Error", "error",
            error$error_SS, design$dfe, error$error_MS, NA_real_
        )
    }

    within_items <- Filter(function(item) item$family %in% c("w", "bw"), items)
    if (length(within_items) > 0L) {
        within_ids <- unique(vapply(within_items, `[[`, integer(1), "within_id"))
        for (within_id in within_ids) {
            group <- Filter(function(item) item$within_id == within_id, within_items)
            main_effects <- Filter(function(item) item$family == "w", group)
            interactions <- Filter(function(item) item$family == "bw", group)
            for (item in c(main_effects, interactions)) {
                statistics <- item$statistics
                source <- if (item$family == "w") {
                    paste0("W", item$within_id)
                } else {
                    paste0("B", item$between_id, "W", item$within_id)
                }
                rows[[length(rows) + 1L]] <- .psyci_legacy_table_row(
                    "Within",
                    item$label,
                    source,
                    "contrast",
                    statistics$SS,
                    1,
                    statistics$MS,
                    if (item$estimable) statistics$F_value else NA_real_,
                    item$table, item$row
                )
            }
            error <- group[[1L]]$statistics
            rows[[length(rows) + 1L]] <- .psyci_legacy_table_row(
                "Within", "", "Error", "error",
                error$error_SS, design$dfe, error$error_MS, NA_real_
            )
        }
    }

    result <- if (length(rows) == 0L) {
        .psyci_legacy_table_row(
            character(), character(), character(), character(),
            numeric(), integer(), numeric(), numeric(), integer(), integer()
        )
    } else {
        do.call(rbind, rows)
    }
    rownames(result) <- NULL
    source_rows <- result[c("table", "row")]
    result[c("table", "row")] <- NULL
    attr(result, "psyr_source_rows") <- source_rows
    class(result) <- c("psyr_legacy_anova", "data.frame")
    result
}

.finalize_psyci_table <- function(table) {
    table$F_crit <- table$cc^2
    table$cc <- NULL
    if ("theta" %in% names(table)) {
        parameters <- attr(table, "psyr_probability_parameters", exact = TRUE)
        parameters$theta <- table$theta
        attr(table, "psyr_probability_parameters") <- parameters
        table$theta <- NULL
    }
    saved_attributes <- attributes(table)
    if ("df" %in% names(table)) {
        table[["df_error"]] <- table[["df"]]
    }
    table[["df"]] <- NULL
    table[["t.ratio"]] <- NULL
    table[["p.value"]] <- NULL

    inferential <- c(
        "estimate", "SE", "df_num", "df_error", "F_value", "F_crit",
        "p_value", "p_value_adjusted", "p_adjustment", "sig", "lower",
        "upper"
    )
    aliases <- c(df_num = "df1", df_error = "df2",
                 p_value_adjusted = "p_adj", p_adjustment = "method")
    matched <- match(names(table), names(aliases))
    names(table)[!is.na(matched)] <- unname(aliases[matched[!is.na(matched)]])
    inferential <- ifelse(inferential %in% names(aliases),
                          aliases[inferential], inferential)
    identifiers <- setdiff(names(table), inferential)
    table <- table[c(identifiers, intersect(inferential, names(table)))]

    retained_attributes <- setdiff(
        names(saved_attributes),
        c("names", "row.names", "adjust", "mesg")
    )
    for (attribute in retained_attributes) {
        attr(table, attribute) <- saved_attributes[[attribute]]
    }
    class(table) <- c("psyr_contrast_table", "data.frame")

    messages <- saved_attributes$mesg
    if (!is.null(messages)) {
        messages <- messages[!startsWith(messages, "P value adjustment:")]
        attr(table, "mesg") <- messages
    }
    attr(table, "adjust") <- NULL
    table
}

#' Print a PsyR contrast table
#'
#' Numeric values print to three decimal places, with integer degrees of
#' freedom and blank unavailable values. Polynomial indicators and the
#' unadjusted p-value are hidden. Interaction component labels are combined
#' into a single printed `contrast` column when their metadata is available.
#' Stored column names `df1`, `df2`, `p_adj`, and `method` print unchanged.
#'
#' @param x A contrast table returned by [psyci()].
#' @param ... Additional arguments passed to [print.data.frame()].
#'
#' @return `x`, invisibly.
#' @export
print.psyr_contrast_table <- function(x, ...) {
    display <- .psyci_display_table(x)
    print.data.frame(display, ...)
    messages <- attr(x, "mesg", exact = TRUE)
    if (length(messages) > 0L) {
        cat("\n", paste(messages, collapse = "\n"), "\n", sep = "")
    }
    invisible(x)
}

.psyci_format_values <- function(value, integer = FALSE) {
    display <- if (is.numeric(value) && !integer) {
        formatC(value, format = "f", digits = 3)
    } else {
        as.character(value)
    }
    display[is.na(value)] <- ""
    display
}

.psyci_display_table <- function(x) {
    display <- as.data.frame(x)
    hidden <- c("p_value", "p.value", "is_polynomial", "is.polynomial")
    if (!"contrast" %in% names(display)) {
        components <- attr(x, "pri.vars", exact = TRUE)
        by_columns <- attr(x, "by.vars", exact = TRUE)
        if (length(components) > 0L &&
            all(components %in% setdiff(names(display), c(hidden, by_columns)))) {
            labels <- lapply(display[components], function(value) {
                paste0("(", as.character(value), ")")
            })
            display$contrast <- do.call(paste, c(labels, sep = " \u00d7 "))
            display <- display[c("contrast", setdiff(names(display),
                                                     c("contrast", components)))]
        }
    }
    display[intersect(hidden, names(display))] <- NULL
    for (column in names(display)) {
        display[[column]] <- .psyci_format_values(
            display[[column]], integer = column %in% c("df1", "df2", "df")
        )
    }
    display
}

#' Print a Psy legacy ANOVA contrast table
#'
#' Displays `F_crit` and `p_adj` from the selected PsyR procedure alongside
#' the legacy statistics. Inference cells for error rows are blank.
#'
#' @param x A legacy ANOVA table attached to a result from [psyci()].
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.psyr_legacy_anova <- function(x, ...) {
    cat("Analysis of Variance Summary Table\n")
    if (nrow(x) == 0L) {
        cat("No supplied contrasts were available.\n")
        return(invisible(x))
    }
    format_value <- function(value, infinite = FALSE) {
        if (is.na(value)) {
            return("")
        }
        if (is.infinite(value)) {
            return(if (infinite) "Infinite" else as.character(value))
        }
        formatC(value, format = "f", digits = 3)
    }
    for (section in unique(x$section)) {
        cat("\n", section, "\n", sep = "")
        selected <- x[x$section == section, , drop = FALSE]
        inference <- attr(x, "psyr_inference", exact = TRUE)
        display <- data.frame(
            Label = selected$label,
            Source = selected$source,
            SS = vapply(selected$SS, format_value, character(1)),
            df = as.character(selected$df),
            MS = vapply(selected$MS, format_value, character(1)),
            F = vapply(
                selected$F_value,
                format_value,
                character(1),
                infinite = TRUE
            ),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        if (!is.null(inference)) {
            rows <- x$section == section
            display$F_crit <- .psyci_format_values(inference$F_crit[rows])
            display$p_adj <- .psyci_format_values(inference$p_adj[rows])
        }
        print(display, row.names = FALSE, right = TRUE, quote = FALSE)
    }
    invisible(x)
}
