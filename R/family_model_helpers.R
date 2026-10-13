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

