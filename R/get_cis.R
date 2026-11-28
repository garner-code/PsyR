#' Calculate CIs for a contrast table (emmeans object) given a critical constant
#'
#' @param contrast_table An original `emmGrid` contrast table from emmeans.
#' @param cc - a critical constrant, to be used to calculate CIs.
#'
#' @returns contrast_table - a summary table of the contrasts with the CI info added
#' @details Standard equally spaced polynomial trends over three or more levels,
#' including polynomial interactions, retain estimates and tests but receive
#' missing confidence limits. This is new PsyR behavior, not a legacy Psy rule.
#' Summarised inputs are rejected because they lack contrast coefficients.
#' @export
#'
#' @examples
#' data(spacing)
#' # fit model
#' mod <- afex::aov_ez("subj","yield", spacing, within = "spacing", between = "group")
#' sum_emm_btwn <- emmeans::emmeans(mod, "group") # get marginal means
#' con <- list( # define some contrasts
#'             "12vs34" = c(0.5, 0.5, -0.5, -0.5),
#'              "1vs2" = c(1, -1, 0, 0),
#'              "3vs4" = c(0, 0, 1, -1)
#'              )
#' btwn_con <- emmeans::contrast(sum_emm_btwn, con) # get contrast table
#' cc <- cc_ind_t(v_e=12, alpha=0.05) # for example
#' get_cis(btwn_con, cc)
get_cis <- function(contrast_table, cc){
    if (!inherits(contrast_table, "emmGrid")) {
        stop("`contrast_table` must be an original emmGrid object.", call. = FALSE)
    }
    extraction <- .extract_emm_family_data(contrast_table)
    if (extraction$provenance == "unavailable") {
        stop(paste(extraction$issues, collapse = " "), call. = FALSE)
    }
    .get_cis(contrast_table, cc, .polynomial_contrast_rows(extraction))
}

.get_cis <- function(contrast_table, cc, is_polynomial) {
    contrast_table <- summary(contrast_table, infer = c(FALSE, TRUE))
    inherited_limits <- intersect(
        c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL"), names(contrast_table)
    )
    contrast_table[inherited_limits] <- NULL
    keep <- !is_polynomial
    width <- rep(NA_real_, nrow(contrast_table))
    if (any(keep)) {
        width[keep] <- compute_contrast_ci(cc, contrast_table$SE[keep])
    }
    contrast_table$is_polynomial <- is_polynomial
    contrast_table$cc <- cc
    contrast_table$lower <- contrast_table$estimate - width
    contrast_table$upper <- contrast_table$estimate + width
    if (any(is_polynomial)) {
        attr(contrast_table, "mesg") <- c(
            attr(contrast_table, "mesg"),
            "PsyR suppressed confidence intervals for polynomial contrasts (is_polynomial = TRUE)."
        )
    }
    class(contrast_table) <- unique(c("psyr_contrast_table", class(contrast_table)))
    contrast_table
}
