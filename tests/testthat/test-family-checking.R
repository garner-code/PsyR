make_family_between_fixture <- function() {
    data("fatigue", package = "PsyR", envir = environment())

    fatigue$Group <- factor(
        fatigue$Group,
        levels = c(4, 1, 3, 2),
        labels = c("Extreme", "Alert", "Moderate", "Mild")
    )

    model <- afex::aov_ez(
        id = "id",
        dv = "Errors",
        data = fatigue,
        between = "Group"
    )
    emmeans <- emmeans::emmeans(model, "Group")
    contrast_table <- emmeans::contrast(
        emmeans,
        method = list("Extreme - Alert" = c(1, -1, 0, 0)),
        adjust = "none"
    )

    list(
        model = model,
        emmeans = emmeans,
        table = contrast_table,
        factor_order = c("Extreme", "Alert", "Moderate", "Mild")
    )
}

make_family_within_fixture <- function() {
    data("priming", package = "PsyR", envir = environment())

    priming$Subject <- factor(priming$Subject)
    priming$Prime <- factor(
        priming$Prime,
        levels = c("Stereotypical", "Atypical", "Neutral")
    )

    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    afex::afex_options(emmeans_model = "multivariate")

    model <- afex::aov_ez(
        id = "Subject",
        dv = "ReactionTime",
        data = priming,
        within = "Prime"
    )
    emmeans <- emmeans::emmeans(model, "Prime")
    contrast_table <- emmeans::contrast(
        emmeans,
        method = list(
            "Stereotypical - Others" = c(-1, 0.5, 0.5),
            "Atypical - Neutral" = c(0, 1, -1)
        ),
        adjust = "none"
    )

    list(
        model = model,
        emmeans = emmeans,
        table = contrast_table,
        factor_order = c("Stereotypical", "Atypical", "Neutral")
    )
}

make_family_mixed_fixture <- function() {
    data("spacing", package = "PsyR", envir = environment())

    spacing$group <- factor(spacing$group, levels = c(1, 2, 3, 4))
    spacing$spacing <- factor(
        spacing$spacing,
        levels = c("TWENTY", "FORTY", "SIXTY")
    )

    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    afex::afex_options(emmeans_model = "multivariate")

    model <- afex::aov_ez(
        id = "subj",
        dv = "yield",
        data = spacing,
        between = "group",
        within = "spacing"
    )

    between_contrasts <- list("Group 1 - Group 2" = c(1, -1, 0, 0))
    within_contrasts <- list("Twenty - Forty" = c(1, -1, 0))

    between_emmeans <- emmeans::emmeans(model, "group")
    within_emmeans <- emmeans::emmeans(model, "spacing")
    interaction_emmeans <- emmeans::emmeans(model, c("group", "spacing"))

    between_table <- emmeans::contrast(
        between_emmeans,
        method = between_contrasts,
        adjust = "none"
    )
    within_table <- emmeans::contrast(
        within_emmeans,
        method = within_contrasts,
        adjust = "none"
    )
    interaction_table <- emmeans::contrast(
        interaction_emmeans,
        interaction = list(between_contrasts, within_contrasts),
        adjust = "none"
    )
    within_by_between <- emmeans::contrast(
        emmeans::emmeans(model, ~ spacing | group),
        method = within_contrasts,
        adjust = "none"
    )
    between_by_within <- emmeans::contrast(
        emmeans::emmeans(model, ~ group | spacing),
        method = between_contrasts,
        adjust = "none"
    )
    within_at_between <- emmeans::contrast(
        emmeans::emmeans(
            model,
            "spacing",
            at = list(group = "1")
        ),
        method = within_contrasts,
        adjust = "none"
    )

    list(
        model = model,
        between = between_table,
        within = within_table,
        interaction = interaction_table,
        within_by_between = within_by_between,
        between_by_within = between_by_within,
        within_at_between = within_at_between
    )
}

test_that("between-subject contrasts are accepted without reordering levels", {
    fixture <- make_family_between_fixture()
    coefficients_before <- stats::coef(fixture$table)

    expect_identical(
        as.character(coefficients_before$Group),
        fixture$factor_order
    )
    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = list(fixture$table),
            family_list = list("b"),
            between_factors = list("Group"),
            mode = "error"
        )
    )
    expect_s3_class(check, "psyr_family_check")
    expect_true(check$valid)
    expect_identical(check$tables[[1]]$detected_family, "b")

    coefficients_after <- stats::coef(fixture$table)
    expect_identical(coefficients_after, coefficients_before)
})

test_that("within-subject contrasts use the actual priming columns and order", {
    fixture <- make_family_within_fixture()
    coefficients_before <- stats::coef(fixture$table)

    expect_identical(
        as.character(coefficients_before$Prime),
        fixture$factor_order
    )
    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = list(fixture$table),
            family_list = list("w"),
            within_factors = list("Prime"),
            mode = "error"
        )
    )
    expect_s3_class(check, "psyr_family_check")
    expect_true(check$valid)
    expect_identical(check$tables[[1]]$detected_family, "w")
    expect_identical(stats::coef(fixture$table), coefficients_before)
})

test_that("standard mixed-design b, w, and bw families are accepted", {
    fixture <- make_family_mixed_fixture()

    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = list(
                fixture$between,
                fixture$within,
                fixture$interaction
            ),
            family_list = list("b", "w", "bw"),
            between_factors = list("group"),
            within_factors = list("spacing"),
            mode = "error"
        )
    )
    expect_s3_class(check, "psyr_family_check")
    expect_true(check$valid)
    expect_identical(
        vapply(
            check$tables,
            function(table) table$detected_family,
            character(1)
        ),
        c("b", "w", "bw")
    )
})

test_that("conditioned mixed-design contrasts are classified as bw", {
    fixture <- make_family_mixed_fixture()

    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = list(
                fixture$within_by_between,
                fixture$between_by_within
            ),
            family_list = list("bw", "bw"),
            between_factors = list("group"),
            within_factors = list("spacing"),
            mode = "error"
        )
    )
    expect_true(check$valid)
    expect_identical(
        vapply(
            check$tables,
            function(table) table$detected_family,
            character(1)
        ),
        c("bw", "bw")
    )
})

test_that("clear family mismatches respect error, warn, and none modes", {
    fixture <- make_family_mixed_fixture()
    args <- list(
        model = fixture$model,
        contrast_tables = list(fixture$between),
        family_list = list("w"),
        between_factors = list("group"),
        within_factors = list("spacing")
    )

    expect_error(
        do.call(check_contrast_families, c(args, list(mode = "error"))),
        regexp = "family|declared|detected",
        ignore.case = TRUE
    )
    expect_warning(
        do.call(check_contrast_families, c(args, list(mode = "warn"))),
        regexp = "family|declared|detected",
        ignore.case = TRUE
    )
    none_report <- expect_silent(
        do.call(check_contrast_families, c(args, list(mode = "none")))
    )
    expect_false(none_report$valid)
    expect_identical(none_report$summary$status, "mismatch")
    expect_true(any(grepl("declared|detected", none_report$issues)))
    expect_true(any(grepl(
        "declared|detected",
        none_report$tables[[1L]]$issues
    )))
})

test_that("fixed factors hidden by at are recovered as conditioning", {
    fixture <- make_family_mixed_fixture()
    report <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = list(fixture$within_at_between),
            family_list = list("bw"),
            between_factors = list("group"),
            within_factors = list("spacing"),
            mode = "error"
        )
    )

    expect_true(report$valid)
    expect_identical(report$tables[[1L]]$detected_family, "bw")
    expect_identical(report$tables[[1L]]$scope, "conditioned")
    expect_identical(
        report$tables[[1L]]$extraction$provenance,
        "reconstructed_full_grid"
    )
})

test_that("averaging over a subset of another domain is unverifiable", {
    fixture <- make_family_mixed_fixture()
    subset_emmeans <- emmeans::emmeans(
        fixture$model,
        "spacing",
        at = list(group = c("1", "2"))
    )
    subset_table <- emmeans::contrast(
        subset_emmeans,
        method = list("Twenty - Forty" = c(1, -1, 0)),
        adjust = "none"
    )

    report <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = subset_table,
            family_list = list("w"),
            mode = "none"
        )
    )

    expect_identical(report$tables[[1L]]$detected_family, "unknown")
    expect_identical(report$summary$status, "unverifiable")
    expect_identical(
        report$tables[[1L]]$extraction$provenance,
        "reconstructed_full_grid"
    )
})

test_that("full-grid recovery respects the emmeans coordinate basis", {
    data("spacing", package = "PsyR", envir = environment())
    spacing$group <- factor(spacing$group)
    model <- afex::aov_ez(
        id = "subj",
        dv = "yield",
        data = spacing,
        between = "group",
        within = "spacing",
        include_aov = TRUE
    )
    full_table <- emmeans::contrast(
        emmeans::emmeans(model, "spacing", model = "univariate"),
        method = list("Twenty - Forty" = c(1, -1, 0)),
        adjust = "none"
    )
    subset_table <- suppressWarnings(emmeans::contrast(
        emmeans::emmeans(
            model,
            "spacing",
            at = list(group = c("1", "2")),
            model = "univariate"
        ),
        method = list("Twenty - Forty" = c(1, -1, 0)),
        adjust = "none"
    ))

    full_report <- expect_silent(
        check_contrast_families(
            model = model,
            contrast_tables = full_table,
            family_list = list("w"),
            mode = "error"
        )
    )
    subset_report <- expect_silent(
        check_contrast_families(
            model = model,
            contrast_tables = subset_table,
            family_list = list("w"),
            mode = "none"
        )
    )

    expect_identical(full_report$tables[[1L]]$detected_family, "w")
    expect_identical(
        full_report$tables[[1L]]$extraction$provenance,
        "reconstructed_full_grid"
    )
    expect_identical(subset_report$tables[[1L]]$detected_family, "unknown")
})

test_that("failed full-grid recovery cannot pass from a reduced grid", {
    fixture <- make_family_mixed_fixture()
    subset_emmeans <- emmeans::emmeans(
        fixture$model,
        "spacing",
        at = list(group = c("1", "2"))
    )
    regridded <- emmeans::regrid(subset_emmeans)
    subset_table <- emmeans::contrast(
        regridded,
        method = list("Twenty - Forty" = c(1, -1, 0)),
        adjust = "none"
    )

    report <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = subset_table,
            family_list = list("w"),
            mode = "none"
        )
    )

    expect_identical(report$tables[[1L]]$detected_family, "unknown")
    expect_identical(report$tables[[1L]]$extraction$provenance, "unavailable")
    expect_identical(report$summary$status, "unverifiable")
    expect_true(any(grepl(
        "incompatible|recover|reconstruct",
        report$issues,
        ignore.case = TRUE
    )))
})

test_that("summary_emm inputs follow the requested validation mode", {
    fixture <- make_family_between_fixture()
    summary_table <- summary(fixture$table)
    args <- list(
        model = fixture$model,
        contrast_tables = list(summary_table),
        family_list = list("b"),
        between_factors = list("Group")
    )

    expect_error(
        do.call(check_contrast_families, c(args, list(mode = "error"))),
        regexp = "summary_emm|coefficient|emmGrid|verify",
        ignore.case = TRUE
    )
    expect_warning(
        do.call(check_contrast_families, c(args, list(mode = "warn"))),
        regexp = "summary_emm|coefficient|emmGrid|verify",
        ignore.case = TRUE
    )
    none_report <- expect_silent(
        do.call(check_contrast_families, c(args, list(mode = "none")))
    )
    expect_identical(
        none_report$tables[[1L]]$contrasts$label,
        as.character(summary_table$contrast)
    )
})

test_that("a scalar family code is recycled across compatible tables", {
    fixture <- make_family_between_fixture()
    second_table <- emmeans::contrast(
        fixture$emmeans,
        method = list("Moderate - Mild" = c(0, 0, 1, -1)),
        adjust = "none"
    )

    expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = list(fixture$table, second_table),
            family_list = list("b"),
            between_factors = list("Group"),
            mode = "error"
        )
    )
})

test_that("family checking does not change valid psyci results", {
    fixture <- make_family_between_fixture()
    psyci_args <- list(
        model = fixture$model,
        contrast_tables = list(fixture$table),
        method = "ph",
        family_list = list("b"),
        between_factors = list("Group")
    )

    result_without_check <- do.call(
        psyci,
        c(psyci_args, list(family_check = "none"))
    )
    result_with_check <- do.call(
        psyci,
        c(psyci_args, list(family_check = "error"))
    )

    expect_equal(result_with_check, result_without_check, tolerance = 1e-12)
})

test_that("checking preserves mixed-design post-hoc results", {
    fixture <- make_family_mixed_fixture()
    psyci_args <- list(
        model = fixture$model,
        contrast_tables = list(
            fixture$between,
            fixture$within,
            fixture$interaction
        ),
        method = "ph",
        family_list = list("b", "w", "bw"),
        between_factors = list("group"),
        within_factors = list("spacing")
    )

    result_without_check <- do.call(
        psyci,
        c(psyci_args, list(family_check = "none"))
    )
    result_with_check <- do.call(
        psyci,
        c(psyci_args, list(family_check = "error"))
    )

    expect_equal(result_with_check, result_without_check, tolerance = 1e-12)
})

test_that("afex roles and emmeans coefficient data are extracted explicitly", {
    fixture <- make_family_mixed_fixture()
    roles <- PsyR:::.get_afex_factor_roles(fixture$model)

    expect_identical(roles$between, "group")
    expect_identical(roles$within, "spacing")
    expect_true(roles$verifiable)

    extracted <- PsyR:::.extract_emm_family_data(fixture$interaction)
    expect_identical(extracted$provenance, "complete")
    expect_true(all(c("group", "spacing") %in% names(extracted$grid)))
    expect_true(is.matrix(extracted$coefficients))
    expect_identical(
        ncol(extracted$coefficients),
        nrow(as.data.frame(summary(fixture$interaction)))
    )
})

test_that("one extracted contrast remains a matrix", {
    fixture <- make_family_between_fixture()
    coefficients <- get_contrast_coefficients(fixture$table)

    expect_length(coefficients, 1L)
    expect_true(is.matrix(coefficients[[1L]]))
    expect_identical(ncol(coefficients[[1L]]), 1L)
})

test_that("the checker trusts coefficient validity while checking its domain", {
    fixture <- make_family_between_fixture()
    nonzero_sum <- emmeans::contrast(
        fixture$emmeans,
        method = list("Trusted coefficients" = c(1, -0.5, 0, 0)),
        adjust = "none"
    )

    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = nonzero_sum,
            family_list = list("b"),
            between_factors = list("Group"),
            mode = "error"
        )
    )

    expect_true(check$valid)
    expect_identical(check$tables[[1L]]$detected_family, "b")
})

test_that("a full mixed grid can still contain a marginal b contrast", {
    fixture <- make_family_mixed_fixture()
    cell_means <- emmeans::emmeans(
        fixture$model,
        c("group", "spacing")
    )
    grid <- cell_means@grid
    coefficient <- ifelse(
        grid$group == "1",
        1 / 3,
        ifelse(grid$group == "2", -1 / 3, 0)
    )
    marginal_b <- emmeans::contrast(
        cell_means,
        method = list("Group 1 - Group 2" = coefficient),
        adjust = "none"
    )

    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = marginal_b,
            family_list = list("b"),
            mode = "error"
        )
    )

    expect_identical(check$tables[[1L]]$detected_family, "b")
    expect_identical(check$tables[[1L]]$scope, "marginal")
})

test_that("unequal averaging weights remain marginal to the other domain", {
    fixture <- make_family_mixed_fixture()
    cell_means <- emmeans::emmeans(
        fixture$model,
        c("group", "spacing")
    )
    grid <- cell_means@grid
    spacing_weights <- c(TWENTY = 0.5, FORTY = 0.3, SIXTY = 0.2)
    group_sign <- ifelse(
        grid$group == "1",
        1,
        ifelse(grid$group == "2", -1, 0)
    )
    coefficient <- group_sign * spacing_weights[as.character(grid$spacing)]
    weighted_marginal <- emmeans::contrast(
        cell_means,
        method = list("Weighted group difference" = coefficient),
        adjust = "none"
    )

    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = weighted_marginal,
            family_list = list("b"),
            mode = "error"
        )
    )

    expect_identical(check$tables[[1L]]$detected_family, "b")
    expect_identical(check$tables[[1L]]$scope, "marginal")
    expect_identical(
        unname(check$tables[[1L]]$participation[[1L]]$factors$status),
        c("varying", "averaged")
    )
})

test_that("non-common averaging profiles reveal both subject domains", {
    fixture <- make_family_mixed_fixture()
    cell_means <- emmeans::emmeans(
        fixture$model,
        c("group", "spacing")
    )
    grid <- cell_means@grid
    coefficient <- numeric(nrow(grid))
    coefficient[grid$group == "1"] <- c(1, 2, 3)
    coefficient[grid$group == "2"] <- c(-3, -2, -1)
    nonseparable <- emmeans::contrast(
        cell_means,
        method = list("Non-common profiles" = coefficient),
        adjust = "none"
    )

    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = nonseparable,
            family_list = list("bw"),
            mode = "error"
        )
    )

    expect_identical(check$tables[[1L]]$detected_family, "bw")
    expect_identical(
        unname(check$tables[[1L]]$participation[[1L]]$factors$status),
        c("varying", "varying")
    )
})

test_that("structural-zero custom averaging is reported as ambiguous", {
    fixture <- make_family_mixed_fixture()
    cell_means <- emmeans::emmeans(
        fixture$model,
        c("group", "spacing")
    )
    grid <- cell_means@grid
    spacing_weights <- c(TWENTY = 0.5, FORTY = 0.5, SIXTY = 0)
    group_sign <- ifelse(
        grid$group == "1",
        1,
        ifelse(grid$group == "2", -1, 0)
    )
    coefficient <- group_sign * spacing_weights[as.character(grid$spacing)]
    partial_average <- emmeans::contrast(
        cell_means,
        method = list("Partial average" = coefficient),
        adjust = "none"
    )

    report <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = partial_average,
            family_list = list("b"),
            mode = "none"
        )
    )

    expect_identical(report$tables[[1L]]$detected_family, "unknown")
    expect_identical(report$summary$status, "unverifiable")
    expect_true(any(grepl("structural zeros", report$issues)))
})

test_that("one-level selectors are classified as conditioning", {
    fixture <- make_family_mixed_fixture()
    cell_means <- emmeans::emmeans(
        fixture$model,
        c("group", "spacing")
    )
    grid <- cell_means@grid
    group_contrast <- ifelse(
        grid$group == "1",
        1,
        ifelse(grid$group == "2", -1, 0)
    )
    spacing_selector <- as.numeric(grid$spacing == "FORTY")
    simple_effect <- emmeans::contrast(
        cell_means,
        method = list(
            "Group 1 - Group 2 at Forty" =
                group_contrast * spacing_selector
        ),
        adjust = "none"
    )

    report <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = simple_effect,
            family_list = list("bw"),
            mode = "error"
        )
    )

    expect_identical(report$tables[[1L]]$detected_family, "bw")
    expect_identical(report$tables[[1L]]$scope, "conditioned")
})

test_that("joint nuisance weights need not factor across nuisance factors", {
    grid <- expand.grid(
        group = c("g1", "g2"),
        B = c("b1", "b2"),
        C = c("c1", "c2"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    roles <- list(
        between = "group",
        within = c("B", "C"),
        levels = list(
            group = c("g1", "g2"),
            B = c("b1", "b2"),
            C = c("c1", "c2")
        )
    )
    joint_key <- interaction(grid$B, grid$C, lex.order = TRUE)
    joint_weights <- c(0.1, 0.2, 0.3, 0.4)
    coefficient <- ifelse(grid$group == "g1", 1, -1) *
        joint_weights[as.integer(joint_key)]
    participation <- PsyR:::.detect_factor_participation(
        coefficients = matrix(coefficient, ncol = 1L),
        grid = grid,
        roles = roles
    )[[1L]]
    classification <- PsyR:::.classify_contrast_family(
        participation,
        roles
    )

    expect_identical(classification$family, "b")
    expect_false(any(participation$factors$participates[
        participation$factors$role == "within"
    ]))
})

test_that("same-domain joint nuisance profiles are treated as averaging", {
    grid <- expand.grid(
        A = c("a1", "a2"),
        B = c("b1", "b2"),
        C = c("c1", "c2"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    joint_key <- interaction(grid$B, grid$C, lex.order = TRUE)
    coefficient <- ifelse(grid$A == "a1", 1, -1) *
        c(0.1, 0.2, 0.3, 0.4)[as.integer(joint_key)]

    classify_role <- function(role) {
        roles <- list(
            between = if (role == "between") c("A", "B", "C") else character(),
            within = if (role == "within") c("A", "B", "C") else character(),
            levels = list(
                A = c("a1", "a2"),
                B = c("b1", "b2"),
                C = c("c1", "c2")
            )
        )
        participation <- PsyR:::.detect_factor_participation(
            coefficients = matrix(coefficient, ncol = 1L),
            grid = grid,
            roles = roles
        )[[1L]]
        list(
            participation = participation,
            classification = PsyR:::.classify_contrast_family(
                participation,
                roles
            )
        )
    }

    between_result <- classify_role("between")
    within_result <- classify_role("within")

    expect_identical(between_result$classification$family, "b")
    expect_identical(within_result$classification$family, "w")
    expect_identical(
        between_result$participation$factors$factor[
            between_result$participation$factors$participates
        ],
        "A"
    )
    expect_identical(
        within_result$participation$factors$factor[
            within_result$participation$factors$participates
        ],
        "A"
    )
    expect_identical(
        unname(between_result$participation$factors$status),
        c("varying", "averaged", "averaged")
    )
})

test_that("other-domain nuisance factors do not block joint averaging", {
    grid <- expand.grid(
        B = c("b1", "b2"),
        C = c("c1", "c2"),
        D = c("d1", "d2"),
        W = c("w1", "w2"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    roles <- list(
        between = c("B", "C", "D"),
        within = "W",
        levels = list(
            B = c("b1", "b2"),
            C = c("c1", "c2"),
            D = c("d1", "d2"),
            W = c("w1", "w2")
        )
    )
    joint_key <- interaction(grid$C, grid$D, lex.order = TRUE)
    coefficient <- as.numeric(grid$B == "b1") *
        c(0.1, 0.2, 0.3, 0.4)[as.integer(joint_key)]

    participation <- PsyR:::.detect_factor_participation(
        coefficients = matrix(coefficient, ncol = 1L),
        grid = grid,
        roles = roles
    )[[1L]]
    classification <- PsyR:::.classify_contrast_family(
        participation,
        roles
    )

    expect_identical(classification$family, "b")
    expect_identical(classification$between_factors, "B")
    expect_identical(
        unname(participation$factors$status),
        c("conditioning", "averaged", "averaged", "inactive")
    )
})

test_that("non-common same-domain joint profiles remain participating", {
    grid <- expand.grid(
        A = c("a1", "a2"),
        B = c("b1", "b2"),
        C = c("c1", "c2"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    joint_key <- interaction(grid$B, grid$C, lex.order = TRUE)
    profile_a1 <- c(0.1, 0.2, 0.3, 0.4)
    profile_a2 <- -rev(profile_a1)
    coefficient <- ifelse(
        grid$A == "a1",
        profile_a1[as.integer(joint_key)],
        profile_a2[as.integer(joint_key)]
    )
    roles <- list(
        between = c("A", "B", "C"),
        within = character(),
        levels = list(
            A = c("a1", "a2"),
            B = c("b1", "b2"),
            C = c("c1", "c2")
        )
    )

    participation <- PsyR:::.detect_factor_participation(
        coefficients = matrix(coefficient, ncol = 1L),
        grid = grid,
        roles = roles
    )[[1L]]

    expect_true(all(participation$factors$participates))
})

test_that("same-domain joint averaging supports factor-coverage checks", {
    design <- expand.grid(
        A = factor(c("a1", "a2")),
        B = factor(c("b1", "b2")),
        C = factor(c("c1", "c2")),
        replicate = seq_len(4L),
        KEEP.OUT.ATTRS = FALSE
    )
    design$id <- factor(seq_len(nrow(design)))
    design$response <- with(
        design,
        as.numeric(A) + 2 * as.numeric(B) + 3 * as.numeric(C) +
            as.numeric(replicate) / 10
    )
    model <- afex::aov_ez(
        id = "id",
        dv = "response",
        data = design,
        between = c("A", "B", "C")
    )
    cell_means <- emmeans::emmeans(model, c("A", "B", "C"))
    grid <- cell_means@grid
    joint_key <- interaction(grid$B, grid$C, lex.order = TRUE)
    coefficient <- ifelse(grid$A == "a1", 1, -1) *
        c(0.1, 0.2, 0.3, 0.4)[as.integer(joint_key)]
    contrast_table <- emmeans::contrast(
        cell_means,
        method = list("A with joint B-C averaging" = coefficient),
        adjust = "none"
    )

    report <- expect_silent(
        check_contrast_families(
            model = model,
            contrast_tables = contrast_table,
            family_list = list("b"),
            between_factors = list("A"),
            mode = "error"
        )
    )

    expect_true(report$valid)
    expect_identical(report$tables[[1L]]$between_factors, "A")
})

test_that("all grid strata contribute to domain detection", {
    fixture <- make_family_mixed_fixture()
    cell_means <- emmeans::emmeans(
        fixture$model,
        c("group", "spacing")
    )
    grid <- cell_means@grid
    coefficient <- numeric(nrow(grid))
    coefficient[grid$group == "1"] <- c(1, 2, 3)
    coefficient[grid$group == "2"] <- c(1, -1, -6)
    mixed_pattern <- emmeans::contrast(
        cell_means,
        method = list("Both-domain pattern" = coefficient),
        adjust = "none"
    )

    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = mixed_pattern,
            family_list = list("bw"),
            mode = "error"
        )
    )

    expect_identical(check$tables[[1L]]$detected_family, "bw")
})

test_that("domain detection is invariant to contrast rescaling", {
    fixture <- make_family_between_fixture()
    scaled_table <- emmeans::contrast(
        fixture$emmeans,
        method = list(
            "Tiny Extreme - Alert" = 1e-12 * c(1, -1, 0, 0)
        ),
        adjust = "none"
    )

    check <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = scaled_table,
            family_list = list("b"),
            between_factors = list("Group"),
            mode = "error"
        )
    )

    expect_identical(check$tables[[1L]]$detected_family, "b")
})

test_that("unresolved model factors make classification unverifiable", {
    roles <- list(
        between = c("A", "B"),
        within = "W",
        levels = list(
            A = c("a1", "a2"),
            B = c("b1", "b2"),
            W = c("w1", "w2")
        )
    )
    grid <- data.frame(
        A = c("a1", "a1", "a2", "a2"),
        B = c("b1", "b1", "b2", "b2"),
        W = c("w1", "w2", "w1", "w2")
    )
    participation <- PsyR:::.detect_factor_participation(
        coefficients = matrix(c(1, -1, 2, -2), ncol = 1L),
        grid = grid,
        roles = roles
    )[[1L]]
    classification <- PsyR:::.classify_contrast_family(
        participation,
        roles
    )

    expect_false(participation$metadata_complete)
    expect_setequal(participation$ambiguous_context, c("A", "B"))
    expect_identical(classification$family, "unknown")
})

test_that("incomplete joint grids are not classified from observed cells", {
    roles <- list(
        between = "A",
        within = "W",
        levels = list(
            A = c("a1", "a2"),
            W = c("w1", "w2")
        )
    )
    grid <- data.frame(
        A = c("a1", "a1", "a2"),
        W = c("w1", "w2", "w1")
    )
    participation <- PsyR:::.detect_factor_participation(
        coefficients = matrix(c(1, 1, -1), ncol = 1L),
        grid = grid,
        roles = roles
    )[[1L]]
    classification <- PsyR:::.classify_contrast_family(
        participation,
        roles
    )

    expect_true(participation$joint_grid_incomplete)
    expect_false(participation$metadata_complete)
    expect_identical(classification$family, "unknown")
})

test_that("same-domain ambiguity does not obscure a known family", {
    grid <- expand.grid(
        A = c("a1", "a2"),
        B = c("b1", "b2", "b3"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    roles <- list(
        between = c("A", "B"),
        within = character(),
        levels = list(
            A = c("a1", "a2"),
            B = c("b1", "b2", "b3")
        )
    )
    coefficient <- ifelse(grid$A == "a1", 1, -1) *
        c(0.5, 0.5, 0)[match(grid$B, c("b1", "b2", "b3"))]
    participation <- PsyR:::.detect_factor_participation(
        coefficients = matrix(coefficient, ncol = 1L),
        grid = grid,
        roles = roles
    )[[1L]]
    classification <- PsyR:::.classify_contrast_family(
        participation,
        roles
    )

    expect_true("B" %in% participation$ambiguous_context)
    expect_identical(classification$family, "b")
    expect_identical(classification$scope, "unknown")
})

test_that("heterogeneous tables are reported instead of majority classified", {
    fixture <- make_family_mixed_fixture()
    cell_means <- emmeans::emmeans(
        fixture$model,
        c("group", "spacing")
    )
    grid <- cell_means@grid
    marginal_b <- ifelse(
        grid$group == "1",
        1 / 3,
        ifelse(grid$group == "2", -1 / 3, 0)
    )
    interaction_bw <- ifelse(
        grid$group == "1" & grid$spacing == "TWENTY",
        1,
        ifelse(
            grid$group == "2" & grid$spacing == "TWENTY",
            -1,
            ifelse(
                grid$group == "1" & grid$spacing == "FORTY",
                -1,
                ifelse(
                    grid$group == "2" & grid$spacing == "FORTY",
                    1,
                    0
                )
            )
        )
    )
    mixed_table <- emmeans::contrast(
        cell_means,
        method = list(
            "Marginal B" = marginal_b,
            "BW interaction" = interaction_bw
        ),
        adjust = "none"
    )

    expect_error(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = mixed_table,
            family_list = list("bw"),
            mode = "error"
        ),
        regexp = "mixed families|separate",
        ignore.case = TRUE
    )
    expect_error(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = mixed_table
        ),
        regexp = "more than one family|separate",
        ignore.case = TRUE
    )
    expect_error(
        psyci(
            model = fixture$model,
            contrast_tables = mixed_table,
            method = "bf",
            family_check = "none"
        ),
        regexp = "one resolved family|split mixed-family",
        ignore.case = TRUE
    )
})

test_that("known mismatches are retained when another contrast is unknown", {
    fixture <- make_family_between_fixture()
    partly_unknown <- emmeans::contrast(
        fixture$emmeans,
        method = list(
            "Extreme - Alert" = c(1, -1, 0, 0),
            "Uninformative" = c(0, 0, 0, 0)
        ),
        adjust = "none"
    )

    expect_error(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = partly_unknown,
            family_list = list("w"),
            mode = "error"
        ),
        regexp = "declared.*w|detected|Extreme - Alert",
        ignore.case = TRUE
    )
})

test_that("factor-role arguments are checked independently of family labels", {
    fixture <- make_family_mixed_fixture()

    expect_error(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = fixture$between,
            family_list = list("b"),
            between_factors = list("spacing"),
            mode = "error"
        ),
        regexp = "within-subject",
        ignore.case = TRUE
    )
})

test_that("supplied factors must cover the factors detected in a table", {
    data("drive_sleepy", package = "PsyR", envir = environment())
    drive_sleepy$sleep_deprivation <- factor(
        drive_sleepy$sleep_deprivation
    )
    drive_sleepy$nvh <- factor(drive_sleepy$nvh)
    model <- afex::aov_ez(
        id = "subID",
        dv = "performance",
        data = drive_sleepy,
        between = c("sleep_deprivation", "nvh")
    )
    sleep_emmeans <- emmeans::emmeans(model, "sleep_deprivation")
    sleep_table <- emmeans::contrast(
        sleep_emmeans,
        method = list("0 - 12 hours" = c(1, -1)),
        adjust = "none"
    )

    expect_error(
        check_contrast_families(
            model = model,
            contrast_tables = sleep_table,
            family_list = list("b"),
            between_factors = list("nvh"),
            mode = "error"
        ),
        regexp = "factor|sleep_deprivation|missing",
        ignore.case = TRUE
    )

    expect_silent(
        check_contrast_families(
            model = model,
            contrast_tables = sleep_table,
            family_list = list("b"),
            between_factors = list("sleep_deprivation", "nvh"),
            mode = "error"
        )
    )
})

test_that("none mode normalizes factor arguments without signalling", {
    fixture <- make_family_between_fixture()
    report <- expect_silent(
        check_contrast_families(
            model = fixture$model,
            contrast_tables = fixture$table,
            family_list = list("b"),
            between_factors = list("Group", "Group", NA),
            mode = "none"
        )
    )

    expect_identical(report$factor_arguments$between, "Group")
})

test_that("single-contrast Bonferroni intervals count one contrast", {
    fixture <- make_family_between_fixture()
    result <- psyci(
        model = fixture$model,
        contrast_tables = fixture$table,
        method = "bf",
        family_list = list("b"),
        between_factors = list("Group")
    )
    error_df <- unique(as.data.frame(summary(fixture$table))$df)

    expect_equal(
        unique(result[[1L]]$cc),
        cc_bonf_t(n_k = 1, v_e = error_df, alpha = 0.05),
        tolerance = 1e-12
    )
})

test_that("psyci infers post-hoc family and factor arguments", {
    fixture <- make_family_between_fixture()
    inferred <- psyci(
        model = fixture$model,
        contrast_tables = fixture$table,
        method = "ph"
    )
    explicit <- psyci(
        model = fixture$model,
        contrast_tables = fixture$table,
        method = "ph",
        family_list = list("b"),
        between_factors = list("Group")
    )

    expect_equal(inferred, explicit, tolerance = 1e-12)
    expect_identical(names(inferred), "b")
    expect_match(
        paste(attr(inferred[[1L]], "mesg"), collapse = " "),
        "between subject factor\\(s\\) are: Group"
    )
})

test_that("automatic family inference resolves mixed-design tables", {
    fixture <- make_family_mixed_fixture()
    tables <- list(fixture$between, fixture$within, fixture$interaction)
    report <- check_contrast_families(
        model = fixture$model,
        contrast_tables = tables
    )
    inferred <- psyci(
        model = fixture$model,
        contrast_tables = tables,
        method = "ph"
    )

    expect_identical(report$families, c("b", "w", "bw"))
    expect_true(report$inferred)
    expect_true(all(is.na(report$summary$declared_family)))
    expect_identical(report$summary$resolved_family, c("b", "w", "bw"))
    expect_identical(names(inferred), c("b", "w", "bw"))
})

test_that("post-hoc dimensions are calculated separately for each table", {
    design <- expand.grid(
        A = factor(c("a1", "a2")),
        B = factor(c("b1", "b2", "b3")),
        replicate = seq_len(5L),
        KEEP.OUT.ATTRS = FALSE
    )
    design$id <- factor(seq_len(nrow(design)))
    design$response <- with(
        design,
        as.numeric(A) + as.numeric(B) + as.numeric(replicate) / 10
    )
    model <- afex::aov_ez(
        id = "id",
        dv = "response",
        data = design,
        between = c("A", "B")
    )
    table_a <- emmeans::contrast(
        emmeans::emmeans(model, "A"),
        list("A1 - A2" = c(1, -1))
    )
    table_b <- emmeans::contrast(
        emmeans::emmeans(model, "B"),
        list("B1 - B2" = c(1, -1, 0))
    )
    result <- psyci(
        model = model,
        contrast_tables = list(table_a, table_b),
        method = "ph"
    )
    error_df <- unique(as.data.frame(summary(table_a))$df)

    expect_equal(
        unique(result[[1L]]$cc),
        cc_ph_b(v_b = 1, v_e = error_df, alpha = 0.05),
        tolerance = 1e-12
    )
    expect_equal(
        unique(result[[2L]]$cc),
        cc_ph_b(v_b = 2, v_e = error_df, alpha = 0.05),
        tolerance = 1e-12
    )
    expect_false(isTRUE(all.equal(
        unique(result[[1L]]$cc),
        unique(result[[2L]]$cc)
    )))
})

test_that("psyci uses normalized factor arguments for post-hoc dimensions", {
    fixture <- make_family_between_fixture()
    standard <- psyci(
        model = fixture$model,
        contrast_tables = fixture$table,
        method = "ph",
        family_list = list("b"),
        between_factors = list("Group")
    )
    duplicated <- psyci(
        model = fixture$model,
        contrast_tables = fixture$table,
        method = "ph",
        family_list = list("b"),
        between_factors = list("Group", "Group", NA)
    )

    expect_equal(duplicated, standard, tolerance = 1e-12)
})

test_that("independent post-hoc intervals require a resolved factor dimension", {
    design <- expand.grid(
        A = factor(c("a1", "a2")),
        B = factor(c("b1", "b2", "b3")),
        replicate = seq_len(4L),
        KEEP.OUT.ATTRS = FALSE
    )
    design$id <- factor(seq_len(nrow(design)))
    design$response <- with(
        design,
        as.numeric(A) + 2 * as.numeric(B) + as.numeric(replicate) / 10
    )
    model <- afex::aov_ez(
        id = "id",
        dv = "response",
        data = design,
        between = c("A", "B")
    )
    cell_means <- emmeans::emmeans(model, c("A", "B"))
    grid <- cell_means@grid
    coefficient <- ifelse(grid$A == "a1", 1, -1) *
        c(b1 = 0.5, b2 = 0.5, b3 = 0)[as.character(grid$B)]
    contrast_table <- emmeans::contrast(
        cell_means,
        method = list("A over part of B" = coefficient),
        adjust = "none"
    )
    check <- check_contrast_families(
        model = model,
        contrast_tables = contrast_table,
        family_list = list("b"),
        between_factors = list("A"),
        mode = "error"
    )
    args <- list(
        model = model,
        contrast_tables = contrast_table,
        method = "ph",
        family_list = list("b"),
        between_factors = list("A")
    )

    expect_identical(check$tables[[1L]]$detected_family, "b")
    expect_identical(check$tables[[1L]]$scope, "unknown")
    expect_error(
        do.call(psyci, args),
        regexp = "participating-factor dimension|post-hoc critical",
        ignore.case = TRUE
    )
    warned <- expect_warning(
        do.call(psyci, c(args, list(family_check = "warn"))),
        regexp = "participating-factor dimension|post-hoc critical",
        ignore.case = TRUE
    )
    silent <- expect_silent(
        do.call(psyci, c(args, list(family_check = "none")))
    )
    expect_type(warned, "list")
    expect_type(silent, "list")
})

test_that("psyci documents strict handling of summary_emm inputs", {
    fixture <- make_family_between_fixture()
    summary_table <- summary(fixture$table)
    args <- list(
        model = fixture$model,
        contrast_tables = summary_table,
        method = "ind",
        family_list = list("b"),
        between_factors = list("Group")
    )

    for (mode in c("error", "warn", "none")) {
        expect_error(
            do.call(psyci, c(args, list(family_check = mode))),
            regexp = "original emmGrid|genuine contrasts|orthogonality",
            ignore.case = TRUE
        )
    }
    expect_error(
        psyci(
            model = fixture$model,
            contrast_tables = summary_table,
            method = "ind"
        ),
        regexp = "original emmGrid|genuine contrasts|orthogonality",
        ignore.case = TRUE
    )
})
