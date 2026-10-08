# Independent reference implementation of the calculations in the original
# Delphi/Pascal Psy program. This file deliberately does not call PsyR's
# confidence-interval or critical-value functions.

legacy_psy_design <- function(data, id, dv, between = character(),
                              within = character()) {
    make_grid <- function(columns) {
        if (length(columns) == 0L) {
            return(data.frame(.cell = "1"))
        }
        values <- lapply(data[columns], function(x) {
            if (is.factor(x)) levels(x) else unique(x)
        })
        names(values) <- columns
        expand.grid(values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    }
    make_key <- function(x) {
        if (ncol(x) == 0L) rep("1", nrow(x)) else do.call(paste, c(x, sep = "\r"))
    }

    group_grid <- make_grid(between)
    repeat_grid <- make_grid(within)
    group_key <- if (length(between)) make_key(data[between]) else rep("1", nrow(data))
    repeat_key <- if (length(within)) make_key(data[within]) else rep("1", nrow(data))
    group_levels <- if (length(between)) make_key(group_grid) else "1"
    repeat_levels <- if (length(within)) make_key(repeat_grid) else "1"
    subjects <- unique(as.character(data[[id]]))
    wide <- matrix(
        NA_real_, nrow = length(subjects), ncol = length(repeat_levels),
        dimnames = list(subjects, repeat_levels)
    )
    groups <- integer(length(subjects))

    for (i in seq_along(subjects)) {
        rows <- which(as.character(data[[id]]) == subjects[[i]])
        subject_groups <- unique(match(group_key[rows], group_levels))
        if (length(subject_groups) != 1L || is.na(subject_groups)) {
            stop("Each subject must belong to exactly one between-subject group.")
        }
        groups[[i]] <- subject_groups
        repeat_cells <- match(repeat_key[rows], repeat_levels)
        if (anyDuplicated(repeat_cells) || anyNA(repeat_cells)) {
            stop("Each subject must have one observation in every repeated cell.")
        }
        wide[i, repeat_cells] <- data[[dv]][rows]
    }
    if (anyNA(wide)) {
        stop("The legacy Psy calculation requires complete repeated measures.")
    }

    list(
        data = wide,
        group = groups,
        group_grid = group_grid,
        repeat_grid = repeat_grid,
        n_group = tabulate(groups, nbins = nrow(group_grid)),
        dfe = nrow(wide) - nrow(group_grid)
    )
}

# GCR.pas: Ksmn and RoyExact. BetaRoy is the unregularised incomplete beta.
.legacy_beta_roy <- function(x, m, n) {
    stats::pbeta(x, m + 1, n + 1) * beta(m + 1, n + 1)
}

.legacy_ksmn <- function(s, m, n) {
    i <- seq_len(round(s))
    exp(
        s / 2 * log(pi) +
            sum(lgamma((2 * m + 2 * n + s + i + 2) / 2)) -
            sum(lgamma((2 * m + i + 1) / 2)) -
            sum(lgamma((2 * n + i + 1) / 2)) -
            sum(lgamma(i / 2))
    )
}

.legacy_roy_exact <- function(s, m, n, x) {
    b <- b0 <- numeric(7)
    k <- pr <- numeric(4)
    k[2] <- .legacy_ksmn(2, m, n)
    b[1] <- .legacy_beta_roy(x, 2 * m + 1, 2 * n + 1)
    b0[1] <- x^(m + 1) * (1 - x)^(n + 1)
    b[2] <- .legacy_beta_roy(x, m, n)
    pr[2] <- k[2] / (m + n + 2) * (2 * b[1] - b0[1] * b[2])
    if (s == 2) return(pr[2])

    k[3] <- .legacy_ksmn(3, m, n)
    b0[2] <- x^(m + 2) * (1 - x)^(n + 1)
    b[3] <- .legacy_beta_roy(x, 2 * m + 3, 2 * n + 1)
    b[4] <- .legacy_beta_roy(x, 2 * m + 2, 2 * n + 1)
    b[5] <- .legacy_beta_roy(x, m + 1, n)
    pr[3] <- k[3] / (m + n + 3) *
        (2 * b[3] * b[2] - 2 * b[4] * b[5] - b0[2] * pr[2] / k[2])
    if (s == 3) return(pr[3])

    k[4] <- .legacy_ksmn(4, m, n)
    b0[3] <- x^(m + 3) * (1 - x)^(n + 1)
    b[6] <- .legacy_beta_roy(x, 2 * m + 5, 2 * n + 1)
    b[7] <- .legacy_beta_roy(x, 2 * m + 4, 2 * n + 1)
    work <- -b0[3] * pr[3] / k[3] + 2 * b[6] * pr[2] / k[2] +
        2 * b[3] / (m + n + 3) * (2 * b[3] - b0[2] * b[5]) -
        2 * b[7] / (m + n + 3) *
            (-b0[2] * b[2] + (m + 2) * pr[2] / k[2] + 2 * b[4])
    k[4] / (m + n + 4) * work
}

# GCR.pas: PillaiApproxOriginal. This path is used by Psy when m < 0.
.legacy_pillai_approx <- function(s, m, n, x) {
    ln_c <- function(s0) {
        i <- seq_len(round(s0))
        s0 / 2 * log(pi) +
            sum(lgamma((2 * m + 2 * n + s0 + i + 2) / 2)) -
            sum(lgamma((2 * m + i + 1) / 2) +
                lgamma((2 * n + i + 1) / 2) + lgamma(i / 2))
    }
    si <- round(s)
    hsx <- if (si %% 2L) stats::pbeta(x, m + 1, n + 1) else 1
    kms <- numeric(si)
    ratio <- exp(ln_c(s) - ln_c(s - 1))
    for (i in seq_len(si - 1L)) {
        product <- if (i == 1L) 1 else prod(
            (2 * m + s - seq_len(i - 1L) + 1) /
                (2 * m + 2 * n + 2 * s - seq_len(i - 1L) + 1)
        )
        previous <- if (i == 1L) 0 else kms[i - 1L]
        kms[i] <- (ratio * choose(s - 1, i - 1) * product -
            (m + s - i + 1) * previous) / (m + n + s - i + 1)
    }
    terms <- (-1)^seq_len(si - 1L) * kms[seq_len(si - 1L)] *
        x^(s - seq_len(si - 1L))
    hsx + exp(m * log(x) + (n + 1) * log1p(-x)) * sum(terms)
}

legacy_gcr_critical <- function(alpha, s, m, n) {
    cdf <- if (m < 0) {
        function(x) .legacy_pillai_approx(s, m, n, x)
    } else if (s <= 4 && n <= 600) {
        function(x) .legacy_roy_exact(s, m, n, x)
    } else {
        function(x) .legacy_pillai_approx(s, m, n, x)
    }
    stats::uniroot(function(x) 1 - cdf(x) - alpha, c(0.001, 0.999),
        tol = 1e-10)$root
}

# SMR.pas: Gamma_ab, Gfnr and the m = 2 Davis expression. The vignette
# analysis requests p = q = 2, so no unused Davis cases are ported here.
.legacy_gamma_ab <- function(alpha_minus_one, rate, x) {
    alpha <- alpha_minus_one + 1
    rate^(-alpha) * gamma(alpha) * stats::pgamma(rate * x, shape = alpha)
}

.legacy_gfnr <- function(r, a, q, b, l) {
    if (r == 0L) {
        ifelse(l == 0, 0, exp(-a * l + log(l) * (a * q + b)))
    } else if (r == 1L) {
        .legacy_gamma_ab(a * q + b, a, l)
    } else {
        stop("Only the Gfnr cases required by Davis(m = 2) are implemented.")
    }
}

.legacy_davis_2 <- function(q, l) {
    exp(-lgamma(q - 1)) * (
        .legacy_gfnr(1L, 1, q, -2, l) -
            0.5 * .legacy_gfnr(0L, 0.5, q, -0.5, l) *
                .legacy_gfnr(1L, 0.5, q, -1.5, l)
    )
}

legacy_smr_critical <- function(p, q, dfe, alpha) {
    if (min(p, q) != 2L) {
        stop("The source-faithful SMR port currently covers p = 2 or q = 2.")
    }
    m <- min(p, q)
    q0 <- max(p, q)
    davis <- Vectorize(function(l) .legacy_davis_2(q0, l))
    cdf <- function(c) {
        stats::integrate(
            function(x) davis(c * x / dfe) * stats::dchisq(x, dfe),
            0, Inf, subdivisions = 500L, rel.tol = 1e-9
        )$value
    }
    lower <- stats::uniroot(function(x) 1 - davis(x) - alpha,
        c(1, 2000), tol = 1e-8)$root
    stats::uniroot(function(x) 1 - cdf(x) - alpha,
        c(lower, 2000), tol = 1e-8)$root
}

legacy_psy_critical <- function(method, alpha, family, dfe, n_k = 1L,
                                groups, repeats, smr = NULL) {
    if (method == "ind") return(sqrt(stats::qf(1 - alpha, 1, dfe)))
    if (method == "bf") return(sqrt(stats::qf(1 - alpha / n_k, 1, dfe)))
    if (method == "smr") {
        return(sqrt(legacy_smr_critical(smr$p, smr$q, dfe, alpha)))
    }
    if (method != "ph") stop("Unknown legacy Psy method.")
    if (family == "b") {
        return(sqrt((groups - 1) * stats::qf(1 - alpha, groups - 1, dfe)))
    }
    if (family == "w") {
        return(sqrt((repeats - 1) * dfe / (dfe - repeats + 2) *
            stats::qf(1 - alpha, repeats - 1, dfe - repeats + 2)))
    }
    s <- min(groups - 1, repeats - 1)
    if (s == 1L) {
        if (groups == 2L) {
            return(legacy_psy_critical("ph", alpha, "w", dfe,
                groups = groups, repeats = repeats))
        }
        return(legacy_psy_critical("ph", alpha, "b", dfe,
            groups = groups, repeats = repeats))
    }
    theta <- legacy_gcr_critical(
        alpha, s, (abs(groups - repeats) - 1) / 2,
        (dfe - repeats) / 2
    )
    sqrt(dfe * theta / (1 - theta))
}

# PsyFile.pas: NormaliseBCoefficients, GetSSErrorValue, CalculateW,
# SSOfBcontrasts, CalculateFinalSampleValue and WriteConfidenceIntervals.
# Each item in contrasts is list(label, family, a, w). Null a/w means A0/W0.
legacy_psy_analyse <- function(design, contrasts, method, alpha = 0.05,
                               n_k = NULL, smr = NULL, critical = NULL) {
    group_count <- nrow(design$group_grid)
    repeat_count <- nrow(design$repeat_grid)
    group_means <- vapply(seq_len(group_count), function(g) {
        colMeans(design$data[design$group == g, , drop = FALSE])
    }, numeric(repeat_count))
    if (group_count == 1L) group_means <- matrix(group_means, ncol = 1L)

    rows <- lapply(seq_along(contrasts), function(i) {
        item <- contrasts[[i]]
        a <- if (is.null(item$a)) rep(1, group_count) else item$a
        w <- if (is.null(item$w)) rep(1, repeat_count) else item$w
        if (length(a) != group_count || length(w) != repeat_count) {
            stop("Contrast coefficients do not match the legacy design dimensions.")
        }
        normalised_scores <- as.vector(design$data %*% (w / sqrt(sum(w^2))))
        fitted <- ave(normalised_scores, design$group, FUN = mean)
        mse <- sum((normalised_scores - fitted)^2) / design$dfe
        estimate <- sum(a * as.vector(crossprod(w, group_means)))
        se <- sqrt(sum(a^2 / design$n_group) * mse * sum(w^2))
        if (is.null(item$a)) {
            estimate <- estimate / group_count
            se <- se / group_count
        }
        if (is.null(item$w)) {
            estimate <- estimate / repeat_count
            se <- se / repeat_count
        }
        cc <- if (!is.null(critical)) {
            if (length(critical) == 1L) critical else critical[[i]]
        } else {
            this_n <- if (is.null(n_k)) 1L else if (length(n_k) == 1L) n_k else n_k[[i]]
            legacy_psy_critical(
                method, alpha, item$family, design$dfe, this_n,
                group_count, repeat_count, smr
            )
        }
        data.frame(
            contrast = item$label, family = item$family,
            estimate = estimate, SE = se, df = design$dfe, cc = cc,
            lower = estimate - cc * se, upper = estimate + cc * se,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

legacy_psy_cross <- function(a = list(), w = list(), include_b = TRUE,
                             include_w = TRUE, include_bw = TRUE) {
    result <- list()
    if (include_b) {
        result <- c(result, lapply(names(a), function(an) {
            list(label = an, family = "b", a = a[[an]], w = NULL)
        }))
    }
    if (include_w) {
        result <- c(result, lapply(names(w), function(wn) {
            list(label = wn, family = "w", a = NULL, w = w[[wn]])
        }))
    }
    if (include_bw) {
        for (wn in names(w)) for (an in names(a)) {
            result[[length(result) + 1L]] <- list(
                label = paste(an, wn, sep = " x "), family = "bw",
                a = a[[an]], w = w[[wn]]
            )
        }
    }
    result
}
