# Internal helpers for the SMR reference distribution.

.fit_smr_reference_distribution <- function(p, q, n, n_sim, seed) {
    restore_seed <- !is.null(seed)
    if (restore_seed) {
        had_random_seed <- exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
        )
        if (had_random_seed) {
            previous_random_seed <- get(
                ".Random.seed",
                envir = .GlobalEnv,
                inherits = FALSE
            )
        }
        on.exit({
            if (had_random_seed) {
                assign(
                    ".Random.seed",
                    previous_random_seed,
                    envir = .GlobalEnv
                )
            } else if (exists(
                ".Random.seed",
                envir = .GlobalEnv,
                inherits = FALSE
            )) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        }, add = TRUE)
    }

    moments <- wishartlr_sample(
        p = p,
        q = q,
        n_sim = n_sim,
        seed = seed
    )

    k <- (moments$mu2 + moments$mu1^2) * (n - 2) /
        ((n - 4) * moments$mu1^2)
    u <- (moments$mu3 + 3 * moments$mu2 * moments$mu1 + moments$mu1^3) *
        (n - 2)^2 /
        (moments$mu1^3 * (n - 4) * (n - 6))

    f_pi <- 2 * (k + 3 * u - 4 * k^2) / (k + u - 2 * k^2)
    f_gamma <- 2 * (f_pi - 2) / ((f_pi - 2) * (k - 1) - 2 * k)
    f_omega <- moments$mu1 * n * (f_pi - 2) / ((n - 2) * f_pi)

    parameters <- c(
        f_gamma = f_gamma,
        f_pi = f_pi,
        f_omega = f_omega
    )
    if (any(!is.finite(parameters)) || any(parameters <= 0)) {
        stop(
            paste(
                "The simulated SMR reference distribution produced invalid",
                "scaled-F parameters. Try increasing `n_sim` or using a",
                "different `seed`."
            ),
            call. = FALSE
        )
    }

    as.list(parameters)
}

.smr_reference_quantile <- function(alpha, distribution) {
    stats::qf(
        1 - alpha,
        df1 = distribution$f_gamma,
        df2 = distribution$f_pi
    ) * distribution$f_omega
}

.smr_reference_p_value <- function(statistic, distribution) {
    stats::pf(
        statistic / distribution$f_omega,
        df1 = distribution$f_gamma,
        df2 = distribution$f_pi,
        lower.tail = FALSE
    )
}
