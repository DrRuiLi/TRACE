distinct_norm_from_random_backgroud <- function(
    x_norm,x_random){


  bg_kde <- density(x_random, n = 1e5)

  f_bg <- function(x) {
    approx(bg_kde$x, bg_kde$y, xout = x,
           rule = 2, ties = mean)$y
  }

  fit_bg_norm_mixture <- function(
    x,
    f_bg,
    max_iter = 1e3,
    tol = 1e-6
  ) {
    n <- length(x)

    pi  <- 0.2
    mu  <- mean(x)
    sd  <- sd(x)

    loglik_old <- -Inf

    for (iter in seq_len(max_iter)) {

      bg_d <- f_bg(x)
      bg_d[bg_d <= 0] <- min(bg_d[bg_d > 0]) * 1e-3

      norm_d <- dnorm(x, mu, sd)

      w <- pi * norm_d / ((1 - pi) * bg_d + pi * norm_d)

      pi <- mean(w)
      mu <- sum(w * x) / sum(w)
      sd <- sqrt(sum(w * (x - mu)^2) / sum(w))

      loglik <- sum(log((1 - pi) * bg_d + pi * norm_d))
      if (abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik
    }

    list(
      pi = pi,
      mu = mu,
      sd = sd,
      posterior = w,
      iter = iter
    )
  }

  fit <- fit_bg_norm_mixture(x_norm, f_bg)
  return(fit)

}



extract_formula_CN <- function(x){
  x <- stringr::str_extract_all(x, "[CN]\\d+") |>
    vapply(paste0, collapse = "", FUN.VALUE = character(1))
  paste0(x, ifelse(grepl("N", x), "", "N0"))
}



