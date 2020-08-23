is_outlier <- function(x, up_prob = 0.995) {
  return(x > 5 * quantile(x, up_prob))
}

empirical_prob_zero <- function(x, rmoutliers = T) {
  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
  }
  invisible(1 - length(which(x > 0)) / length(x))

}

fit_poi_with_scalefactors <- function(x, s, rmoutliers = T) {
  if (length(x) != length(s)) {
    stop("unequal length of x and s.")
  }

  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
    s <- s[!outliers]
  }
  sumcnts <- sum(x)
  sumfactors <- sum(s)
  if (sumfactors < 1) {
    stop("all of s are zeros ?")
  }

  estimate <- sumcnts / sumfactors
  sds <- sqrt(estimate / sumfactors)
  names(estimate) <- names(sds) <- "lambda"
  vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("lambda",
    "lambda"))
  loglik <- sum(dpois(x, lambda = s * estimate, log = TRUE))
  invisible(structure(list(estimate = estimate, sd = sds,
    vcov = vc, n = length(x), loglik = loglik), class = "fitdistr"))
}

prob_zero_poi_scale <- function(x, s, rmoutliers = T) {
  thefit <- fit_poi_with_scalefactors(x, s, rmoutliers)
  scaled_mean <- thefit$estimate
  p0 <- exp(-scaled_mean * s)
  invisible(list(fitpoiscl = thefit,
    mnp0 = mean(p0),
    mdnp0 = median(p0)))
}

prob_zero_poi <- function(x, rmoutliers = T) {
  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
  }

  thefit <- MASS::fitdistr(x, densfun = "poisson")
  p0 <- exp(-thefit$estimate)
  invisible(list(fitpoi = thefit, p0 = p0))
}

prob_zero_nb <- function(x, rmoutliers = T) {
  if (rmoutliers) {
    outliers <- is_outlier(x)
    x <- x[!outliers]
  }
  nbfit <- MASS::fitdistr(x, densfun = "negative binomial")
  r <- nbfit$size
  mu <- nbfit$mu
  prob <- r / (r + mu)
  thevar <- mu + mu^2 / r
}
