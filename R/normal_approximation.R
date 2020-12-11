#' Returns corresponding normalization method
#'
#' Given a string, returns the corresponding normalization function
#' @param name string
#'
#' @return function
#' @export
#'
#' @examples
#' get_normalization_method("alfers_dinges")
get_normalization_method <- function(name) {
  return(switch(name,
    "alfers_dinges" = alfers_dinges_approx,
    "wise" = wise_approx,
    "peizer_pratt" = peizer_pratt_approx
  ))
}



#' Normal approximation for Beta distribution
#'
#' Takes in a vector Y containing data from a beta distribution and returns data that is approximately normal.
#' Uses approximation found in "A Normal Approximation for Beta and Gamma Tail Probabilities"
#' (Alfers and Dinges)
#' @param Y vector
#'
#' @return vector
#' @export
#'
#' @examples
#' alfers_dinges_approx(rbeta(1000, 1, 10))
alfers_dinges_approx <- function(Y) {
  A <- function(alpha, p) {
    q <- 1 - p
    beta <- 1 - alpha

    pi <- function(alpha, p) {
      h <- function(x) {
        return(1 / x + (1 - x) * log(1 - x) / x^2 - 0.5)
      }

      return(2 * p * h(1 - (beta / q)) + 2 * q * h(1 - (alpha / p)))
    }

    return(((alpha - p) / sqrt(p * q)) * sqrt(1 + pi(alpha, p)))
  }

  alpha <- mean(Y)
  return(A(alpha, Y))
}



#' Normal approximation for Beta distribution
#'
#' Takes in a vector Y containing data from a beta distribution and returns data that is approximately normal.
#' Uses approximation found in "A Normal Approximation for Binomial, F, Beta, and Other Common, Related Tail Probabilities, I"
#' (Peizer and Pratt)
#' @param Y vector
#'
#' @return vector
#' @export
#'
#' @examples
#' peizer_pratt_approx(rbeta(1000, 1, 10))
peizer_pratt_approx <- function(Y) {
  estBetaParams <- function(mu, variance) {
    alpha <- ((1 - mu) / variance - 1 / mu) * mu^2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }

  d1 <- function() {
    return(S + 1 / 6 - (n + 1 / 3) * p)
  }

  d2 <- function() {
    return(d1() + 0.02 * ((q / (S + 0.5)) - (p / (T + 0.5)) + ((q - 0.5) / (n + 1))))
  }

  g <- function(x) {
    return((1 - x)^(-2) * (1 - x^2 + 2 * x * log(x)))
  }

  params <- estBetaParams(mean(Y), var(Y))

  S <- params$beta - 0.5
  T <- params$alpha - 0.5
  n <- params$alpha + params$beta - 1
  p <- 1 - Y
  q <- Y

  return(d2() * ((1 + q * g(S / (n * p)) + p * g(T / (n * q))) / ((n + 1 / 6) * p * q))^(1 / 2))
}



#' Normal approximation for Beta distribution
#'
#' Takes in a vector Y containing data from a beta distribution and returns data that is approximately normal.
#' Uses Wise's transformation: (-logY)^(1/3)
#' @param Y vector
#'
#' @return vector
#' @export
#'
#' @examples
#' wise_approx(rbeta(1000, 1, 10))
wise_approx <- function(Y) {
  return((-log(as.numeric(Y)))^(1 / 3))
}
