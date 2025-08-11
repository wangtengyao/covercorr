#' Coverage-based Dependence Measure with Optional Visualisation
#'
#' Computes the coverage correlation coefficient between input \code{x} and \code{y}, as introduced in \href{https://arxiv.org/abs/2508.06402}{the arXiv preprint}. This coefficient measures the dependence between two random variables or vectors.
#'
#' The procedure is as follows:
#' \enumerate{
#'   \item Calculate the rank transformations \eqn{(r_x, r_y)} of the inputs \code{x} and \code{y}.
#'   \item Construct small cubes (in 2D, squares) of volume \eqn{n^{-1}} centered at each rank-transformed point.
#'   \item Compute the total area of the union of these cubes, intersected with \eqn{[0,1]^d} where \eqn{d = d_x + d_y}.
#' }
#' The coverage correlation coefficient is then calculated based on this union area.
#'
#' For more details, please refer to the original paper: \href{https://arxiv.org/abs/2508.06402}{the arXiv preprint}.
#'
#' @param x Numeric vector or matrix.
#' @param y Numeric vector or matrix with the same number of rows as \code{x}.
#' @param visualise Logical; if \code{TRUE}, displays a scatter plot of the
#'   rank-transformed points with overlaid rectangles to illustrate the coverage
#'   calculation. The default is \code{FALSE} (no plot). If set to \code{TRUE}
#'   but either \code{x} or \code{y} has more than one column, a warning is
#'   issued and \code{visualise} is reset to \code{FALSE}.
#' @param method Character string specifying the computation method. Options are \code{"auto"}, \code{"exact"}, or \code{"approx"}. See Details.
#' @param M Integer; Number of Monte Carlo integration sample points (used when \code{method = "approx"}). Optional.
#' @param na.rm Logical; if \code{TRUE}, remove \code{NA} values before computation.
#'
#' @details
#' The \code{method} argument controls how the computation is performed:
#' \itemize{
#'   \item \code{"exact"}: Computes the exact value.
#'   \item \code{"approx"}: Uses a Monte Carlo approximation with \code{M} sample points.
#'   \item \code{"auto"}: Automatically selects a method based on the total number of columns in \code{x} and \code{y}: if more than 6, \code{"approx"} is used (with \code{M = nrow(x)^{1.5}} if \code{M} is not provided); otherwise, \code{"exact"} is used.
#' }
#'
#' @return A list with four elements:
#' \itemize{
#'   \item \code{stat} – The numeric value of the coverage correlation coefficient.
#'   \item \code{pval} – The p-value, calculated using the exact variance under the null hypothesis of independence between \code{x} and \code{y}.
#'   \item \code{method} – A character string indicating the computation method used.
#'   \item \code{mc_se} – A numeric value. If method "approx" was used \code{mc_se} is the standard error of the Monte Carlo approximation, otherwise it is 0.
#' }
#'
#' @importFrom stats runif
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- runif(n)
#' y <- sin(3*x) + runif(n) * 0.01
#' coverage_correlation(x, y, visualise = TRUE)
#'
#' @export

coverage_correlation <- function(x, y, visualise = FALSE,
                                 method = c('auto', 'exact', 'approx'),
                                 M = NULL,
                                 na.rm = TRUE){

  # if inputs are not in proper matrix format change if possible
  if(!is.matrix(x)) {
    x = as.matrix(x)
  }
  if(!is.matrix(y)) {
    y = as.matrix(y)
  }

  if (nrow(y) != nrow(x)) stop("Number of rows of x and y should be equal.")


  if (na.rm == TRUE) {
    # NAs are removed here:
    ok = complete.cases(x, y)
    x = as.matrix(x[ok,])
    y = as.matrix(y[ok,])
  }

  n <- nrow(x)
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")

  dx <- ncol(x)
  dy <- ncol(y)

  # Visualise is only supported in 2D (first two coords) — warn and disable otherwise
  if (isTRUE(visualise) && (dx > 1 || dy > 1)) {
    warning("visualise=TRUE requires univariate x and y; turning visualise=FALSE.")
    visualise <- FALSE
  }

  d <- dx + dy

  method <- match.arg(method)
  if (method == 'auto') method <- ifelse(d <= 6, 'exact', 'approx')

  # MK ranks
  u <- matrix(runif(n * dx), n)
  v <- matrix(runif(n * dy), n)
  x_rank <- MK_rank(x, u)
  y_rank <- MK_rank(y, v)
  eps <- n^(-1/d) / 2
  zmin <- cbind(x_rank - eps, y_rank - eps)
  zmax <- cbind(x_rank + eps, y_rank + eps)

  # Wrap around [0,1]^d (split rectangles that cross boundaries); in C
  ret <- split_rectangles(zmin, zmax)

  zmin_s <- ret$zmin
  zmax_s <- ret$zmax

  # --- optional visualisation ---#
  if (visualise){
    xmin <- zmin_s[, 1]; ymin <- zmin_s[, 2]
    xmax <- zmax_s[,1]; ymax <- zmax_s[, 2]
    plot(x_rank, y_rank, pch = 20, cex = 0.3, asp = 1)
    plot_rectangles(xmin, xmax, ymin, ymax, add = TRUE)
  }

  # --- Covered volume --- #
  if (method == 'exact'){
    total_volume <- covered_volume_partitioned(zmin_s, zmax_s)
  } else {
    if (is.null(M)) {
      M <- as.integer(ceiling(n^(1.5)))
    }
    cv <- covered_volume_mc(zmin_s, zmax_s, M)
    mc_vol <- cv$volume
    mc_se  <- cv$se
    total_volume <- mc_vol
  }

  excess_vacancy <- 1 - exp(-1) - total_volume
  kappa <- excess_vacancy / (1 - exp(-1))
  sd <- sqrt(variance_formula(n, d)) # use variance formula to compute exact variance
  Z <- excess_vacancy * sqrt(n) / sd # standardised statistic
  pval <- 1 - pnorm(Z)
  out <- list(stat = kappa, pval = pval, method = method)
  if (method == 'approx') {
    out$mc_se <- mc_se
  }else {
    out$mc_se <- 0
  }
  return(out)
}
