#' Compute Monge-Kantorovich ranks
#' @param X data points
#' @param U reference points
#' @export
MK_rank <- function(X, U){
  n <- nrow(X)
  if (ncol(X) == 1){
    return(sort(U)[rank(X, tie.method='random')])
  }

  cost <- matrix(0, n, n)
  for (j in seq_len(ncol(X))) {
    cost <- cost + (outer(X[,j], U[,j], "-"))^2
  }
  
  a <- b <- rep(1/n, n)
  ret <- transport::transport(a = a, b = b, costm=cost, p = 2, method='networkflow')
  assignment <- ret$to
  return(U[assignment, ])
}

#' exact formula for variance of test stat
#' @param n sample size
#' @param d dimension
#' @export
variance_formula <- function(n, d){
  V <- ((1 - 2 / n)^n - (1 - 1 / n)^(2 * n)) * n
  multiplier <- 1.0
  inv_factorial <- 1.0
  
  if (n >= 1) {
    for (s in 1:n) {
      multiplier <- multiplier * (1 - (s - 1) / n)
      inv_factorial <- inv_factorial / s
      V <- V + multiplier * inv_factorial * (1 - 2 / n)^(n - s) * (2 / (s + 1))^d
    }
  }
  return(V)
}

#' Split rectangels by wrapping them around edges of [0,1]^d
#' @param zmin n x d matrix of bottomleft coordinates, one row per rectangle
#' @param zmax n x d matrix of topright coordinates, one row per rectangle
#' @return a list of zmin and zmax, describing the bottomleft and topright 
#' coordinates of splitted rectangles
#' @details This is a wrapper of the C_split_rectangles function implemented in C
#' @export
split_rectangles <- function(zmin, zmax){
  splitted <- .Call(C_split_rectangles, as.double(zmin), as.double(zmax),
                    as.integer(nrow(zmin)), as.integer(ncol(zmin)))
  zmin_s <- matrix(splitted[[1L]], ncol = ncol(zmin), byrow = FALSE)
  zmax_s <- matrix(splitted[[2L]], ncol = ncol(zmin), byrow = FALSE)
  return(list(zmin=zmin_s, zmax=zmax_s))
}

#' Total volume of union of rectangles using volume hashing
#' @param zmin n x d matrix of bottomleft coordinates, one row per rectangle
#' @param zmax n x d matrix of topright coordinates, one row per rectangle
#' @returns a numeric value of the volume of the union
#' @details This is a wrapper of the C_covered_volume_partitioned function in C
#' @export
covered_volume_partitioned <- function(zmin, zmax){
  .Call(C_covered_volume_partitioned, as.double(zmin), as.double(zmax),
        as.integer(nrow(zmin)), as.integer(ncol(zmin)))
}

#' Total volume of union of rectangles
#' @param zmin n x d matrix of bottomleft coordinates, one row per rectangle
#' @param zmax n x d matrix of topright coordinates, one row per rectangle
#' @returns a numeric value of the volume of the union
#' @details This is a wrapper of the C_covered_volume_partitioned function in C
#' @export
covered_volume <- function(zmin, zmax){
  .Call(C_covered_volume, as.double(zmin), as.double(zmax),
        as.integer(nrow(zmin)), as.integer(ncol(zmin)))
}

#' Total volume of union of rectangles using Monte Carlo integration
#' @param zmin n x d matrix of bottomleft coordinates, one row per rectangle
#' @param zmax n x d matrix of topright coordinates, one row per rectangle
#' @param M number of Monte Carlo integration sample points
#' @returns a list of the estimated volume of the union and its standard error
#' @details This is a wrapper of the C_covered_volume_mc function in C
#' @export
covered_volume_mc <- function(zmin, zmax, M){
  .Call(C_covered_volume_mc,
               as.double(zmin_s), as.double(zmax_s),
               as.integer(nrow(zmin_s)), as.integer(ncol(zmin_s)),
               M)
}


#' Plot a collection of axis-aligned rectangles in the unit square
#'
#' Draws rectangles specified by their \code{xmin}, \code{xmax}, \code{ymin},
#' and \code{ymax}, optionally adding them to an existing plot. When
#' \code{add = FALSE}, a fresh \eqn{[0,1]\times[0,1]} plot with a grid and
#' equal aspect ratio is created.
#'
#' @param xmin Numeric vector of left x-coordinates.
#' @param xmax Numeric vector of right x-coordinates (same length as \code{xmin}).
#' @param ymin Numeric vector of bottom y-coordinates (same length as \code{xmin}).
#' @param ymax Numeric vector of top y-coordinates (same length as \code{xmin}).
#' @param add Logical; if \code{TRUE}, add to an existing plot. Default \code{FALSE}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of plotting.
#' @export
plot_rectangles <- function(xmin, xmax, ymin, ymax, add=FALSE) {
  
  # Set up plot with equal aspect and a grid
  if (!add) {
    plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "", asp = 1)
    grid()  # add grid lines
  }
  
  # Draw rectangles
  for (i in seq_along(xmin)) {
    rect(xleft = xmin[i], ybottom = ymin[i],
         xright = xmax[i], ytop = ymax[i],
         border = "gray", lwd = 1, col = adjustcolor("gray", alpha.f = 0.2))
  }
  # Draw the [0,1] x [0,1] box
  rect(0, 0, 1, 1, border = "black", lwd = 1, col = NA)
}

#' Coverage-based dependence measure with optional visualisation
#'
#' Constructs small squares around the rank-transformed samples
#' \eqn{(r_x, r_y)} on \eqn{[0,1]^2}, splits them by shifts in
#' \eqn{\{-1,0,1\}^d} and clips to \eqn{[0,1]^d}, and returns the union area of
#' the resulting rectangles. Optionally plots the points and rectangles.
#'
#' @param x Numeric vector or matrix. 
#' @param y Numeric vector or matrix of the same number of rows as \code{x}.
#'   Will be coerced to a matrix and ranked.
#' @param visualise Logical; if \code{TRUE}, scatter the rank points and overlay
#'   rectangles. Default \code{FALSE}.
#'
#' @details Let \eqn{n} be the number of rows after coercion, and
#' \eqn{d = d_x + d_y} the total number of columns across \code{x} and \code{y}.
#' Square half-widths are set to \eqn{\varepsilon = \frac{1}{2} n^{-1/d}} in each
#' dimension. Ranks are scaled to \eqn{(0,1)} using independent sorted uniforms,
#' producing a smoothed rank transform.
#'
#' @return A \code{units} scalar giving the area (coverage) of the union of the
#'   rectangles; use \code{as.numeric()} to obtain a plain numeric.
#'
#' @importFrom stats runif
#'
#' @seealso \code{\link{split_rectangles}}, \code{\link{plot_rectangles}},
#'   \code{\link{union_area}}
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' x <- runif(n)
#' y <- sin(3*x) + runif(n) * 0.01
#' coverage_correlation(x, y, TRUE)
#'
#' @export
coverage_correlation <- function(x, y, visualise=FALSE, 
                                 method=c('auto', 'exact', 'approx'), M=NULL){
  x <- as.matrix(x); y <- as.matrix(y)
  n <- nrow(x); d <- ncol(x) + ncol(y)
  stopifnot(n == nrow(y))
  
  method <- match.arg(method)
  if (method == 'auto') method <- ifelse(d <= 6, 'exact', 'approx')
  
  u <- matrix(runif(n*ncol(x)), n)
  v <- matrix(runif(n*ncol(y)), n)
  x_rank <- MK_rank(x, u)
  y_rank <- MK_rank(y, v)
  eps <- n^(-1/d) / 2
  zmin <- cbind(x_rank - eps, y_rank - eps)
  zmax <- cbind(x_rank + eps, y_rank + eps)
  
  # Wrap around [0,1]^d (split rectangles that cross boundaries); in C
  ret <- split_rectangles(zmin, zmax)
  zmin_s <- ret$zmin
  zmax_s <- ret$zmax

  # --- optional visualisation (first two coord only) ---#
  if (visualise){
    xmin <- zmin_s[, 1]; ymin <- zmin_s[, 2]
    xmax <- zmax_s[ ,1]; ymax <- zmax_s[, 2]
    plot(x_rank, y_rank, pch=20, cex=0.3, asp=1)
    plot_rectangles(xmin, xmax, ymin, ymax, add=TRUE)
  }
  
  # --- Covered volume --- #
  if (method == 'exact'){
    total_volume <- covered_volume_partitioned(zmin_s, zmax_s)
  } else {
    if (is.null(M)) M <- as.integer(ceiling(n^(1.5)))
    ret <- covered_volume_mc(zmin_s, zmax_s, M)
    mc_vol <- ret$volume
    mc_se  <- ret$se
    total_volume <- mc_vol
    attr(total_volume, "mc_se") <- mc_se
  }
  
  kappa <- 1 - (1 - 1/n)^n - total_volume
  sd <- sqrt(variance_formula(n, d)) # use variance formula to compute exact variance
  Z <- kappa * sqrt(n) / sd # standardised statistic
  pval <- 1 - pnorm(Z)
  out <- list(stat = kappa, pval = pval, method = method)
  if (!is.null(attr(total_volume, "mc_se"))) {
    out$mc_se <- attr(total_volume, "mc_se")
  }
  return(out)
}

