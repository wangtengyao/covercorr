#' @useDynLib covercorr, .registration = TRUE
NULL

#' @importFrom grDevices adjustcolor
#' @importFrom graphics grid rect
#' @importFrom stats complete.cases pnorm

# Suppress NOTES for data.frame columns
utils::globalVariables(c("zmin_s", "zmax_s"))


#' Monge–Kantorovich ranks (uniform OT via squared distances)
#'
#' Computes the optimal matching that maps each observation in \code{X} to a
#' reference point in \code{U} using uniform weights and squared Euclidean cost.
#' Internally uses \code{transport::transport(method = "networkflow", p = 2)}.
#' In 1D, this reduces to a rank-based matching
#' \code{sort(U)[rank(X, ties.method = "random")]}.
#'
#' @param X Numeric vector of length \eqn{n}, or numeric matrix with \eqn{n}
#'   rows and \eqn{d} columns. If not a matrix, it is coerced with
#'   \code{as.matrix()}.
#' @param U Numeric vector of length \eqn{n}, or numeric matrix with \eqn{n}
#'   rows and \eqn{d} columns. If not a matrix, it is coerced with
#'   \code{as.matrix()}. Must have the same number of rows as \code{X}.
#'
#' @return If \code{ncol(X) == 1}, a numeric vector of length \eqn{n}
#'   containing the entries of \code{U} reordered to match the ranks of
#'   \code{X}. Otherwise, a numeric \eqn{n \times d} matrix whose \eqn{i}-th row
#'   is the matched row of \code{U} corresponding to the \eqn{i}-th row of
#'   \code{X}.
#'
#' @details
#' - Rows must match: \code{nrow(X) == nrow(U)} (otherwise an error is thrown).
#' - Columns must match: \code{ncol(X) == ncol(U)} (otherwise an error is thrown).
#' - Weights are uniform (\eqn{1/n}) and the cost matrix is the sum of squared
#'   coordinate differences across columns.
#' - In 1D, ties in \code{X} are broken at random via
#'   \code{ties.method = "random"}; use \code{set.seed()} for reproducibility.
#'
#' @section Dependencies:
#' Requires the \pkg{transport} package.
#'
#' @examples
#' # 1D example (set seed for reproducible tie-breaking)
#' set.seed(1)
#' x <- rnorm(10)
#' u <- seq(0, 1, length.out = 10)
#' MK_rank(x, u)
#'
#' # 2D example
#' set.seed(42)
#' X <- matrix(rnorm(200), ncol = 2)   # 100 x 2
#' U <- matrix(runif(200),  ncol = 2)  # 100 x 2
#' R <- MK_rank(X, U)
#' dim(R)  # 100 2
#'
#' @importFrom transport transport
#' @export
MK_rank <- function(X, U){
  # if inputs are not in proper matrix format change if possible
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(U)) {
    U = as.matrix(U)
  }

  if (nrow(X) != nrow(U)) stop("Number of rows of U and X should be equal.")
  if (ncol(X) != ncol(U)) stop("Number of columns of U and X should be equal.")

  n <- nrow(X)
  if (ncol(X) == 1){
    return(sort(U)[rank(X, ties.method='random')])
  }

  cost <- matrix(0, n, n)
  for (j in seq_len(ncol(X))) {
    cost <- cost + (outer(X[,j], U[,j], "-"))^2
  }

  a <- b <- rep(1/n, n)
  ret <- transport::transport(a = a, b = b, costm=cost, p = 2, method = 'networkflow')
  assignment <- ret$to
  return(U[assignment, ])
}

#' Variance of the the excess vacancy
#'
#' Exact formula for \eqn{n} times the variance of the excess vacancy.
#' For independent \eqn{X} and \eqn{Y}, the variance of the coverage correlation
#' coefficient is obtained by dividing the returned value by \eqn{n(1 - e^{-1})^2}.
#' check \href{https://arxiv.org/abs/2508.06402}{the arXiv preprint} for more details
#' @param n sample size
#' @param d dimension \eqn{(X, Y)}
#' @return variance formula in paper
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

#' Split rectangles by wrapping them around edges of \eqn{[0,1]^d}
#' @param zmin n x d matrix of bottom-left coordinates, one row per rectangle
#' @param zmax n x d matrix of top-right coordinates, one row per rectangle
#' @return a list of zmin and zmax, describing the bottom-left and top-right
#' coordinates of splitted rectangles
#' @details This is a wrapper of the C_split_rectangles function implemented in C
#' @export
split_rectangles <- function(zmin, zmax){
  splitted <- .Call("C_split_rectangles", as.double(zmin), as.double(zmax),
                    as.integer(nrow(zmin)), as.integer(ncol(zmin)))
  zmin_s <- matrix(splitted[[1L]], ncol = ncol(zmin), byrow = FALSE)
  zmax_s <- matrix(splitted[[2L]], ncol = ncol(zmin), byrow = FALSE)
  return(list(zmin = zmin_s, zmax = zmax_s))
}

#' Total volume of union of rectangles using volume hashing
#' @param zmin n x d matrix of bottomleft coordinates, one row per rectangle
#' @param zmax n x d matrix of topright coordinates, one row per rectangle
#' @returns a numeric value of the volume of the union
#' @details This is a wrapper of the C_covered_volume_partitioned function in C
#' @export
covered_volume_partitioned <- function(zmin, zmax){
  .Call("C_covered_volume_partitioned", as.double(zmin), as.double(zmax),
        as.integer(nrow(zmin)), as.integer(ncol(zmin)))
}

#' Total volume of union of rectangles
#' @param zmin n x d matrix of bottomleft coordinates, one row per rectangle
#' @param zmax n x d matrix of topright coordinates, one row per rectangle
#' @returns a numeric value of the volume of the union
#' @details This is a wrapper of the C_covered_volume_partitioned function in C
#' @export
covered_volume <- function(zmin, zmax){
  .Call("C_covered_volume", as.double(zmin), as.double(zmax),
        as.integer(nrow(zmin)), as.integer(ncol(zmin)))
}

#' Total volume of union of rectangles using Monte Carlo integration
#' @param zmin_s n x d matrix of bottomleft coordinates, one row per rectangle
#' @param zmax_s n x d matrix of topright coordinates, one row per rectangle
#' @param M number of Monte Carlo integration sample points
#' @returns a list of the estimated volume of the union and its standard error
#' @details This is a wrapper of the C_covered_volume_mc function in C
#' @export
covered_volume_mc <- function(zmin_s, zmax_s, M){
  .Call("C_covered_volume_mc",
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
#' @return Invisibly returns \code{NULL}. Use this function for its plotting output, not for a returned value.
#' @export
plot_rectangles <- function(xmin, xmax, ymin, ymax, add = FALSE) {

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
