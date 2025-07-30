library(sf)

#' Split axis-aligned rectangles mod [0,1]^d
#' Given lower and upper corners of axis-aligned rectangles in R^d, this 
#' function generates the image of these rectangles modulo [0,1]^d (i.e. 
#' viewing [0,1]^d with periodic boundaries)
#'
#' @param zmin A numeric matrix of size \eqn{n \times d} containing the lower
#'   bounds (one rectangle per row) before shifting.
#' @param zmax A numeric matrix of size \eqn{n \times d} containing the upper
#'   bounds (one rectangle per row) before shifting. Must satisfy
#'   \code{zmin < zmax} elementwise.
#'
#' @return A list with two numeric matrices:
#' \describe{
#'   \item{\code{zmin_splitted}}{Matrix of lower bounds of the non-empty
#'   intersections with \eqn{[0,1]^d}.}
#'   \item{\code{zmax_splitted}}{Matrix of upper bounds of the non-empty
#'   intersections with \eqn{[0,1]^d}.}
#' }
#'
#' @export
split_rectangles <- function(zmin, zmax) {
  d <- ncol(zmin)
  
  # all 3^d shifts in {-1, 0, 1}^d
  shifts_df <- expand.grid(rep(list(c(-1, 0, 1)), d))
  
  zmin_splitted <- zmax_splitted <- matrix(nrow=0, ncol=d)
  
  for (i in seq_len(nrow(shifts_df))) {
    shift <- as.numeric(shifts_df[i, ])
    
    # shift rectangles and intersect with [0,1]^d
    zmin_shifted <- pmax(sweep(zmin, 2, shift, `+`), 0)
    zmax_shifted <- pmin(sweep(zmax, 2, shift, `+`), 1)
    
    # keep rows where intersection is non-empty (strictly positive side lengths)
    idx <- apply(zmin_shifted < zmax_shifted, 1, all)
    zmin_splitted <- rbind(zmin_splitted, zmin_shifted[idx, , drop=FALSE])
    zmax_splitted <- rbind(zmax_splitted, zmax_shifted[idx, , drop=FALSE])
  }
  
  list(zmin_splitted = zmin_splitted, zmax_splitted = zmax_splitted)
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

#' Area of the union of axis-aligned rectangles via \pkg{sf}
#'
#' Builds simple polygons from the supplied rectangles, unions them, and returns
#' their total area.
#'
#' @param xmin Numeric vector of left x-coordinates.
#' @param xmax Numeric vector of right x-coordinates (same length as \code{xmin}).
#' @param ymin Numeric vector of bottom y-coordinates (same length as \code{xmin}).
#' @param ymax Numeric vector of top y-coordinates (same length as \code{xmin}).
#'
#' @return Area of the union of the rectangles. 
#'
#' @importFrom sf st_polygon st_area st_union st_sfc
#' @export
union_area <- function(xmin, xmax, ymin, ymax){
  polys <- lapply(seq_along(xmin), function(i) {
    m <- matrix(
      c(xmin[i], ymin[i],
        xmax[i], ymin[i],
        xmax[i], ymax[i],
        xmin[i], ymax[i],
        xmin[i], ymin[i]),
      ncol = 2, byrow = TRUE
    )
    st_polygon(list(m))
  })
  
  st_area(st_union(st_sfc(polys)))
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
coverage_correlation <- function(x, y, visualise=FALSE){
  x <- as.matrix(x); y <- as.matrix(y)
  n <- nrow(x); d <- ncol(x) + ncol(y)

  u <- sort(runif(n))
  v <- sort(runif(n))
  rx <- u[rank(x)]
  ry <- v[rank(y)]
  
  eps_x <- eps_y <- n^(-1/d) / 2
  zmin <- cbind(rx - eps_x, ry - eps_y)
  zmax <- cbind(rx + eps_x, ry + eps_y)
  
  tmp <- split_rectangles(zmin, zmax)
  
  xmin <- tmp$zmin_splitted[, 1]
  ymin <- tmp$zmin_splitted[, 2]
  xmax <- tmp$zmax_splitted[ ,1]
  ymax <- tmp$zmax_splitted[, 2]
  
  if (visualise){
    plot(rx, ry, pch=20, cex=0.3, asp=1)
    plot_rectangles(xmin, xmax, ymin, ymax, add=TRUE)
  }
  
  union_area(xmin, xmax, ymin, ymax)
}
