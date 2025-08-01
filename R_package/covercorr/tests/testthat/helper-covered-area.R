# tests/testthat/helper-covered-area.R
reference_area_exact <- function(xmin, xmax, ymin, ymax) {
  stopifnot(length(xmin) == length(xmax),
            length(ymin) == length(ymax),
            length(xmin) == length(ymin))
  n <- length(xmin)
  if (n == 0L) return(0)
  
  xs <- sort(unique(c(xmin, xmax)))
  ys <- sort(unique(c(ymin, ymax)))
  if (length(xs) < 2L || length(ys) < 2L) return(0)
  
  area <- 0
  for (ix in seq_len(length(xs) - 1L)) {
    x0 <- xs[ix]; x1 <- xs[ix + 1L]
    wx <- x1 - x0
    if (wx <= 0) next
    # Which rectangles fully cover this vertical strip?
    covered_x <- (xmin <= x0) & (xmax >= x1)
    
    if (!any(covered_x)) next
    ymn <- ymin[covered_x]; ymx <- ymax[covered_x]
    # 1D union length along y inside this strip, exact on the y-edge grid
    for (iy in seq_len(length(ys) - 1L)) {
      y0 <- ys[iy]; y1 <- ys[iy + 1L]
      wy <- y1 - y0
      if (wy <= 0) next
      # Cell [x0,x1] x [y0,y1] is covered iff any rect fully covers it
      if (any(ymn <= y0 & ymx >= y1)) {
        area <- area + wx * wy
      }
    }
  }
  area
}
