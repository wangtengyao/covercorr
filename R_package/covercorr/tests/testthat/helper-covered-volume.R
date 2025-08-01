# Exact union volume in d dims by grid partition on unique edges.
# zmin, zmax are n x d numeric matrices.
reference_volume_exact <- function(zmin, zmax) {
  stopifnot(is.matrix(zmin), is.matrix(zmax),
            nrow(zmin) == nrow(zmax), ncol(zmin) == ncol(zmax))
  n <- nrow(zmin); d <- ncol(zmin)
  if (n == 0L) return(0)
  
  # Build unique, sorted edge coordinates for each dimension
  edges <- vector("list", d)
  for (j in seq_len(d)) {
    ej <- sort(unique(c(zmin[, j], zmax[, j])))
    if (length(ej) < 2L) return(0)  # nothing with positive width
    edges[[j]] <- ej
  }
  
  # Precompute per-dim widths (adjacent differences)
  widths <- lapply(edges, function(e) diff(e))
  
  # Generate the index grid for all elementary hypercells
  # For dim j there are (length(edges[[j]]) - 1) cells.
  idx_lists <- lapply(edges, function(e) seq_len(length(e) - 1L))
  # Keep it small; expand.grid is fine for d up to 4 in unit tests
  grid <- do.call(expand.grid, idx_lists)
  # grid is (#cells) x d; each row is (i1, i2, ..., id)
  
  # For each cell, check if covered by any hyperrectangle.
  vol <- 0
  for (r in seq_len(nrow(grid))) {
    # cell bounds per dimension
    i <- as.integer(grid[r, ])
    # lower/upper bounds in each dim
    lows  <- mapply(function(e, k) e[k],       edges, i)
    highs <- mapply(function(e, k) e[k + 1L],  edges, i)
    
    # Positive widths?
    w <- mapply(function(wd, k) wd[k], widths, i)
    if (any(w <= 0)) next
    
    # Cell is covered if any rectangle fully covers the cell in ALL dims:
    # zmin[k, j] <= lows[j] and zmax[k, j] >= highs[j] for all j
    # Vectorize: for each rect k, test all dims
    covered <- FALSE
    for (k in seq_len(n)) {
      if (all(zmin[k, ] <= lows & zmax[k, ] >= highs)) {
        covered <- TRUE
        break
      }
    }
    if (covered) vol <- vol + prod(w)
  }
  vol
}
