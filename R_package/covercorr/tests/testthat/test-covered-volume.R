test_that("C_covered_volume: basic 3D single box", {
  d <- 3L
  zmin <- matrix(c(0,0,0), ncol = d)
  zmax <- matrix(c(2,3,5), ncol = d)
  n <- nrow(zmin)
  ref <- reference_volume_exact(zmin, zmax)
  got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

test_that("C_covered_volume: 3D disjoint boxes sum", {
  d <- 3L
  zmin <- rbind(c(0,0,0),  c(3,0,0))
  zmax <- rbind(c(1,2,4),  c(4,1,2))
  n <- nrow(zmin)
  ref <- reference_volume_exact(zmin, zmax)
  got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

test_that("C_covered_volume: 3D overlapping boxes", {
  d <- 3L
  zmin <- rbind(c(0,0,0), c(0.5, 1,  2))
  zmax <- rbind(c(2,3,4), c(1.5, 2.5,3.5))
  n <- nrow(zmin)
  ref <- reference_volume_exact(zmin, zmax)
  got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref, tolerance = 1e-12)
})

test_that("C_covered_volume: 3D touching faces (no double-count)", {
  d <- 3L
  zmin <- rbind(c(0,0,0), c(1,0,0))
  zmax <- rbind(c(1,1,1), c(2,1,1))  # touch at x=1 plane
  n <- nrow(zmin)
  ref <- reference_volume_exact(zmin, zmax)
  got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

test_that("C_covered_volume: 3D degenerate boxes (zero width ignored)", {
  d <- 3L
  zmin <- rbind(c(0,0,0), c(2,2,2), c(5, 1, 1))
  zmax <- rbind(c(2,3,4), c(2,3,4), c(6, 1, 3))  # middle has zero x-width; third has zero y-width segment
  n <- nrow(zmin)
  ref <- reference_volume_exact(zmin, zmax)
  got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

set.seed(42)
test_that("C_covered_volume: randomized small 3D cases agree with reference", {
  d <- 3L
  rnd_case <- function(n) {
    a <- matrix(runif(n * d, -1, 2), ncol = d)
    b <- matrix(runif(n * d, -1, 2), ncol = d)
    zmin <- pmin(a, b)
    zmax <- pmax(a, b)
    ref <- reference_volume_exact(zmin, zmax)
    got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
    expect_equal(as.numeric(got), ref, tolerance = 1e-12)
  }
  for (n in c(10L, 20L, 30L, 50L)) rnd_case(n)
})

# ---------------- 4D tests ----------------

test_that("C_covered_volume: basic 4D single box", {
  d <- 4L
  zmin <- matrix(c(0,0,0,0), ncol = d)
  zmax <- matrix(c(1,2,3,4), ncol = d)
  n <- nrow(zmin)
  ref <- reference_volume_exact(zmin, zmax)
  got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

test_that("C_covered_volume: 4D overlapping boxes", {
  d <- 4L
  zmin <- rbind(c(0,0,0,0), c(0.5, 1,  1,  2))
  zmax <- rbind(c(2,3,4,5), c(1.5, 2.5,3.5,3.2))
  n <- nrow(zmin)
  ref <- reference_volume_exact(zmin, zmax)
  got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref, tolerance = 1e-12)
})

test_that("C_covered_volume: 4D disjoint boxes", {
  d <- 4L
  zmin <- rbind(c(0,0,0,0), c(2,2,2,2))
  zmax <- rbind(c(1,1,1,1), c(3,3,3,3))
  n <- nrow(zmin)
  ref <- reference_volume_exact(zmin, zmax)
  got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

set.seed(7)
test_that("C_covered_volume: randomized small 4D cases agree with reference", {
  d <- 4L
  rnd_case <- function(n) {
    a <- matrix(runif(n * d, -0.5, 1.5), ncol = d)
    b <- matrix(runif(n * d, -0.5, 1.5), ncol = d)
    zmin <- pmin(a, b)
    zmax <- pmax(a, b)
    ref <- reference_volume_exact(zmin, zmax)
    got <- .Call("C_covered_volume", zmin, zmax, as.integer(n), as.integer(d), PACKAGE = "covercorr")
    expect_equal(as.numeric(got), ref, tolerance = 1e-12)
  }
  for (n in c(10L, 20L)) rnd_case(n)  # keep small to control grid size
})
