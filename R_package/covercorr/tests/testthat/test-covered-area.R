test_that("C_covered_area: single rectangle is its area", {
  xmin <- 1; xmax <- 3; ymin <- 2; ymax <- 5
  got <- .Call("C_covered_area", xmin, xmax, ymin, ymax, PACKAGE = "covercorr")
  expect_equal(as.numeric(got), (xmax - xmin) * (ymax - ymin))
})

test_that("C_covered_area: disjoint rectangles sum", {
  xmin <- c(0, 10); xmax <- c(2, 12)
  ymin <- c(0, 10); ymax <- c(3, 12)
  ref <- reference_area_exact(xmin, xmax, ymin, ymax)
  got <- .Call("C_covered_area", xmin, xmax, ymin, ymax, PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

test_that("C_covered_area: overlapping rectangles", {
  xmin <- c(0, 1); xmax <- c(3, 4)
  ymin <- c(0, 1); ymax <- c(3, 4)
  ref <- reference_area_exact(xmin, xmax, ymin, ymax)
  got <- .Call("C_covered_area", xmin, xmax, ymin, ymax, PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

test_that("C_covered_area: nested rectangles", {
  xmin <- c(0, 0.5); xmax <- c(5, 4)
  ymin <- c(0, 0.5); ymax <- c(5, 4)
  ref <- reference_area_exact(xmin, xmax, ymin, ymax)
  got <- .Call("C_covered_area", xmin, xmax, ymin, ymax, PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

test_that("C_covered_area: touching edges (no double-count)", {
  xmin <- c(0, 2); xmax <- c(2, 4)
  ymin <- c(0, 0); ymax <- c(2, 2)
  ref <- reference_area_exact(xmin, xmax, ymin, ymax)
  got <- .Call("C_covered_area", xmin, xmax, ymin, ymax, PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

test_that("C_covered_area: zero-width/height rectangles ignored", {
  xmin <- c(0, 1, 2); xmax <- c(2, 1, 5)   # middle one has width 0
  ymin <- c(0, 1, 2); ymax <- c(2, 4, 2)   # last one has height 0
  ref <- reference_area_exact(xmin, xmax, ymin, ymax)
  got <- .Call("C_covered_area", xmin, xmax, ymin, ymax, PACKAGE = "covercorr")
  expect_equal(as.numeric(got), ref)
})

set.seed(1)
test_that("C_covered_area: randomized rectangles agree with reference", {
  rnd_case <- function(n) {
    a <- runif(n, -1, 2); b <- runif(n, -1, 2)
    c <- runif(n, -1, 2); d <- runif(n, -1, 2)
    xmin <- pmin(a, b); xmax <- pmax(a, b)
    ymin <- pmin(c, d); ymax <- pmax(c, d)
    ref <- reference_area_exact(xmin, xmax, ymin, ymax)
    got <- .Call("C_covered_area", xmin, xmax, ymin, ymax, PACKAGE = "covercorr")
    expect_equal(as.numeric(got), ref, tolerance = 1e-12)
  }
  # A few sizes
  for (n in c(10L, 20L, 50L, 200L)) rnd_case(n)
})
