data_generation <- function(n, type, noise=1, dx=1){
  eps_x <- matrix(rnorm(n*dx), n, dx)
  eps_y <- matrix(rnorm(n*dx), n, dx)
  if (type == 'sinusoidal'){
    x <- matrix(runif(n*dx, -1, 1), n, dx)
    y <- cos(8 * pi * x) + eps_y * noise
  } else if (type == 'zigzag'){
    x <- matrix(runif(n*dx, -1, 1), n, dx)
    y <- abs(x - 0.5 * sign(x)) + eps_y * noise
  } else if (type == 'circle'){
    u <- matrix(runif(n*dx, 0, 2 * pi), n, dx)
    x <- cos(u) + 0.5 * eps_x * noise
    y <- sin(u) + 0.5 * eps_y * noise
  } else if (type == 'spiral'){
    u <- matrix(runif(n*dx, 0, 1), n, dx)
    x <- u * sin(10 * pi * u) + eps_x * noise * 0.15
    y <- u * cos(10 * pi * u) + eps_y * noise * 0.15
  } else if (type == 'Lissajous'){
    u <- matrix(runif(n*dx, 0, 20), n, dx)
    x = sin(3 * u + pi / 2) + eps_x * noise * 0.1
    y = sin(4 * u) + eps_y * noise * 0.1
  } else if (type == 'local'){
    x <- matrix(rnorm(n*dx), n, dx) * 0.25
    y <- matrix(rnorm(n*dx), n, dx) * 0.25
    ind <- x > 0 & y > 0
    y[ind] <- x[ind]
    x <- x + eps_x * noise * 0.2
    y <- y + eps_y * noise * 0.2
  } else if (type == 'indep'){
    x <- rnorm(n)
    y <- rnorm(n)
  }
  return(list(x=x, y=y))
}

