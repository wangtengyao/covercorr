#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* Monte Carlo approximation of covered volume in [0,1]^d.
   Inputs:
     - Rzmin, Rzmax: numeric matrices (n x d), column-major
     - Rn, Rd, RM: integers (n, d, M) — number of rectangles, dimension, and MC sample size
   Returns:
     list(volume = p_hat, se = sqrt(p_hat*(1-p_hat)/M))
*/
SEXP C_covered_volume_mc(SEXP Rzmin, SEXP Rzmax, SEXP Rn, SEXP Rd, SEXP RM)
{
  const int n = INTEGER(Rn)[0];
  const int d = INTEGER(Rd)[0];
  const int M = INTEGER(RM)[0];

  const double *zmin = REAL(Rzmin);
  const double *zmax = REAL(Rzmax);

  double *u = (double*) R_alloc(d, sizeof(double));  // temporary buffer

  long hits = 0;

  GetRNGstate();
  for (int m = 0; m < M; ++m) {
    for (int j = 0; j < d; ++j)
      u[j] = unif_rand();

    int covered = 0;
    for (int i = 0; i < n; ++i) {
      int inside = 1;
      for (int j = 0; j < d; ++j) {
        const double a = zmin[i + n*j];
        const double b = zmax[i + n*j];
        if (u[j] < a || u[j] > b) {
          inside = 0;
          break;
        }
      }
      if (inside) {
        covered = 1;
        break;
      }
    }
    if (covered) ++hits;
  }
  PutRNGstate();

  const double p_hat = (double)hits / (double)M;
  const double se = sqrt(p_hat * (1.0 - p_hat) / (double)M);

  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, ScalarReal(p_hat));
  SET_VECTOR_ELT(out, 1, ScalarReal(se));
  SEXP nm = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(nm, 0, mkChar("volume"));
  SET_STRING_ELT(nm, 1, mkChar("se"));
  setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}
