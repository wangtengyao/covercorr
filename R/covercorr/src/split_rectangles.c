#include <R.h>
#include <Rinternals.h>
#include <math.h>

/* Enumerate shifts in {-1,0,1}^d; clip into [0,1]^d; keep non-empty intersections */
SEXP C_split_rectangles(SEXP Rzmin, SEXP Rzmax, SEXP Rn, SEXP Rd) {
  int n = INTEGER(Rn)[0], d = INTEGER(Rd)[0];
  const double *zmin = REAL(Rzmin), *zmax = REAL(Rzmax);

  // worst-case number of pieces is 3^d * n
  int maxp = n;
  for (int i=0;i<d;++i) maxp *= 3;

  SEXP RZmin = PROTECT(allocMatrix(REALSXP, maxp, d));
  SEXP RZmax = PROTECT(allocMatrix(REALSXP, maxp, d));
  double *Zmin = REAL(RZmin), *Zmax = REAL(RZmax);
  int count = 0;

  // iterate all shifts encoded in base-3
  int totalShift = 1; for (int i=0;i<d;++i) totalShift *= 3;
  int *shift = (int*) R_alloc(d, sizeof(int));

  for (int s=0; s<totalShift; ++s) {
    int t = s;
    for (int j=0;j<d;++j) { shift[j] = (t % 3) - 1; t /= 3; }

    for (int i=0;i<n;++i) {
      int ok = 1;
      double *tmpmin = (double*) R_alloc(d, sizeof(double));
      double *tmpmax = (double*) R_alloc(d, sizeof(double));
      for (int j=0;j<d;++j) {
        double a = zmin[i + j*n] + shift[j];
        double b = zmax[i + j*n] + shift[j];
        if (a < 0.0) a = 0.0;
        if (b > 1.0) b = 1.0;
        if (a >= b) { ok = 0; break; }
        tmpmin[j] = a; tmpmax[j] = b;
      }
      if (ok) {
        for (int j=0;j<d;++j) {
          Zmin[count + j*maxp] = tmpmin[j];
          Zmax[count + j*maxp] = tmpmax[j];
        }
        ++count;
      }
    }
  }

  // shrink to actual count
  SEXP RZmin2 = PROTECT(allocMatrix(REALSXP, count, d));
  SEXP RZmax2 = PROTECT(allocMatrix(REALSXP, count, d));
  double *Zmin2 = REAL(RZmin2), *Zmax2 = REAL(RZmax2);
  for (int j=0;j<d;++j) {
    for (int i=0;i<count;++i) {
      Zmin2[i + count*j] = Zmin[i + maxp*j];
      Zmax2[i + count*j] = Zmax[i + maxp*j];
    }
  }

  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, RZmin2);
  SET_VECTOR_ELT(out, 1, RZmax2);
  UNPROTECT(5);
  return out;
}
