#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "covered_area.h"

/* ---------------------------
   Internal types / comparators
   --------------------------- */

typedef struct {
  double val;  /* coordinate value on the sweep axis */
  int    id;   /* rectangle index 0..n-1 */
  int    open; /* 1 = entering (zmin), 0 = leaving (zmax) */
} evt_t;

/* Sort by coordinate; for ties, process openings before closings. */
static int cmp_evt(const void *a, const void *b) {
  const evt_t *pa = (const evt_t*)a, *pb = (const evt_t*)b;
  if (pa->val < pb->val) return -1;
  if (pa->val > pb->val) return 1;
  /* same coordinate: open before close */
  if (pa->open > pb->open) return -1;
  if (pa->open < pb->open) return 1;
  /* tie-breaker on id just for determinism */
  return (pa->id < pb->id) ? -1 : (pa->id > pb->id);
}

/* ------------------------------------
   Recursive union volume (no SEXP here)
   zmin, zmax: column-major n x d arrays
   ------------------------------------ */
static double covered_volume_rec(const double *zmin,
                                 const double *zmax,
                                 int n, int d)
{
  if (n <= 0) return 0.0;

  if (d == 2) {
    const double *xmin = zmin + 0 * n;
    const double *xmax = zmax + 0 * n;
    const double *ymin = zmin + 1 * n;
    const double *ymax = zmax + 1 * n;
    return covered_area_native(xmin, xmax, ymin, ymax, n);
  }

  /* ---------- d >= 3: sweep on the first coord ---------- */

  const int m2 = 2 * n;
  evt_t *ev = (evt_t*) R_alloc(m2, sizeof(evt_t));

  /* Build events: zmin[,0] are openings; zmax[,0] are closings */
  const double *x0min = zmin + 0 * n;
  const double *x0max = zmax + 0 * n;
  for (int i = 0; i < n; ++i) {
    ev[i]     = (evt_t){ x0min[i], i, 1 };
    ev[i + n] = (evt_t){ x0max[i], i, 0 };
  }
  qsort(ev, m2, sizeof(evt_t), cmp_evt);

  /* Active set flags for current slice */
  int *active = (int*) R_alloc(n, sizeof(int));
  for (int i = 0; i < n; ++i) active[i] = 0;

  double volume = 0.0;

  for (int k = 0; k < m2 - 1; ++k) {
    /* Apply event k (toggle active set) */
    active[ev[k].id] = ev[k].open ? 1 : 0;

    /* Compute width to next event */
    double width = ev[k + 1].val - ev[k].val;
    if (width <= 0.0) continue;

    /* Count active rectangles */
    int na = 0;
    for (int i = 0; i < n; ++i) if (active[i]) ++na;
    if (na == 0) continue;

    /* Build slice arrays smin, smax in (d-1) dims, column-major (na x (d-1)) */
    double *smin = (double*) R_alloc(na * (d - 1), sizeof(double));
    double *smax = (double*) R_alloc(na * (d - 1), sizeof(double));

    int pos = 0;
    /* For each original rectangle i that is active, copy its coords j=1..d-1 */
    for (int i = 0; i < n; ++i) if (active[i]) {
      for (int j = 1; j < d; ++j) {
        /* destination index: (pos) + (j-1)*na */
        smin[pos + (j - 1) * na] = zmin[i + j * n];
        smax[pos + (j - 1) * na] = zmax[i + j * n];
      }
      ++pos;
    }

    /* Recurse on (d-1)-dimensional slice */
    double area = covered_volume_rec(smin, smax, na, d - 1);
    volume += area * width;
  }

  return volume;
}

/* ------------------------------
   R-visible .Call() entry point
   ------------------------------ */
SEXP C_covered_volume(SEXP Rzmin, SEXP Rzmax, SEXP Rn, SEXP Rd)
{
  const int n = INTEGER(Rn)[0];
  const int d = INTEGER(Rd)[0];

  /* Expect zmin, zmax as numeric matrices with column-major layout n x d */
  const double *zmin = REAL(Rzmin);
  const double *zmax = REAL(Rzmax);

  double vol = covered_volume_rec(zmin, zmax, n, d);

  SEXP Rout = PROTECT(allocVector(REALSXP, 1));
  REAL(Rout)[0] = vol;
  UNPROTECT(1);
  return Rout;
}
