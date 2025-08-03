#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "covered_area.h"
#include <math.h>

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
  evt_t *ev = (evt_t*) malloc((size_t)m2 * sizeof(evt_t));
  if (!ev) error("covered_volume_rec: failed to allocate events buffer");

  /* Build events: zmin[,0] are openings; zmax[,0] are closings */
  const double *x0min = zmin + 0 * n;
  const double *x0max = zmax + 0 * n;
  for (int i = 0; i < n; ++i) {
    ev[i]     = (evt_t){ x0min[i], i, 1 };
    ev[i + n] = (evt_t){ x0max[i], i, 0 };
  }
  qsort(ev, (size_t)m2, sizeof(evt_t), cmp_evt);

  /* Active set flags for current slice */
  int *active = (int*) malloc((size_t)n * sizeof(int));
  if (!active) {
      free(ev);
      error("covered_volume_rec: failed to allocate active flags");
    }
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
	size_t slice_len = (size_t)na * (size_t)(d - 1);
    double *smin = (double*) malloc(slice_len * sizeof(double));
	if (!smin) {
	  free(active);
	  free(ev);
	  error("covered_volume_rec: failed to allocate smin slice");
	}
    double *smax = (double*) malloc(slice_len * sizeof(double));
    if (!smax) {
      free(smin);
      free(active);
      free(ev);
      error("covered_volume_rec: failed to allocate smax slice");
    }

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
	
    /* free per-slice buffers before moving to next interval */
    free(smin);
    free(smax);
  }

  /* free per-call buffers */
  free(active);
  free(ev);
  
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

// === Partitioned version === //

/* Linked list node for bucketing */
typedef struct node_t {
  int idx;
  struct node_t *next;
} node_t;

static void powb_fill(int *powb, int b, int d){
  powb[0] = 1;
  for (int j=1; j<d; ++j) powb[j] = powb[j-1] * b;
}

static int cell_floor(double x, double h, int b){
  int v = (int) floor(x / h);
  if (v < 0) v = 0;
  if (v >= b) v = b-1;
  return v;
}

static int linearize(const int *cur, const int *powb, int d){
  int id = 0;
  for (int j=0; j<d; ++j) id += cur[j] * powb[j];
  return id;
}

static void decode_id(int id, int b, int d, int *idx){
  for (int j=0; j<d; ++j){
    idx[j] = id % b;
    id /= b;
  }
}

/* Main entry point: .Call interface */
SEXP C_covered_volume_partitioned(SEXP Rzmin, SEXP Rzmax, SEXP Rn, SEXP Rd)
{
  const int n = INTEGER(Rn)[0];
  const int d = INTEGER(Rd)[0];
  const double *zmin = REAL(Rzmin);
  const double *zmax = REAL(Rzmax);

  int b = (int)floor(pow((double)n, 1.0 / (double)d));

  const double h = 1.0 / (double)b;
  const int B = (int)pow(b, d);  /* total number of blocks */

  node_t **heads = (node_t**) R_alloc(B, sizeof(node_t*));
  for (int i=0; i<B; ++i) heads[i] = NULL;

  int *powb = (int*) R_alloc(d, sizeof(int));
  powb_fill(powb, b, d);

  int *lo_idx = (int*) R_alloc(d, sizeof(int));
  int *hi_idx = (int*) R_alloc(d, sizeof(int));
  int *cur    = (int*) R_alloc(d, sizeof(int));

  /* Insert each rectangle into relevant buckets */
  for (int i=0; i<n; ++i){
    for (int j=0; j<d; ++j){
      double L = zmin[i + j*n];
      double H = zmax[i + j*n];
      int a_lo = cell_floor(L, h, b);
      double tiny = 1e-15;
      int a_hi = cell_floor(H - tiny, h, b);
      if (a_hi < a_lo) a_hi = a_lo;
      lo_idx[j] = a_lo;
      hi_idx[j] = a_hi;
      cur[j] = a_lo;
    }
    while (1){
      int id = linearize(cur, powb, d);
      node_t *nd = (node_t*) R_alloc(1, sizeof(node_t));
      nd->idx = i;
      nd->next = heads[id];
      heads[id] = nd;
      int pos = 0;
      while (pos < d){
        if (cur[pos] < hi_idx[pos]){ cur[pos]++; break; }
        cur[pos] = lo_idx[pos];
        pos++;
      }
      if (pos == d) break;
    }
  }

  double total = 0.0;
  int *idx_vec = (int*) R_alloc(d, sizeof(int));

  for (int id=0; id<B; ++id){
    node_t *p = heads[id];
    if (!p) continue;

    decode_id(id, b, d, idx_vec);
    double *blo = (double*) R_alloc(d, sizeof(double));
    double *bhi = (double*) R_alloc(d, sizeof(double));
    for (int j=0; j<d; ++j){
      blo[j] = (double) idx_vec[j] * h;
      bhi[j] = blo[j] + h;
      if (bhi[j] > 1.0) bhi[j] = 1.0;
    }

    int m_cap = 0;
    for (node_t *q=p; q; q=q->next) ++m_cap;
    if (m_cap == 0) continue;

    double *tmin = (double*) R_alloc(m_cap * d, sizeof(double));
    double *tmax = (double*) R_alloc(m_cap * d, sizeof(double));

    int m = 0;
    for (node_t *q=p; q; q=q->next){
      int i = q->idx;
      int ok = 1;
      for (int j=0; j<d; ++j){
        double L = zmin[i + j*n];
        double H = zmax[i + j*n];
        double Lc = (L > blo[j]) ? L : blo[j];
        double Hc = (H < bhi[j]) ? H : bhi[j];
        if (Hc <= Lc){ ok = 0; break; }
        tmin[m + j*m_cap] = Lc;
        tmax[m + j*m_cap] = Hc;
      }
      if (ok) ++m;
    }
    if (m == 0) continue;

    double *umin = (double*) R_alloc(m * d, sizeof(double));
    double *umax = (double*) R_alloc(m * d, sizeof(double));
    for (int j=0; j<d; ++j){
      for (int r=0; r<m; ++r){
        umin[r + j*m] = tmin[r + j*m_cap];
        umax[r + j*m] = tmax[r + j*m_cap];
      }
    }

    double vol = covered_volume_rec(umin, umax, m, d);
    total += vol;
  }

  SEXP out = PROTECT(allocVector(REALSXP, 1));
  REAL(out)[0] = total;
  UNPROTECT(1);
  return out;
}