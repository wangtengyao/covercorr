#include "covered_area.h"
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>

/* ----- Segment tree over elementary y-segments ----- */
typedef struct {
  int m;            /* number of leaves as power-of-two */
  int *covered;     /* length 2m */
  double *weights;  /* length 2m */
  double *scores;   /* length 2m */
  int *start;       /* length 2m */
  int *end;         /* length 2m */
} segtree_t;

typedef struct { double val; int id; } xv_t;

int cmp_xv(const void *a, const void *b){
  const xv_t *pa=(const xv_t*)a,*pb=(const xv_t*)b;
  if (pa->val < pb->val) return -1;
  if (pa->val > pb->val) return 1;
  return (pa->id < pb->id) ? -1 : 1;
}

static int next_pow2(int n) {
  int p = 1; while (p < n) p <<= 1; return p;
}

static segtree_t* st_build(const double *w, int n) {
  segtree_t *st = (segtree_t*) R_alloc(1, sizeof(segtree_t));
  int m = next_pow2(n);
  st->m = m;
  st->covered = (int*) R_alloc(2*m, sizeof(int));
  st->weights = (double*) R_alloc(2*m, sizeof(double));
  st->scores  = (double*) R_alloc(2*m, sizeof(double));
  st->start   = (int*) R_alloc(2*m, sizeof(int));
  st->end     = (int*) R_alloc(2*m, sizeof(int));
  memset(st->covered, 0, 2*m*sizeof(int));
  memset(st->scores,  0, 2*m*sizeof(double));
  for (int node = 2*m-1; node >= 1; --node) {
    if (node >= m) { 
      int idx = node - m;
      st->weights[node] = (idx < n) ? w[idx] : 0.0;
      st->start[node] = idx;
      st->end[node]   = idx + 1;
    } else {
      st->weights[node] = st->weights[2*node] + st->weights[2*node+1];
      st->start[node]   = st->start[2*node];
      st->end[node]     = st->end[2*node+1];
    }
  }
  return st;
}

static void st_update(segtree_t *st, int node, int l, int r, int mult) {
  // update covered count 
  if (st->start[node] >= r || st->end[node] <= l) return;
  if (st->start[node] >= l && st->end[node] <= r) {
    st->covered[node] += mult;
  } else {
    st_update(st, 2*node,   l, r, mult);
    st_update(st, 2*node+1, l, r, mult);
  }
  // update score 
  if (st->covered[node] != 0) {
    st->scores[node] = st->weights[node];
  } else if (node >= st->m) {
    st->scores[node] = 0.0;
  } else {
    st->scores[node] = st->scores[2*node] + st->scores[2*node+1];
  }
}

double covered_area_native(const double *xmin, const double *xmax,
                           const double *ymin, const double *ymax, int n)
{
  int m2 = 2 * n;

  double *x      = (double*) malloc(m2 * sizeof(double));
  int    *etype  = (int*) malloc(m2 * sizeof(int));
  int    *idx    = (int*) malloc(m2 * sizeof(int));
  xv_t   *xv     = (xv_t*)   malloc(m2 * sizeof(xv_t));
  double *x_sorted = (double*) malloc(m2 * sizeof(double));
  int    *order  = (int*) malloc(m2 * sizeof(int));
  double *ycomb  = (double*) malloc(m2 * sizeof(double));
  xv_t   *yv     = (xv_t*) malloc(m2 * sizeof(xv_t));
  int    *ranky  = (int*) malloc(m2 * sizeof(int));
  double *ysorted = (double*) malloc(m2 * sizeof(double));
  double *elelen = (double*) malloc(m2 * sizeof(double));

  if (!x || !etype || !idx || !xv || !x_sorted || !order || !ycomb ||
      !yv || !ranky || !ysorted || !elelen)
  {
    error("covered_area_native: memory allocation failed");
  }

  for (int i = 0; i < n; ++i) {
    x[i] = xmin[i];
    etype[i] = +1;
    idx[i] = i;
  }
  for (int i = 0; i < n; ++i) {
    x[i + n] = xmax[i];
    etype[i + n] = -1;
    idx[i + n] = i;
  }

  for (int i = 0; i < m2; ++i) {
    xv[i].val = x[i];
    xv[i].id = i;
  }
  qsort(xv, m2, sizeof(xv_t), cmp_xv);
  for (int i = 0; i < m2; ++i) {
    x_sorted[i] = xv[i].val;
    order[i] = xv[i].id;
  }

  for (int i = 0; i < n; ++i) {
    ycomb[i] = ymin[i];
    ycomb[i + n] = ymax[i];
  }

  for (int i = 0; i < m2; ++i) {
    yv[i].val = ycomb[i];
    yv[i].id = i;
  }
  qsort(yv, m2, sizeof(xv_t), cmp_xv);
  for (int r = 0; r < m2; ++r) {
    ranky[yv[r].id] = r;
  }

  int *rank_ymin = ranky;       // First n
  int *rank_ymax = ranky + n;   // Last n (see below note)

  for (int i = 0; i < m2; ++i) ysorted[i] = yv[i].val;
  for (int i = 0; i < m2 - 1; ++i) elelen[i] = ysorted[i + 1] - ysorted[i];
  elelen[m2 - 1] = 0.0;

  segtree_t *st = st_build(elelen, m2);  // assume it allocates internally and will be freed elsewhere

  double area = 0.0;
  for (int i = 0; i < m2 - 1; ++i) {
    int j = order[i];
    int id = idx[j];
    int mult = etype[j];

    int l = (j < n) ? rank_ymin[id] : rank_ymin[id];  // ymin index always in first n
    int r = (j < n) ? rank_ymax[id] : rank_ymax[id];

    st_update(st, 1, l, r, mult);
    double width = x_sorted[i + 1] - x_sorted[i];
    area += st->scores[1] * width;
  }

  /* Free temp memory (except for st if st_build handles it) */
  free(x);
  free(etype);
  free(idx);
  free(xv);
  free(x_sorted);
  free(order);
  free(ycomb);
  free(yv);
  free(ranky);
  free(ysorted);
  free(elelen);

  return area;
}


/* wrapper to call from R for tests */
SEXP C_covered_area(SEXP Rxmin, SEXP Rxmax, SEXP Rymin, SEXP Rymax) {
  int n = LENGTH(Rxmin);
  double area = covered_area_native(REAL(Rxmin), REAL(Rxmax),
                                    REAL(Rymin), REAL(Rymax), n);
  SEXP ans = PROTECT(allocVector(REALSXP, 1));
  REAL(ans)[0] = area;
  UNPROTECT(1);
  return ans;
}
