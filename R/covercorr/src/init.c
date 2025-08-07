#include <R.h>
#include <Rinternals.h>

SEXP C_covered_area(SEXP xmin, SEXP xmax, SEXP ymin, SEXP ymax);
SEXP C_covered_volume(SEXP zmin, SEXP zmax, SEXP Rn, SEXP Rd);
SEXP C_split_rectangles(SEXP zmin, SEXP zmax, SEXP Rn, SEXP Rd);
SEXP C_covered_volume_mc(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_covered_volume_partitioned(SEXP zmin, SEXP zmax, SEXP Rn, SEXP Rd);

static const R_CallMethodDef CallEntries[] = {
  {"C_covered_area",      (DL_FUNC) &C_covered_area,     4},
  {"C_covered_volume",    (DL_FUNC) &C_covered_volume,   4},
  {"C_split_rectangles",  (DL_FUNC) &C_split_rectangles, 4},
  {"C_covered_volume_mc", (DL_FUNC) &C_covered_volume_mc, 5},
  {"C_covered_volume_partitioned",(DL_FUNC) &C_covered_volume_partitioned,4},
  {NULL, NULL, 0}
};

void R_init_covercorr(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
