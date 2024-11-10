#include <R.h>
#include <Rinternals.h>

SEXP xor_merge(SEXP vec1, SEXP vec2) {
  // Lengths and data of inputs
  int n1 = LENGTH(vec1);
  int n2 = LENGTH(vec2);
  int *v1 = INTEGER(vec1);
  int *v2 = INTEGER(vec2);

  // Workspace
  SEXP work = PROTECT(allocVector(INTSXP, n1 + n2));
  int *wk = INTEGER(work);

  int i = 0, j = 0, k = 0;

  // Merge v1 and v2 dropping common elements
  while(i < n1 && j < n2) {
    if (v1[i] < v2[j]) {
      wk[k++] = v1[i++];
    } else if(v1[i] > v2[j]) {
      wk[k++] = v2[j++];
    } else {
      // Skip both if they are equal
      i++;
      j++;
    }
  }

  // Append remaining elements
  while(i < n1) wk[k++] = v1[i++];
  while(j < n2) wk[k++] = v2[j++];

  // Shrink to actual size
  SEXP result = PROTECT(lengthgets(work, k));
  UNPROTECT(2);
  return result;
}


static const R_CallMethodDef CallEntries[] = {
  {"xor_merge", (DL_FUNC) &xor_merge, 2},
  {NULL, NULL, 0}
};

void R_init_sparseZ2(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);  // Disable dynamic lookup for security
}
