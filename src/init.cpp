// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rversion.h>

static const R_CallMethodDef CallEntries[] = {
  {NULL, NULL, 0}
};

void R_init_eggCounts(DllInfo* info) {
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
