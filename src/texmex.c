#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* accurately compute (exp(x) - 1) / x */

static inline
double true_dexprl(const double x) {
  if (R_FINITE(x)) {
    /* actual value */
    if (x != 0) {
      return expm1(x) / x;
    } else {
      return 1;
    }
  }
  else if (ISNA(x)) {
    return NA_REAL;
  }
  else if (ISNAN(x)) {
    return R_NaN;
  }
  else if (x == R_PosInf) {
    return R_PosInf;
  }
  else if (x == R_NegInf) {
    return 0;
  } else {
    /* can't happen */
    return NA_REAL;
  }
}

/* wrap that as an R function */

static
SEXP dexprl(SEXP x) {
  double *xa, *ra;
  SEXP res;
  int n;
  PROTECT(x = coerceVector(x, REALSXP));
  n = length(x);
  PROTECT(res = allocVector(REALSXP, n));
  xa = REAL(x);
  ra = REAL(res);
  for (;n;--n) {
    *ra = true_dexprl(*xa);
    ++xa; ++ra;
  }
  UNPROTECT(2);
  return res;
}

/* accurately compute log(1+x) / x */

static inline
double true_log1prel(const double x) {
  if (R_FINITE(x)) {
    /* actual value */
    if (x != 0) {
      return log1p(x) / x;
    } else {
      return 1;
    }
  }
  else if (ISNA(x)) {
    return NA_REAL;
  }
  else if (ISNAN(x)) {
    return R_NaN;
  }
  else if (x == R_PosInf) {
    return 0;
  }
  else if (x == R_NegInf) {
    return R_NaN;
  } else {
    /* can't happen */
    return NA_REAL;
  }
}

/* wrap that as an R function */

static
SEXP log1prel(SEXP x) {
  double *xa, *ra;
  SEXP res;
  int n;
  PROTECT(x = coerceVector(x, REALSXP));
  n = length(x);
  PROTECT(res = allocVector(REALSXP, n));
  xa = REAL(x);
  ra = REAL(res);
  for (;n;--n) {
    *ra = true_log1prel(*xa);
    ++xa; ++ra;
  }
  UNPROTECT(2);
  return res;
}


/*
   if someone could tell me why this isn't a part of base R
   I'd be insanely grateful
*/

static inline
double true_log1mexp(const double x) {
  /* accurately compute log(1-exp(x)) */
  if (R_FINITE(x)) {
    if (x > -M_LN2) {
      return log(-expm1(x));
    } else {
      return log1p(-exp(x));
    }
  }
  else if (ISNA(x)) {
    return NA_REAL;
  }
  else if (ISNAN(x)) {
    return R_NaN;
  }
  else if (x == R_PosInf) {
    return R_NaN;
  }
  else if (x == R_NegInf) {
    return 0;
  } else {
    /* can't happen */
    return NA_REAL;
  }
}

static
SEXP log1mexp(SEXP x) {
  double *xa, *ra;
  SEXP res;
  int n;
  PROTECT(x = coerceVector(x, REALSXP));
  n = length(x);
  PROTECT(res = allocVector(REALSXP, n));
  xa = REAL(x);
  ra = REAL(res);
  for (;n;--n) {
    *ra = true_log1mexp(*xa);
    ++xa; ++ra;
  }
  UNPROTECT(2);
  return res;
}

static R_CallMethodDef callMethods[] = {
  {".c.exprel", (DL_FUNC) &dexprl, 1},
  {".c.log1prel", (DL_FUNC) &log1prel, 1},
  {".c.log1mexp", (DL_FUNC) &log1mexp, 1},
  NULL
};


void R_init_texmex(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
