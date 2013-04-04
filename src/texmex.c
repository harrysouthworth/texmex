#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* herewith sampling from the GPD distribution */

static
void rgpd(double* res,
	  const int* num,
	  const double* sigma,
	  const double* xi,
	  const double* u) {
  int it;
  GetRNGstate();
  for (it=*num; it!=0; --it) {
    if (*xi != 0) {
      *res = *u + (*sigma / *xi) * (R_pow(unif_rand(), -*xi) - 1);
    } else {
      *res = *u + exp_rand() * (*sigma);
    }
    ++res;
    ++sigma;
    ++xi;
    ++u;
  }
  PutRNGstate();
}

static
void dgpd(double* res,
	  const int* num,
	  const double *x,
	  const double *sigma,
	  const double *xi,
	  const double* u,
	  const int* log_res) {
  int it;
  const int want_exp = !*log_res;
  for (it=*num; it!=0; --it) {
    const double my_x = *x;
    const double my_sigma = *sigma;
    const double my_u = *u;
    const double my_xi = *xi;
    if (!(ISNA(my_x) || ISNA(my_sigma) || ISNA(my_u) || ISNA(my_xi))) {
      const double q = (my_x - my_u) / my_sigma;
      double tmp;
      if (my_xi != 0) {
	const double work = my_xi * q;
	if (work <= -1) {
	  tmp = R_NegInf;
	} else {
	  tmp = log1p(work) * (-1 / my_xi - 1) - log(my_sigma);
	}
      } else {
	/* this is the limit as xi -> 0 of the above */
	tmp = -q - log(my_sigma);
      }
      if (want_exp) {
	*res = exp(tmp);
      } else {
	*res = tmp;
      }
    } else {
      *res = NA_REAL;
    }
    ++res; ++x; ++sigma; ++xi; ++u;
  }
}

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


static R_CMethodDef cMethods[] = {
  {".c.rgpd", (DL_FUNC) &rgpd, 5},
  {".c.dgpd", (DL_FUNC) &dgpd, 7},
  NULL
};

static R_CallMethodDef callMethods[] = {
  {".c.dexprl", (DL_FUNC) &dexprl, 1},
  NULL
};

void R_init_texmex(DllInfo *info) {
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}
