#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
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
    double tmp;
    const double my_sigma = *sigma;
    const double q = (*x - *u) / my_sigma;
    const double my_xi = *xi;
    if (my_xi != 0) {
      tmp = log1p(my_xi * q) * (-1 / my_xi - 1) - log(my_sigma);
    } else {
      /* this is the limit as xi -> 0 of the above */
      tmp = -q - log(my_sigma);
    }
    if (want_exp) {
      *res = exp(tmp);
    } else {
      *res = tmp;
    }
    ++res; ++x; ++sigma; ++xi; ++u;
  }
}

static inline
double gpd_negative_log_likelihood(const double* param,
				   const int ndata, const double* data,
				   const int nphi, const double* phi,
				   const int nxi, const double* xi,
				   double* phiwork, double* xiwork) {
  const double one = 1.0;
  const double zero = 0.0;
  const int one_i = 1;
  int itcount;
  double res = 0.0;

  /* phiwork <- X.phi .* param[:n_phi] */
  F77_CALL(dgemv)("N", &ndata, &nphi,
		  &one, phi, &ndata,
		  param, &one_i,
		  &zero, phiwork, &one_i);

  /* xiwork <- X.xi .* param[n_phi:] */
  F77_CALL(dgemv)("N", &ndata, &nxi,
		  &one, xi, &ndata,
		  param + nphi, &one_i,
		  &zero, xiwork, &one_i);

  for (itcount = ndata; itcount != 0; --itcount) {
    const double my_xi = *xiwork;
    const double my_phi = *phiwork;
    const double sigma = exp(my_phi);
    const double datum = *data;

    if (my_xi != 0) {
      const double item = my_xi * datum / sigma;
      if (item <= -1.0) {
	return HUGE_VAL;
      }
      res += my_phi + log1p(item) * (1 / my_xi + 1);
    } else {
      res += my_phi + datum / sigma;
    }
    ++data; ++phiwork; ++xiwork;
  }

  return res;

}

static
void R_gpd_negative_log_likelihood(double* result,
				   const double* param,
				   const int* ndata, const double* data,
				   const int* nphi, const double* phi,
				   const int* nxi, const double* xi,
				   double* work) {
  const int my_ndata = *ndata;
  *result = gpd_negative_log_likelihood(param,
					my_ndata, data,
				        *nphi, phi,
					*nxi, xi,
					work, work + my_ndata);
}



static R_CMethodDef cMethods[] = {
  {".c.rgpd", (DL_FUNC) &rgpd, 5},
  {".c.dgpd", (DL_FUNC) &dgpd, 7},
  {".c.gpdNegativeLogLikelihood", (DL_FUNC) &R_gpd_negative_log_likelihood, 9},
  NULL
};

void R_init_texmex(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
