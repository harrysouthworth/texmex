#include <Rcpp.h>

#include "safe_product.h"

namespace {
  inline
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
  
  inline
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
  
  /*
    if someone could tell me why this isn't a part of base R
    I'd be insanely grateful
  */
  
  inline
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

}

//' Accurately compute (exp(x) - 1) / x
//' @param x numeric vector
//' @return numeric vector
// [[Rcpp::export(name='.exprel', rng=FALSE)]]
Rcpp::NumericVector exprel(const Rcpp::NumericVector& x) {
  return Rcpp::sapply(x, true_dexprl);
}

//' Accurately compute log(1 + x) / x
//' @param x numeric vector
//' @return numeric vector
// [[Rcpp::export(name='.log1prel', rng=FALSE)]]
Rcpp::NumericVector log1prel(const Rcpp::NumericVector& x) {
  return Rcpp::sapply(x, true_log1prel);
}

//' Accurately compute log(1-exp(x))
//' @param x numeric vector
//' @return a numeric vector
// [[Rcpp::export(name=".log1mexp", rng=FALSE)]]
Rcpp::NumericVector log1mexp(const Rcpp::NumericVector& x) {
  return Rcpp::sapply(x, true_log1mexp);
}

//' Compute pmax(x y, -1) in such a way that zeros in x beat
//' infinities in y.
//'
//' This is a common pattern in much of the distribution code, so it's
//' worth factoring out.
//' @param x a numeric vector
//' @param y a numeric vector
//' @return an appropriate numeric vector
// [[Rcpp::export(name=".specfun.safe.product", rng=FALSE)]]
Rcpp::NumericVector _safe_product(const Rcpp::NumericVector &x,
				  const Rcpp::NumericVector &y) {
  return Rcpp::mapply(x, y, safe_product);
}
