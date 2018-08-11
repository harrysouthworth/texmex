#ifndef SAFE_PRODUCT_H
#define SAFE_PRODUCT_H

#include <algorithm>

#include <Rmath.h>

namespace {
  inline double safe_product(const double x,
			     const double y) {
    const double xy = x * y;
    if ( (x==0) && ( (y == R_PosInf) || (y == R_NegInf)) ) {
      return 0;
    } else {
      return std::max(xy, -1.0);
    }
  }
}

#endif
