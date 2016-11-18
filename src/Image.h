#include <Rcpp.h>
using namespace Rcpp;

#ifndef IM_C_H_
#define IM_C_H_

class Image {
  IntegerVector dim;
  double eps;
  double xstep, ystep, half_xstep, half_ystep;
  NumericMatrix bbox;
  NumericMatrix values;
  double (Image::*getValue_p)(double, double);
  public:
    Image(List);
    virtual ~Image();
    double getValue(double, double);
    double getValue_by_coordinate(double x, double y);
    double getValue_const_0(double, double);
};

#endif
