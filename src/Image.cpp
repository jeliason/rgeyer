#include "Rcpp.h"
#include "Image.h"

/********************************************************************************************/
Image::~Image() {
}

/********************************************************************************************/
Image::Image(List im) {
  if(!im.inherits("im")){ //  check Nill
    getValue_p = &Image::getValue_const_0;
  }
  else{
    // read components
    dim = im["dim"];
    bbox = NumericMatrix(2,2);
    NumericVector s = im["xrange"];
    bbox(_,0) = s;
    s = im["yrange"];
    bbox(_,1) = s;
    xstep = im["xstep"];
    ystep = im["ystep"];
    s = im["xcol"];
    half_xstep = s(0);
    s = im["yrow"];
    half_ystep = s(0);
    NumericMatrix v = im["v"];
    values = v;
    getValue_p = &Image::getValue_by_coordinate;
    eps = 2.220446e-16;
  }
}

double Image::getValue(double x, double y){
  return (this->*getValue_p)(x,y);
}

double Image::getValue_const_0(double x, double y) { return 0.0;}

double Image::getValue_by_coordinate(double x, double y){
  if(x + eps < bbox(0,0) | x - eps > bbox(1,0) | y + eps < bbox(0,1) | y -eps > bbox(1,1)) return 0;
  int dx =  fmax(1, fmin(round( (x - half_xstep)/xstep +1), dim(1)));
  int dy =  fmax(1, fmin(round( (y - half_ystep)/ystep +1), dim(0)));
  return values(dy-1, dx-1);
}


/* Spatstat:
nr <- im$dim[1]
nc <- im$dim[2]
cc <- round(1 + (x - im$xcol[1])/im$xstep)
rr <- round(1 + (y - im$yrow[1])/im$ystep)
cc <- pmax.int(1,pmin.int(cc, nc))
rr <- pmax.int(1,pmin.int(rr, nr))

 */
