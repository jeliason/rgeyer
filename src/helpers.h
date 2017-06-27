#include <Rcpp.h>
#include "Pplite.h"
#include "Graphnested.h"
#include "Graphnested_biv.h"


int sample_j(int n);
int imin(int, int);
int roundup(double);
double potential_multi(Pplite *x,
                       Rcpp::NumericVector theta0,
                       std::vector<std::vector<double > > * thetas1,
                       std::vector<std::vector<double > > * thetas2,
                       std::vector<std::vector<int > > * saturations1,
                       std::vector<std::vector<int > > * saturations2,
                       std::vector<Graphnested *> *graphs1,
                       std::vector<Graphnested *> *graphs2);

double maxv(std::vector<double> v);
