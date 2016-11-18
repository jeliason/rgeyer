#include <Rcpp.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rtest_graph_c(NumericMatrix x0, NumericMatrix bbox, NumericVector r) {
  int i, j, it;

  // init point pattern window here so don't need dependence on Rcpp in Pplite
  std::vector<double> win(bbox.ncol()*2);
  for(i=0; i < bbox.ncol(); i++) {
  win.at(2*i) = bbox(0, i);
  win.at(2*i+1) = bbox(1, i);
  }
  int dbg = 100;
  int toroidal = 1;
  // init point pattern object
  if(dbg>10) Rprintf("Initialising Pplite");
  Pplite x(win, toroidal);
  for(i = 0; i < x0.nrow(); i++) x.push_back(x0(i,0), x0(i,1), (int) x0(i,2));
  if(dbg>10) Rprintf("ok\n");
  // point pattern initialised

  std::vector<double > rvec(r.size());
  for(i=0; i < r.size(); i++) rvec.at(i) = r(i);

  // Now lets check the different graph options:

  // Within types
  NumericVector types = unique(x0(_,2));
  types.sort();
  int ntypes = types.length();

  Rprintf("Found %i types\n", ntypes);

  std::vector<Graphnested *> graphs1(ntypes);
  std::vector<int > target(1);
  if(dbg>10) Rprintf("Intra-graphs");

  for(i = 0; i < ntypes; i++) {
    graphs1.at(i) = new Graphnested(&x, rvec);
    target.at(0) = types(i);
    graphs1.at(i)->set_pars_i(target);
    graphs1.at(i)->set_pair_suitable_method(&Graphnested::pair_suitable_same_mark_but_of_given_set);
    // compute
    graphs1.at(i)->compute_edges();
  }
  if(dbg>10) Rprintf(" ok\n");

  // betweentypes
  // Within types
  int npairs = ntypes * (ntypes-1) / 2;
  std::vector<Graphnested *> graphs2(npairs);
  target.resize(2);
  if(dbg>10) Rprintf("Inter");
  int k = 0;
  for(i = 0; i < ntypes-1; i++) {
    target.at(0) = types(i);
    for(j = i+1 ; j < ntypes; j++){
      graphs2.at(k) = new Graphnested(&x, rvec);
      target.at(1) = types(j);
      graphs2.at(k)->set_pars_i(target);
      graphs2.at(k)->set_pair_suitable_method(&Graphnested::pair_suitable_different_mark_but_of_given_set);
      graphs2.at(k)->compute_edges();
      k++;
    }
  }
  if(dbg>10) Rprintf(" ok\n");

  // done, return coordinates and neighbour counts per graph

  if(dbg>10) Rprintf("building result object");
  List out1(ntypes);
  NumericMatrix *M;
  int l;
  for(l = 0; l < ntypes; l++){
    M = new NumericMatrix(x.size(), r.size());
    for(k = 0; k < r.size(); k++){
      for(i = 0; i  < x.size(); i++) {
        M->at(i,k) = graphs1.at(l)->neighbour_count_of_i_in_graph_k(k,i);
      }
    }
    out1(l) = *M;
  }
  List out2(npairs);
  for(l = 0; l < npairs; l++){
    M = new NumericMatrix(x.size(), r.size());
    for(k = 0; k < r.size(); k++){
      for(i = 0; i  < x.size(); i++) {
        M->at(i,k) = graphs2.at(l)->neighbour_count_of_i_in_graph_k(k,i);
      }
    }
    out2(l) = *M;
  }

  if(dbg>10) Rprintf(" ok\n");

  return List::create(out1, out2);
}
