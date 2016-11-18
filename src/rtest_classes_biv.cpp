#include <Rcpp.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested_biv.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rtest_classes_biv_c(NumericVector theta,
                                  NumericVector r,
                                  NumericMatrix bbox,
                                  int iter,
                                  NumericMatrix x0,
                                  int dbg,
                                  int toroidal) {
  int i, j, it;

  // init point pattern window here so don't need dependence on Rcpp in Pplite
  std::vector<double> win(bbox.ncol()*2);
  for(i=0; i < bbox.ncol(); i++) {
  win.at(2*i) = bbox(0, i);
  win.at(2*i+1) = bbox(1, i);
  }

  // init point pattern object
  if(dbg>10) Rprintf("Initialising Pplite");
  Pplite x(win, toroidal);
  for(i = 0; i < x0.nrow(); i++)   x.push_back(x0(i,0), x0(i,1), (int) x0(i,2));
  if(dbg>10) Rprintf("ok\n");
  // point pattern initialised

  // initialise graphs
  if(dbg>10) Rprintf("Initialising Graphnested");
  int K = r.size();
  std::vector<double> rvec(K);
  for(i=0; i < K; i++) rvec.at(i) = r(i);
  Graphnested_biv graphs(&x, rvec);
  // compute the initial graph
  if(dbg>10) Rprintf(", computing edges");
  graphs.compute_edges();
  if(dbg>10) Rprintf("ok\n");

  // for testing: remove one point
  if(dbg>10) Rprintf("removing second item");
  i = 1;
  x.remove(&i);
  graphs.update_edges_after_delete(&i);
  if(dbg>10) Rprintf("ok\n");

  // for testing: add one point
  if(dbg>10) Rprintf("adding a new point");
  x.push_back(0.5, 0.5, 0);
  graphs.update_edges_after_addition();
  if(dbg>10) Rprintf("ok\n");


  // for testing: move one point
  if(dbg>10) Rprintf("moving first point");
  i = 0;
  x.move(&i, 0.5, 0.51);
  graphs.update_edges_after_move(&i);
  if(dbg>10) Rprintf("ok\n");


  // done, return coordinates and neighbour counts

  if(dbg>10) Rprintf("building result object");
  NumericMatrix out(x.size(), 3 + K);
  for(i = 0; i < x.size(); i++){
  out(i,0) = x.getX(&i);
  out(i,1) = x.getY(&i);
  out(i,2) = x.getType(&i);
  for(j = 0; j < K; j++) out(i,j+3) = graphs.neighbour_count_of_i_in_graph_k(j,i);
  }
  if(dbg>10) Rprintf("ok\n");

  return out;
}
