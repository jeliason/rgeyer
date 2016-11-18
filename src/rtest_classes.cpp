#include <Rcpp.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rtest_classes_c(NumericVector theta,
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
  for(i = 0; i < x0.nrow(); i++) x.push_back(x0(i,0), x0(i,1));
  if(dbg>10) Rprintf("ok\n");
  // point pattern initialised

  // initialise graphs
  if(dbg>10) Rprintf("Initialising Graphnested");
  int K = r.size();
  std::vector<double> rvec(K);
  for(i=0; i < K; i++) rvec.at(i) = r(i);
  Graphnested graphs(&x, rvec);
  // compute the initial graph
  if(dbg>10) Rprintf(", computing edges");
  graphs.compute_edges();
  if(dbg>10) Rprintf("ok\n");




  // the main loop
  for(it = 0; it < iter; it++ ) {

  if(dbg) Rprintf("       \r[%i/%i]", it+1, iter);
  }// eof mainloop
  if(dbg) Rprintf("\n");


  // for testing: remove one point
  i = 2;
  x.remove(&i);
  graphs.update_edges_after_delete(&i);

  // for testing: add one point
  x.push_back(0.5, 0.5);
  graphs.update_edges_after_addition();


  // for testing: move one point
  i = 0;
  x.move(&i, 0.5, 0.51);
  graphs.update_edges_after_move(&i);


  // done, return coordinates and neighbour counts

  if(dbg>10) Rprintf("building result object");
  NumericMatrix out(x.size(), 2 + K);
  for(i = 0; i < x.size(); i++){
  out(i,0) = x.getX(&i);
  out(i,1) = x.getY(&i);
  for(j = 0; j < K; j++) out(i,j+2) = graphs.neighbour_count_of_i_in_graph_k(j,i);
  }
  if(dbg>10) Rprintf("ok\n");

  return out;
}
