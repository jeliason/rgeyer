#include <Rcpp.h>
#include <unistd.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested.h"
#include "helpers.h"

using namespace Rcpp;




// [[Rcpp::export]]
NumericVector rstepper_log_papangelou_c(
                                    NumericMatrix from,
                                    NumericMatrix to,
                                    NumericVector theta,
                                    NumericVector r,
                                    NumericVector sat,
                                    NumericMatrix bbox,
                                    int dbg,
                                    int toroidal) {

  RNGScope scope;
  int i, j, k;
  double xn, yn,  newpot, pot, change;
  NumericVector log_papangelou(to.nrow());

  // init point pattern window here so don't need dependence on Rcpp in Pplite
  std::vector<double> win(bbox.ncol()*2);
  double Area=1;
  for(i=0; i < bbox.ncol(); i++) {
    win.at(2*i) = bbox(0, i);
    win.at(2*i+1) = bbox(1, i);
    Area *= bbox(1,i) - bbox(0,i);
  }

  // init point pattern object
  if(dbg>10) Rprintf("Initialising Pplite");
  Pplite x(win, toroidal);
  for(i = 0; i < from.nrow(); i++) x.push_back(from(i,0), from(i,1));
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

  /* compute current potential */
  pot = 0;
  for(i=0; i < x.size(); i++){
    pot += theta[0];
    for(k=0; k < K; k++) {
      pot += theta[k+1] * imin(sat(k), graphs.neighbour_count_of_i_in_graph_k(k,i));
    }
  }

  /* here begins the evaluation */

  // the main loop
  for(i = 0; i < to.nrow(); i++ ) {
    // new location
    xn = to(i,0);
    yn = to(i,1);
    // add new point and compute the potential
    x.push_back(xn, yn);
    graphs.update_edges_after_addition();
    newpot = 0;
    for(j=0; j < x.size(); j++){
      newpot += theta[0];
      for(k=0; k < K; k++) {
        newpot += theta[k+1] * imin(sat(k), graphs.neighbour_count_of_i_in_graph_k(k,j));
      }
    }
    change = newpot - pot;
    log_papangelou(i) = change;
    // restore
    j = x.size()-1;
    x.pop_back();
    graphs.update_edges_after_delete(&j);
    if(dbg) {
      Rprintf("       \r[%i/%i]", i+1, to.nrow(), x.size());
    }
  }// eof mainloop


  if(dbg) Rprintf("\n");
  return log_papangelou;
}
