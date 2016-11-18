#include <Rcpp.h>
#include <unistd.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested.h"
#include "helpers.h"

using namespace Rcpp;


// For non-data points

// [[Rcpp::export]]
NumericMatrix rstepper_components_c(
                                    NumericMatrix from,
                                    NumericMatrix to,
                                    NumericVector r,
                                    NumericMatrix bbox,
                                    int dbg,
                                    int toroidal) {

  RNGScope scope;
  int i, j, l, k;
  double xn, yn,  newpot, pot, change;

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


  NumericMatrix components(to.nrow(), K);

  NumericMatrix original_factors(from.nrow(), K);

  for(j = 0; j < x.size(); j++)
    for(k = 0; k < K; k++)
      original_factors(j, k) = graphs.has_neighbours_i_in_graph_k(k, j);

  /* here begins the evaluation */

  // the main loop
  for(i = 0; i < to.nrow(); i++ ) {
    // new location
    xn = to(i,0);
    yn = to(i,1);
    // add new point and compute the unweighted potential difference
    x.push_back(xn, yn);
    l = x.size()-1;
    graphs.update_edges_after_addition();

    for(k = 0; k < K; k++) {
      change = graphs.has_neighbours_i_in_graph_k(k, l); // the new point
      for(j = 0; j < l; j++) { // and others
        change += graphs.has_neighbours_i_in_graph_k(k, j) - original_factors(j, k);
      }
      components(i,k) = change;
    }

    // restore
    x.pop_back();
    graphs.update_edges_after_delete(&l);
    if(dbg) {
      Rprintf("       \r[%i/%i]", i+1, to.nrow(), x.size());
    }
  }// eof mainloop


  if(dbg) Rprintf("\n");
  return components;
}




/// For data points

// [[Rcpp::export]]
NumericMatrix rstepper_components_at_data_c(
    NumericMatrix from,
    NumericVector r,
    NumericMatrix bbox,
    int dbg,
    int toroidal) {

  RNGScope scope;
  int i, j, l, k;
  double xn, yn, change;

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
  graphs.compute_edges(); // this gives the upper part of the ratio exp(x cup X)/exp(X)
  if(dbg>10) Rprintf("ok\n");


  NumericMatrix components(from.nrow(), K);

  NumericMatrix original_factors(from.nrow(), K);

  for(j = 0; j < x.size(); j++)
    for(k = 0; k < K; k++)
      original_factors(j, k) = graphs.has_neighbours_i_in_graph_k(k, j);

  /* here begins the evaluation */

  // the main loop
  for(i = 0; i < from.nrow(); i++ ) {
    for(k = 0; k < K; k++) {
      change = -graphs.has_neighbours_i_in_graph_k(k, i); // the old point
      for(j = 0; j < graphs.neighbour_count_of_i_in_graph_k(k, i); j++) { // check if the point is the only neighbour
        if(graphs.neighbour_count_of_neighbour_j_of_i_in_graph_k(k, i, j) == 1) change--;
      }
      components(i,k) = change;
    }

    if(dbg) {
      Rprintf("       \r[%i/%i]", i+1, x.size());
    }
  }// eof mainloop


  if(dbg) Rprintf("\n");
  return components;
}

