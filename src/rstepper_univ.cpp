#include <Rcpp.h>
#include <unistd.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested.h"
#include "helpers.h"
#include "Image.h"

using namespace Rcpp;




// [[Rcpp::export]]
NumericMatrix rstepper_univ_c(NumericVector theta,
                              NumericVector r,
                              NumericMatrix bbox,
                              int iter,
                              NumericMatrix x0,
                              int dbg,
                              int toroidal,
                              List trend_im) {
  RNGScope scope;
  int i, j, k, it;
  int b_or_d;
  double xn, yn, pap, change, newpot, pot;
  int acc[2] = {0,0};
  int dbg_step = (int) (1.0*iter/100.0);
  if(dbg_step <1) dbg_step = 1;
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

  // init trend
  Image trend(trend_im);


  /* compute current potential */
  pot = 0;
  for(i=0; i < x.size(); i++){
    pot += theta[0] + trend.getValue(x.getX(&i), x.getY(&i));
    for(k=0; k < K; k++) {
      pot += theta[k+1] * graphs.has_neighbours_i_in_graph_k(k, i);
    }
  }

  /* here begins the simulation */

  // the main loop
  for(it = 0; it < iter; it++ ) {
    // flip a coin, birth or death
    if(runif(1)(0) < 0.5) { // death
      b_or_d = 1;
      // sample the dying point
      j = sample_j(x.size());
      // keep the old location
      xn = x.getX(&j);
      yn = x.getY(&j);
      // delete and compute the potential
      x.remove(&j);
      graphs.update_edges_after_delete(&j);
      newpot = 0;
      for(i=0; i < x.size(); i++){
        pot += theta[0] + trend.getValue(x.getX(&i), x.getY(&i));
        for(k=0; k < K; k++) {
          newpot += theta[k+1] * graphs.has_neighbours_i_in_graph_k(k, i);
        }
      }
      change = newpot - pot;
      pap = (1 + 1.0*x.size()) * exp(change) / Area ;
      pap = fmin(1.0, pap);
      if(runif(1)(0) < pap ){ // accept
        pot = newpot;
        acc[1]++;
      }
      else{ // reject
        // need to restore
        x.push_back(xn, yn);
        graphs.update_edges_after_addition();
      }
    }
    // birth
    else{
      b_or_d = 0;
      // sample new point
      xn = runif(1, bbox(0,0),bbox(1,0))(0);
      yn = runif(1, bbox(0,1),bbox(1,1))(0);
      // add new point and compute the potential
      x.push_back(xn, yn);
      graphs.update_edges_after_addition();
      newpot = 0;
      for(i=0; i < x.size(); i++){
        newpot += theta[0] + trend.getValue(x.getX(&i), x.getY(&i));
        for(k=0; k < K; k++) {
          newpot += theta[k+1] * graphs.has_neighbours_i_in_graph_k(k, i);
        }
      }
      change = newpot - pot;
      pap = Area * exp(change) / (1.0*x.size()) ;
      pap = fmin(1.0, pap);
      if(runif(1)(0) < pap) { // accept
        pot = newpot;
        acc[0]++;
      }
      else{ // reject
        j = x.size()-1;
        x.pop_back();
        graphs.update_edges_after_delete(&j);
      }
    }

    if(dbg) if(it % dbg_step == 0){
      Rprintf("       \r[%i/%i n=%i]", it+1, iter, x.size());
      if(dbg>10) Rprintf("[pot%10.3f newpot %10.3f b%i d%i]", pot, newpot, acc[0],acc[1]);
      if(dbg > 1000) usleep(100);
    }
  }// eof mainloop


  if(dbg) Rprintf("\n");

  /* here ends the simulation */

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






// [[Rcpp::export]]
NumericMatrix rstepper_univ_fixed_c(NumericVector n,
                                    NumericVector theta,
                                    NumericVector r,
                                    NumericMatrix bbox,
                                    int iter,
                                    NumericMatrix x0,
                                    int dbg,
                                    int toroidal,
                                    List trend_im) {
  RNGScope scope;
  int i, j, k, it;
  int b_or_d;
  double xn, yn, xo, yo, pap, change, newpot, pot;
  int acc[2] = {0,0};
  int dbg_step = (int) (1.0*iter/100.0);
  if(dbg_step <1) dbg_step = 1;
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

  if(dbg>10){
    Rprintf("Parameters:\n * Ranges\n");
    for(i=0; i < K; i++) {
      Rprintf("%f ", rvec.at(i));
    }
    Rprintf("\n * Thetas\n");
    for(i=0; i < K; i++) {
      Rprintf("%f ", theta[i]);
    }
    Rprintf("\n");
  }

  // trend
  Image trend(trend_im);

  /* compute current potential */
  pot = 0;
  for(i=0; i < x.size(); i++){
    for(k=0; k < K; k++) {
      pot += theta[k] * graphs.has_neighbours_i_in_graph_k(k, i); // no intercepts here
    }
  }

  /* here begins the simulation */

  // the main loop
  for(it = 0; it < iter; it++ ) {
    // sample the moving point
    j = sample_j(x.size());
    // keep the old location
    xo = x.getX(&j);
    yo = x.getY(&j);
    // sample new location
    xn = runif(1, bbox(0,0),bbox(1,0))(0);
    yn = runif(1, bbox(0,1),bbox(1,1))(0);
    // move
    x.move(&j, xn, yn);
    // delete and compute the potential
    graphs.update_edges_after_move(&j);
    newpot = 0;
    for(i=0; i < x.size(); i++){
      for(k=0; k < K; k++) {
        newpot += theta[k] * graphs.has_neighbours_i_in_graph_k(k, i);
      }
    }
    change = newpot - pot + (trend.getValue(xn, yn) - trend.getValue(xo,yo));
    pap = exp(change);
    pap = fmin(1.0, pap);
    if(runif(1)(0) < pap ){ // accept move
      pot = newpot;
      acc[0]++;
    }
    else{ // reject
      // need to restore
      x.move(&j, xo, yo);
      graphs.update_edges_after_move(&j); // could try to remember these?
    }

    if(dbg)
      if(it % dbg_step == 0){
        Rprintf("       \r[%i/%i acc %i %4.2f%%]", it, iter, acc[0], 100.0*acc[0]/it);
        if(dbg>10) Rprintf("[pot%10.3f newpot %10.3f]", pot, newpot);
        if(dbg > 1000) usleep(50000);
    }
  }// eof mainloop


  if(dbg) Rprintf("       \r[%i/%i acc %i %4.2f%%]\n", it, iter, acc[0], 100.0*acc[0]/it);
  /* here endeth ye simulation */

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
