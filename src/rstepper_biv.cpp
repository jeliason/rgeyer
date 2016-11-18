#include <Rcpp.h>
#include <unistd.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested_biv.h"
#include "helpers.h"

using namespace Rcpp;




// [[Rcpp::export]]
NumericMatrix rstepper_biv_c(NumericVector theta,
                             NumericVector r,
                             NumericMatrix bbox,
                             int iter,
                             NumericMatrix x0,
                             int dbg,
                             int toroidal) {
  RNGScope scope;
  int i, j, k, it, mn;
  int b_or_d;
  double xn, yn, pap, change, newpot, pot;
  int acc[2] = {0,0}, n[2] = {0,0};

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

  // init point pattern object.
  NumericVector type(x0.nrow());
  if(dbg>10) Rprintf("Initialising Pplite");
  Pplite x(win, toroidal);
  for(i = 0; i < x0.nrow(); i++){
    mn = (int) x0(i,2);
    x.push_back(x0(i,0), x0(i,1), mn);
    n[mn]++;
  }
  if(dbg>10) Rprintf("ok\n");
  // point pattern initialised

  // initialise graphs
  if(dbg>10) Rprintf("Initialising Graphnested");
  int K = r.size();
  std::vector<double> rvec(K);
  for(i=0; i < K; i++) rvec.at(i) = r(i);
  Graphnested_biv   graphs(&x, rvec);
  // compute the initial graph
  if(dbg>10) Rprintf(", computing edges");
  graphs.compute_edges();
  if(dbg>10) Rprintf("ok\n");

  /* compute current potential */
  pot = 0;
  for(i=0; i < x.size(); i++){
    pot += theta[x.getType(&i)];
    for(k=0; k < K; k++) {
      pot += theta[k+2] * graphs.has_neighbours_i_in_graph_k(k, i);
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
      mn = x.getType(&j);
      // delete and compute the potential
      x.remove(&j);
      graphs.update_edges_after_delete(&j);
      newpot = 0;
      for(i=0; i < x.size(); i++){
        newpot += theta[x.getType(&i)];
        for(k=0; k < K; k++) {
          newpot += theta[k+2] * graphs.has_neighbours_i_in_graph_k(k, i);
        }
      }
      change = newpot - pot;
      pap = 0.5 * (1 + 1.0*x.size()) * exp(change) / Area ; // change here 0.5 -> prob(m==1)
      pap = fmin(1.0, pap);
      if(runif(1)(0) < pap ){ // accept
        pot = newpot;
        acc[1]++;
        n[mn]--;
      }
      else{ // reject
        // need to restore
        x.push_back(xn, yn, mn);
        graphs.update_edges_after_addition();
      }
    }
    // birth
    else{
      b_or_d = 0;
      // sample new point
      xn = runif(1, bbox(0,0),bbox(1,0))(0);
      yn = runif(1, bbox(0,1),bbox(1,1))(0);
      mn = (int) rbinom(1,1,0.5)(0);
      // add new point and compute the potential
      x.push_back(xn, yn, mn);
      graphs.update_edges_after_addition();
      newpot = 0;
      for(i=0; i < x.size(); i++){
        newpot += theta[x.getType(&i)];
        for(k=0; k < K; k++) {
            newpot += theta[k+2] * graphs.has_neighbours_i_in_graph_k(k, i);
        }
      }
      change = newpot - pot;
      pap = 2.0 * Area * exp(change) / (1.0*x.size()) ; // change here 2 -> 1/prob(m==1)
      pap = fmin(1.0, pap);
      if(runif(1)(0) < pap) { // accept
        pot = newpot;
        acc[0]++;
        n[mn]++;
      }
      else{ // reject
        j = x.size()-1;
        x.pop_back();
        graphs.update_edges_after_delete(&j);
      }
    }

    if(dbg)
      if(it % dbg_step == 0)
      {
      Rprintf("       \r[%i/%i n=(%i,%i)]", it+1, iter, n[0], n[1]);
      if(dbg>10) Rprintf("[pot%10.3f newpot %10.3f b%i d%i]", pot, newpot, acc[0],acc[1]);
      if(dbg > 1000) usleep(50000);
    }
  }// eof mainloop


  if(dbg) Rprintf("\n");

  /* here ends the simulation */

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
