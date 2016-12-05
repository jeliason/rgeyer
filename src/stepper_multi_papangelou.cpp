#include <Rcpp.h>
#include <unistd.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested.h"
#include "helpers.h"

using namespace Rcpp;


// No first order effects here, add them R-side

// [[Rcpp::export]]
NumericVector rstepper_multi_log_papangelou_c(List theta1_list,
                                              List theta2_list,
                                              NumericMatrix from,
                                              NumericMatrix to,
                                              NumericMatrix bbox,
                                              int dbg,
                                              int toroidal) {
  RNGScope scope;
  int i, j, k, l, t, t1, t2, mn;
  NumericVector log_papangelou(to.nrow());

  int ntypes = theta1_list.size();
  int npairs = ntypes * (ntypes-1) / 2;

  // Parse the parameters, change the storage form
  NumericVector theta0(ntypes);
  std::vector<std::vector<double> > ranges1(ntypes), ranges2(npairs);
  std::vector<std::vector<double> > thetas1(ntypes), thetas2(npairs);
  std::vector<std::vector<int> > saturations1(ntypes), saturations2(npairs);
  // intra-type parameters
  List lt;
  SEXP lts;
  NumericVector y;
  for(i=0; i < ntypes; i++) {
    lt = theta1_list(i);
    lts = lt["r"];
    y = lts;
    for(j = 0; j < y.size(); j++)ranges1.at(i).push_back(y(j));
    lts = lt["theta"];
    y = lts;
    for(j = 0; j < y.size(); j++)thetas1.at(i).push_back(y(j));
    lts = lt["c"];
    y = lts;
    for(j = 0; j < y.size(); j++)saturations1.at(i).push_back(y(j));
  }

  // Parse
  // the pair-wise parameters. Pairwise upper triangle i.e. 01,02,...,12,13,...,(ntypes-1)(ntypes).
  IntegerMatrix pairs(npairs, 2); // tracking which marks are in which pair.
  if(dbg>100)Rprintf("Parsing %i inter-parameters\n", npairs);
  k = 0; // over pairs
  for(t1=0; t1 < ntypes-1; t1++)
    for(t2 = t1+1; t2 < ntypes; t2++) {
      pairs(k,0) = t1;
      pairs(k,1) = t2;
      lt = theta2_list(k);
      lts = lt["r"]; // casting
      y = lts;
      for(l = 0; l < y.size(); l++) ranges2.at(k).push_back(y(l));
      lts = lt["theta"];
      y = lts;
      for(l = 0; l < y.size(); l++) thetas2.at(k).push_back(y(l));
      lts = lt["c"];
      y = lts;
      for(l = 0; l < y.size(); l++)saturations2.at(k).push_back(y(l));
      k++;
    }

    // init point pattern window here so don't need dependence on Rcpp in Pplite
    std::vector<double> win(bbox.ncol()*2);
  double Area=1;
  for(i=0; i < bbox.ncol(); i++) {
    win.at(2*i) = bbox(0, i);
    win.at(2*i+1) = bbox(1, i);
    Area *= bbox(1,i) - bbox(0,i);
  }

  // init data point pattern object.
  NumericVector type(from.nrow());

  if(dbg>10) Rprintf("Initialising Pplite");
  Pplite x(win, toroidal);
  for(i = 0; i < from.nrow(); i++){
    mn = (int) from(i,2);
    x.push_back(from(i,0), from(i,1), mn);
  }


  // initialise graphs. ntypes + npairs of them.
  std::vector<Graphnested *> graphs1(ntypes);
  std::vector<Graphnested *> graphs2(npairs);

  // first the intra-type graphs
  if(dbg>10) Rprintf("Initialising Graphnested, intra (%i)", graphs1.size());
  std::vector<int> target(1);
  for(t = 0; t < ntypes; t++) {
    graphs1.at(t) = new Graphnested(&x, ranges1.at(t));
    target.at(0) = t; // note that types are 0,..., ntypes-1 integers
    graphs1.at(t)->set_pars_i(target);
    graphs1.at(t)->set_pair_suitable_method(&Graphnested::pair_suitable_same_mark_but_of_given_set);
    graphs1.at(t)->compute_edges();
  }

  // then the inter-type graphs
  if(dbg>10) Rprintf("\nInitialising Graphnested, inter (%i)", graphs2.size());
  target.resize(2);
  for(t = 0; t < npairs; t++) {
    graphs2.at(t) = new Graphnested(&x, ranges2.at(t));
    target.at(0) = pairs(t,0);
    target.at(1) = pairs(t,1);
    graphs2.at(t)->set_pars_i(target);
    graphs2.at(t)->set_pair_suitable_method(&Graphnested::pair_suitable_different_mark_but_of_given_set);
    graphs2.at(t)->compute_edges();
  }

  // the potential
  double pot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);
  double newpot;

  ////////
  // Main loop

  double xn, yn;
  int nm;

  for(i = 0; i < to.nrow(); i++) {
    xn = to.at(i,0);
    yn = to.at(i,1);
    mn = (int) to.at(i,2);
    // add the new location
    x.push_back(xn, yn, mn);
    // recompute the graphs:
    for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_addition();
    for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_addition();
    // potential
    newpot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);
    // the difference is the
    log_papangelou.at(i) = newpot - pot;
    // restore
    j = x.size() - 1;
    x.pop_back();
    for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_delete(&j);
    for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_delete(&j);
    if(dbg) {
      Rprintf("       \r[%i/%i]", i+1, to.nrow());
    }
  }
  if(dbg) Rprintf("\n");
  return log_papangelou;
}




// [[Rcpp::export]]
NumericVector rstepper_multi_log_papangelou_at_data_c(List theta1_list,
                                              List theta2_list,
                                              NumericMatrix from,
                                              NumericMatrix bbox,
                                              int dbg,
                                              int toroidal) {
  RNGScope scope;
  int i, j, k, l, t, t1, t2, nm;
  NumericVector log_papangelou(from.nrow());

  int ntypes = theta1_list.size();
  int npairs = ntypes * (ntypes-1) / 2;

  // Parse the parameters, change the storage form
  NumericVector theta0(ntypes);
  std::vector<std::vector<double> > ranges1(ntypes), ranges2(npairs);
  std::vector<std::vector<double> > thetas1(ntypes), thetas2(npairs);
  std::vector<std::vector<int> > saturations1(ntypes), saturations2(npairs);
  // intra-type parameters
  List lt;
  SEXP lts;
  NumericVector y;
  for(i=0; i < ntypes; i++) {
    lt = theta1_list(i);
    lts = lt["r"];
    y = lts;
    for(j = 0; j < y.size(); j++)ranges1.at(i).push_back(y(j));
    lts = lt["theta"];
    y = lts;
    for(j = 0; j < y.size(); j++)thetas1.at(i).push_back(y(j));
    lts = lt["c"];
    y = lts;
    for(j = 0; j < y.size(); j++)saturations1.at(i).push_back(y(j));
  }

  // Parse
  // the pair-wise parameters. Pairwise upper triangle i.e. 01,02,...,12,13,...,(ntypes-1)(ntypes).
  IntegerMatrix pairs(npairs, 2); // tracking which marks are in which pair.
  if(dbg>100)Rprintf("Parsing %i inter-parameters\n", npairs);
  k = 0; // over pairs
  for(t1=0; t1 < ntypes-1; t1++)
    for(t2 = t1+1; t2 < ntypes; t2++) {
      pairs(k,0) = t1;
      pairs(k,1) = t2;
      lt = theta2_list(k);
      lts = lt["r"]; // casting
      y = lts;
      for(l = 0; l < y.size(); l++) ranges2.at(k).push_back(y(l));
      lts = lt["theta"];
      y = lts;
      for(l = 0; l < y.size(); l++) thetas2.at(k).push_back(y(l));
      lts = lt["c"];
      y = lts;
      for(l = 0; l < y.size(); l++)saturations2.at(k).push_back(y(l));
      k++;
    }

    // init point pattern window here so don't need dependence on Rcpp in Pplite
    std::vector<double> win(bbox.ncol()*2);
  double Area=1;
  for(i=0; i < bbox.ncol(); i++) {
    win.at(2*i) = bbox(0, i);
    win.at(2*i+1) = bbox(1, i);
    Area *= bbox(1,i) - bbox(0,i);
  }

  // init data point pattern object.
  NumericVector type(from.nrow());

  if(dbg>10) Rprintf("Initialising Pplite");
  Pplite x(win, toroidal);
  for(i = 0; i < from.nrow(); i++){
    nm = (int) from(i,2);
    x.push_back(from(i,0), from(i,1), nm);
  }



  // initialise graphs. ntypes + npairs of them.
  std::vector<Graphnested *> graphs1(ntypes);
  std::vector<Graphnested *> graphs2(npairs);

  // first the intra-type graphs
  if(dbg>10) Rprintf("Initialising Graphnested, intra (%i)", graphs1.size());
  std::vector<int> target(1);
  for(t = 0; t < ntypes; t++) {
    graphs1.at(t) = new Graphnested(&x, ranges1.at(t));
    target.at(0) = t; // note that types are 0,..., ntypes-1 integers
    graphs1.at(t)->set_pars_i(target);
    graphs1.at(t)->set_pair_suitable_method(&Graphnested::pair_suitable_same_mark_but_of_given_set);
    graphs1.at(t)->compute_edges();
  }

  // then the inter-type graphs
  if(dbg>10) Rprintf("\nInitialising Graphnested, inter (%i)", graphs2.size());
  target.resize(2);
  for(t = 0; t < npairs; t++) {
    graphs2.at(t) = new Graphnested(&x, ranges2.at(t));
    target.at(0) = pairs(t,0);
    target.at(1) = pairs(t,1);
    graphs2.at(t)->set_pars_i(target);
    graphs2.at(t)->set_pair_suitable_method(&Graphnested::pair_suitable_different_mark_but_of_given_set);
    graphs2.at(t)->compute_edges();
  }

  // the potential
  if(dbg>10) Rprintf("\nInitial potential");
  double pot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);
  double newpot;
  if(dbg>10) Rprintf(" ok\n");
  ////////
  // Main loop

  double xn, yn, mn;
  j = 0;
  // stack, always remove first and return to back.
  for(i = 0; i < from.nrow(); i++) {
    xn = x.getX(&j);
    yn = x.getY(&j);
    mn = x.getType(&j);
    // // delete
    x.remove(&j);
    // recompute the graphs:
    for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_delete(&j);
    for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_delete(&j);
    // potential
    newpot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);
    // and f(X u x) / f(X)
    log_papangelou.at(i) = pot - newpot;
    // put the point back
    x.push_back(xn, yn, (int) mn);
    for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_addition();
    for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_addition();

    if(dbg) {
      Rprintf("       \r[%i/%i]", i+1, from.nrow());
    }
  }
  if(dbg) Rprintf("\n");
  return log_papangelou;
}



