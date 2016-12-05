#include <Rcpp.h>
#include <unistd.h>
#include <vector>
#include "Pplite.h"
#include "Graphnested.h"
#include "Graphnested_biv.h"
#include "helpers.h"
#include "Image.h"
using namespace Rcpp;

void print_input_parameters(std::vector<std::vector<double > > * ranges1,
                            std::vector<std::vector<double > > * ranges2,
                            std::vector<std::vector<double > > * thetas1,
                            std::vector<std::vector<double > > * thetas2,
                            std::vector<std::vector<int > > * saturations1,
                            std::vector<std::vector<int > > * saturations2){
  int ntypes = ranges1->size();
  int i,j;
  int npairs = ntypes*(ntypes-1)/2;
  int k = 0;
  IntegerMatrix pairs(npairs,2);
  for(i=0; i < ntypes-1; i++)
    for(j = i+1; j < ntypes; j++) {
      pairs(k,0) = i;
      pairs(k,1) = j;
      k++;
  }

  Rprintf("Here's what we got: \n ** Ranges, Intra:\n");
  for(i=0; i < ntypes; i++){
    Rprintf("%i: ", i);
    for(j = 0; j < ranges1->at(i).size(); j ++){
      Rprintf("%f ", ranges1->at(i).at(j));
    }
    Rprintf("\n");
  }
  Rprintf("\n ** Thetas, Intra:\n");
  for(i=0; i < ntypes; i++){
    Rprintf("%i: ", i);
    for(j = 0; j < thetas1->at(i).size(); j ++){
      Rprintf("%f ", thetas1->at(i).at(j));
    }
    Rprintf("\n");
  }
  Rprintf("\n ** Saturations, Intra:\n");
  for(i=0; i < ntypes; i++){
    Rprintf("%i: ", i);
    for(j = 0; j < saturations1->at(i).size(); j ++){
      Rprintf("%i ", saturations1->at(i).at(j));
    }
    Rprintf("\n");
  }

  Rprintf("\n ** Ranges, Inter:\n");
  for(i=0; i < npairs; i++){
    Rprintf("%i-%i: ", pairs(i,0), pairs(i,1));
    for(j = 0; j < ranges2->at(i).size(); j ++){
      Rprintf("%f ", ranges2->at(i).at(j));
    }
    Rprintf("\n");
  }
  Rprintf("\n ** Thetas, Inter:\n");
  for(i=0; i < npairs; i++){
    Rprintf("%i-%i: ", pairs(i,0), pairs(i,1));
    for(j = 0; j < thetas2->at(i).size(); j ++){
      Rprintf("%f ", thetas2->at(i).at(j));
    }
    Rprintf("\n");
  }
  Rprintf("\n ** Saturations, Inter:\n");
  for(i=0; i < npairs; i++){
    Rprintf("%i-%i: ", pairs(i,0), pairs(i,1));
    for(j = 0; j < saturations2->at(i).size(); j ++){
      Rprintf("%f ", saturations2->at(i).at(j));
    }
    Rprintf("\n");
  }
}


// double potential_multi(Pplite *x, NumericVector theta0,
//                        std::vector<std::vector<double > > * thetas1,
//                        std::vector<std::vector<double > > * thetas2,
//                        std::vector<std::vector<int > > * saturations1,
//                        std::vector<std::vector<int > > * saturations2,
//                        std::vector<Graphnested *> *graphs1,
//                        std::vector<Graphnested *> *graphs2
// ){
//   int i, k, t;
//   double pot = 0;
//   for(i=0; i < x->size(); i++){
//     pot += theta0[x->getType(&i)];
//     // intra
//     for(t=0; t < thetas1->size(); t++) {
//       // over ranges
//       for(k = 0; k < thetas1->at(t).size(); k++)
//         //pot += thetas1->at(t).at(k) * graphs1->at(t)->has_neighbours_i_in_graph_k(k, i);
//         pot += thetas1->at(t).at(k) * imin(graphs1->at(t)->neighbour_count_of_i_in_graph_k(k,i), saturations1->at(t).at(k));
//     }
//     // inter
//     for(t=0; t < thetas2->size(); t++) {
//       // over ranges
//       for(k = 0; k < thetas2->at(t).size(); k++)
//         //pot += thetas2->at(t).at(k) * graphs2->at(t)->has_neighbours_i_in_graph_k(k, i);
//         pot += thetas2->at(t).at(k) * imin(graphs2->at(t)->neighbour_count_of_i_in_graph_k(k,i), saturations2->at(t).at(k));
//     }
//   }
//   return pot;
// }


// [[Rcpp::export]]
NumericMatrix rstepper_multi_c(NumericVector theta0,
                               List theta1_list,
                               List theta2_list,
                               NumericMatrix bbox,
                               int iter,
                               NumericMatrix x0,
                               int dbg,
                               int toroidal,
                               List trend_ims) {
  RNGScope scope;
  int i, j, l, k, it, mn, t, t1, t2;
  int b_or_d;
  double xn, yn, pap, change, newpot, pot;
  int acc[2] = {0,0};

  int dbg_step = (int) (1.0*iter/100.0);
  if(dbg_step <1) dbg_step = 1;

  int ntypes = theta0.size();
  int npairs = ntypes * (ntypes-1) / 2;

  // Parse the parameters, change the storage form
  std::vector<std::vector<double> > ranges1(ntypes), ranges2(npairs);
  std::vector<std::vector<double> > thetas1(ntypes), thetas2(npairs);
  std::vector<std::vector<int> > saturations1(ntypes), saturations2(npairs);
  if(dbg>100)Rprintf("Parsing %i intra-parameters\n", ntypes);
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

    // gotta check it worked


    if(dbg>100){
      print_input_parameters(&ranges1,&ranges2, &thetas1,&thetas2, &saturations1, &saturations2);
    }

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
    }
    if(dbg>10) Rprintf("ok\n");
    // point pattern initialised


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

    // trends
    std::vector<Image *> trends(trend_ims.size());
    for(t = 0; t < ntypes; t++) {
      List f = trend_ims(t);
      trends.at(t) = new Image(f);
    }
    if(dbg>10) Rprintf("\nInitial potential");

    /* compute current potential */

    pot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);

    /* here begins the simulation */


    if(dbg > 10) Rprintf("\nStarting simulation\n");

    // the main loop
    for(it = 0; it < iter; it++ ) {
      Rcpp::checkUserInterrupt();
      // flip a coin, birth or death
      if(runif(1)(0) < 0.5) { // death
        b_or_d = 1;
        // sample the dying point
        j = sample_j(x.size());
        // keep the old location
        xn = x.getX(&j);
        yn = x.getY(&j);
        mn = x.getType(&j);
        // delete
        x.remove(&j);
        // update graphs
        for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_delete(&j);
        for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_delete(&j);
        // new potential
        newpot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);
        // change and MH-ratio
        change = newpot - pot;
        // add trend effects
        change -= trends.at(mn)->getValue(xn,yn);
        pap = (1 + 1.0*x.size()) * exp(change) / (ntypes * Area); // change here 1/ntypes -> prob(m==1)
        pap = fmin(1.0, pap);
        if(runif(1)(0) < pap ){ // accept the deletion
          // Rprintf("d-a\n");
          pot = newpot;
          acc[1]++;
        }
        else{ // reject
          // Rprintf("b-a\n");
          // need to recompute
          x.push_back(xn, yn, mn);
          for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_addition();
          for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_addition();
        }
      }
      // birth
      else{
        b_or_d = 0;
        // sample new point
        xn = runif(1, bbox(0,0),bbox(1,0))(0);
        yn = runif(1, bbox(0,1),bbox(1,1))(0);
        mn = (int) sample_j(ntypes);
        // add new point and compute the potential
        x.push_back(xn, yn, mn);
        for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_addition();
        for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_addition();
        newpot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);
        change = newpot - pot;
        // trend effect
        change += trends.at(mn)->getValue(xn,yn);
        pap = ntypes * Area * exp(change) / (1.0 * x.size()) ; // change here ntypes -> 1/prob(m==1)
        pap = fmin(1.0, pap);
        if(runif(1)(0) < pap) { // accept addition
          // Rprintf("b-a\n");
          pot = newpot;
          acc[0]++;
        }
        else{ // reject
          // Rprintf("b-r\n");
          j = x.size()-1;
          x.pop_back();
          for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_delete(&j);
          for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_delete(&j);
        }
      }

      if(dbg)
        if(it % dbg_step == 0)
        {
          Rprintf("       \r[%i/%i n=%i]", it+1, iter, x.size());
          if(dbg>10) Rprintf("[pot%10.3f newpot %10.3f b%i d%i]", pot, newpot, acc[0],acc[1]);
          if(dbg > 1000) usleep(50000);
        }
    }// eof mainloop

    if(dbg) Rprintf("\n");

    /* here ends the simulation */

    // done, return coordinates and neighbour counts
    if(dbg>10) Rprintf("building result object");
    NumericMatrix out(x.size(), 3);
    for(i = 0; i < x.size(); i++){
      out(i,0) = x.getX(&i);
      out(i,1) = x.getY(&i);
      out(i,2) = x.getType(&i);
    }
    if(dbg>10) Rprintf(" ok\n");

    return out;
}


/* Multi-stepper simulation with fixed point counts */

// [[Rcpp::export]]
NumericMatrix rstepper_multi_fixed_c(List theta1_list,
                                     List theta2_list,
                                     NumericMatrix bbox,
                                     int iter,
                                     NumericMatrix x0,
                                     int dbg,
                                     int toroidal,
                                     List trend_ims) {
  RNGScope scope;
  int i, j, l, k, it, mn, t, t1, t2;
  double xn, yn, xo, yo, pap, change, newpot, pot;
  int acc[2] = {0,0};

  int dbg_step = (int) (1.0*iter/100.0);
  if(dbg_step <1) dbg_step = 1;

  int ntypes = theta1_list.size();
  int npairs = ntypes * (ntypes-1) / 2;

  NumericVector theta0(ntypes);

  // Parse the parameters, change the storage form
  std::vector<std::vector<double> > ranges1(ntypes), ranges2(npairs);
  std::vector<std::vector<double> > thetas1(ntypes), thetas2(npairs);
  std::vector<std::vector<int> > saturations1(ntypes), saturations2(npairs);
  if(dbg>100)Rprintf("Parsing %i intra-parameters\n", ntypes);
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

    // gotta check it worked


    if(dbg>100){
      print_input_parameters(&ranges1,&ranges2, &thetas1,&thetas2, &saturations1, &saturations2);
    }

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
    }
    if(dbg>10) Rprintf("ok\n");
    // point pattern initialised


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

    // trends
    std::vector<Image *> trends(trend_ims.size());
    for(t = 0; t < ntypes; t++) {
      List f = trend_ims(t);
      trends.at(t) = new Image(f);
    }

    if(dbg>10) Rprintf("\nInitial potential");

    /* compute current potential */

    pot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);

    /* here begins the simulation */


    if(dbg > 10) Rprintf("\nStarting fixed n simulation\n");

    // the main loop
    for(it = 0; it < iter; it++ ) {
      Rcpp::checkUserInterrupt();
      // move a random point.
      j = sample_j(x.size());
      // keep the old location
      xo = x.getX(&j);
      yo = x.getY(&j);
      mn = x.getType(&j);
      // sample new location
      xn = runif(1, bbox(0,0),bbox(1,0))(0);
      yn = runif(1, bbox(0,1),bbox(1,1))(0);
      // move
      x.move(&j, xn, yn);
      // update graphs
      for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_move(&j);
      for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_move(&j);
      // new potential
      newpot = potential_multi(&x, theta0, &thetas1, &thetas2, &saturations1, &saturations2, &graphs1, &graphs2);
      change = newpot - pot;
      // trend
      change += trends.at(mn)->getValue(xn,yn) - trends.at(mn)->getValue(xo,yo);
      pap = exp(change);
      pap = fmin(1.0, pap);
      if(runif(1)(0) < pap) { // accept move
        pot = newpot;
        acc[0]++;
      }
      else{ // reject move
        x.move(&j, xo, yo);
        for(t=0; t < ntypes; t++) graphs1.at(t)->update_edges_after_move(&j);
        for(t=0; t < npairs; t++) graphs2.at(t)->update_edges_after_move(&j);
      }
      if(dbg)
        if(it % dbg_step == 0)
        {
          Rprintf("       \r[%i/%i acc %i %4.2f%%]", it, iter, acc[0], 100.0*acc[0]/it);
          if(dbg>10) Rprintf("[pot%10.3f newpot %10.3f]", pot, newpot);
          if(dbg > 1000) usleep(50000);
        }
    }// eof mainloop

    if(dbg) Rprintf("       \r[%i/%i acc %i %4.2f%%]\n", it, iter, acc[0], 100.0*acc[0]/it);

    /* here ends the simulation */

    // done, return coordinates and neighbour counts
    if(dbg>10) Rprintf("building result object");
    NumericMatrix out(x.size(), 3);
    for(i = 0; i < x.size(); i++){
      out(i,0) = x.getX(&i);
      out(i,1) = x.getY(&i);
      out(i,2) = x.getType(&i);
    }
    if(dbg>10) Rprintf(" ok\n");

    return out;
}
