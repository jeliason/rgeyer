// #include <Rcpp.h>
// #include <unistd.h>
// #include <vector>
// #include "Pplite.h"
// #include "Graphnested.h"
// #include "helpers.h"
//
// using namespace Rcpp;
//
//
// // For non-data points
//
// // [[Rcpp::export]]
// NumericMatrix rstepper_components_multi_grid_c(
//                                               NumericMatrix from,
//                                               NumericMatrix grid,
//                                               List ranges1l,
//                                               List ranges2l,
//                                               List sat1l,
//                                               List sat2l,
//                                               NumericMatrix bbox,
//                                               int dbg,
//                                               int toroidal) {
//
//   RNGScope scope;
//   int i, j, l, k, mn, t, t1, t2;
//   double xn, yn, newpot, pot, change;
//
//   int ntypes = ranges1l.size();
//   int npairs = ntypes * (ntypes-1) / 2;
//
//
//   // init point pattern window here so don't need dependence on Rcpp in Pplite
//   std::vector<double> win(bbox.ncol()*2);
//   double Area=1;
//   for(i=0; i < bbox.ncol(); i++) {
//     win.at(2*i) = bbox(0, i);
//     win.at(2*i+1) = bbox(1, i);
//   }
//
//   // Parse the parameters, change the storage form
//   std::vector<std::vector<double> > ranges1(ntypes), ranges2(npairs);
//   std::vector<std::vector<int> > saturations1(ntypes), saturations2(npairs);
//   // intra-type parameters
//   List lt;
//   SEXP lts;
//   NumericVector y;
//   for(i=0; i < ntypes; i++) {
//     lts = ranges1l(i);
//     y = lts;
//     for(j = 0; j < y.size(); j++)ranges1.at(i).push_back(y(j));
//     lts = sat1l(i);
//     y = lts;
//     for(j = 0; j < y.size(); j++)saturations1.at(i).push_back(y(j));
//   }
//   IntegerMatrix pairs(npairs, 2); // tracking which marks are in which pair.
//   if(dbg>100)Rprintf("Parsing %i inter-parameters\n", npairs);
//   k = 0; // over pairs
//   for(t1=0; t1 < ntypes-1; t1++)
//     for(t2 = t1+1; t2 < ntypes; t2++) {
//       pairs(k,0) = t1;
//       pairs(k,1) = t2;
//       lts = ranges2l(i); // casting
//       y = lts;
//       for(l = 0; l < y.size(); l++) ranges2.at(k).push_back(y(l));
//       lts = sat2l(i);
//       y = lts;
//       for(l = 0; l < y.size(); l++)saturations2.at(k).push_back(y(l));
//       k++;
//     }
//
//
//   // init point pattern object
//   if(dbg>10) Rprintf("Initialising Pplite");
//   Pplite x(win, toroidal);
//   for(i = 0; i < from.nrow(); i++) x.push_back(from(i,0), from(i,1), (int) from(i,2));
//   if(dbg>10) Rprintf("ok\n");
//   // point pattern initialised
//
//   // initialise graphs
//
//   // initialise graphs. ntypes + npairs of them.
//   std::vector<Graphnested *> graphs1(ntypes);
//   std::vector<Graphnested *> graphs2(npairs);
//
//
//   // first the intra-type graphs
//   if(dbg>10) Rprintf("Initialising Graphnested, intra (%i)", graphs1.size());
//   std::vector<int> target(1);
//   for(t = 0; t < ntypes; t++) {
//     graphs1.at(t) = new Graphnested(&x, ranges1.at(t));
//     target.at(0) = t; // note that types are 0,..., ntypes-1 integers
//     graphs1.at(t)->set_pars_i(target);
//     graphs1.at(t)->set_pair_suitable_method(&Graphnested::pair_suitable_same_mark_but_of_given_set);
//     graphs1.at(t)->compute_edges();
//   }
//
//
//
//
//
//   NumericMatrix components(to.nrow(), K);
//
//   NumericMatrix original_factors(from.nrow(), K);
//
//   for(j = 0; j < x.size(); j++)
//     for(k = 0; k < K; k++)
//       original_factors(j, k) = imin(sat(k), graphs.neighbour_count_of_i_in_graph_k(k, j));
//
//   /* here begins the evaluation */
//
//   // the main loop
//   for(i = 0; i < to.nrow(); i++ ) {
//     // new location
//     xn = to(i,0);
//     yn = to(i,1);
//     mn = to(i,2);
//     // add new point and compute the unweighted potential difference
//     x.push_back(xn, yn, mn);
//     l = x.size()-1;
//     graphs.update_edges_after_addition();
//
//     for(k = 0; k < K; k++) {
//       //change = graphs.has_neighbours_i_in_graph_k(k, l); // the new point
//       change = imin(sat(k), graphs.neighbour_count_of_i_in_graph_k(k, l));
//       for(j = 0; j < l; j++) { // and others
//         //change += graphs.has_neighbours_i_in_graph_k(k, j) - original_factors(j, k);
//         change += imin(sat(k), graphs.neighbour_count_of_i_in_graph_k(k, j)) - original_factors(j, k);
//       }
//       components(i,k) = change;
//     }
//
//     // restore
//     x.pop_back();
//     graphs.update_edges_after_delete(&l);
//     if(dbg) {
//       Rprintf("       \r[%i/%i]", i+1, to.nrow(), x.size());
//     }
//   }// eof mainloop
//
//
//   if(dbg) Rprintf("\n");
//   return components;
// }
//
