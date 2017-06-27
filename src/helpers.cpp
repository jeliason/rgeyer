#include "helpers.h"

int sample_j(int n){
  return (int) ( Rcpp::runif(1)(0) * n);
}
int imin(int a, int b){
  if(a < b) return a;
  return b;
}

int roundup(double f){
  int k = (int) ceil(f);
  if(k == 0) return 1;
  return k;
}


double potential_multi(Pplite *x, Rcpp::NumericVector theta0,
                       std::vector<std::vector<double > > * thetas1,
                       std::vector<std::vector<double > > * thetas2,
                       std::vector<std::vector<int > > * saturations1,
                       std::vector<std::vector<int > > * saturations2,
                       std::vector<Graphnested *> *graphs1,
                       std::vector<Graphnested *> *graphs2
){
  int i, k, t;
  double pot = 0;
  for(i=0; i < x->size(); i++){
    pot += theta0[x->getType(&i)];
    // intra
    for(t=0; t < thetas1->size(); t++) {
      // over ranges
      for(k = 0; k < thetas1->at(t).size(); k++)
        //pot += thetas1->at(t).at(k) * graphs1->at(t)->has_neighbours_i_in_graph_k(k, i);
        pot += thetas1->at(t).at(k) * imin(graphs1->at(t)->neighbour_count_of_i_in_graph_k(k,i), saturations1->at(t).at(k));
    }
    // inter
    for(t=0; t < thetas2->size(); t++) {
      // over ranges
      for(k = 0; k < thetas2->at(t).size(); k++)
        //pot += thetas2->at(t).at(k) * graphs2->at(t)->has_neighbours_i_in_graph_k(k, i);
        pot += thetas2->at(t).at(k) * imin(graphs2->at(t)->neighbour_count_of_i_in_graph_k(k,i), saturations2->at(t).at(k));
    }
  }
  return pot;
}


double maxv(std::vector<double> v){
  double R=0;
  for(int i=0;i<v.size();i++)if(v.at(i)>R)R=v.at(i);
  return R;
}
