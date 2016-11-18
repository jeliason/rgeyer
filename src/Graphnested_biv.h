#include <vector>
#include "Pplite.h"
#include "Node.h"

#ifndef GNB_H_
#define GNB_H_

class Graphnested_biv {
  std::vector<std::vector<Node *> > nodes;
  std::vector<double> r; // increasing range vector for geometric graph
  Pplite *x;
  bool is_empty;
  public:
    Graphnested_biv(Pplite *x, std::vector<double> r);

  virtual ~Graphnested_biv();
  void compute_edges();
  void update_edges_after_addition(); // after a point has been added
  void update_edges_after_move(int *i);  // after a point has been moved
  void update_edges_after_delete(int *i); // after deletion of a point
  int neighbour_count_of_i_in_graph_k(int k, int i);
  int has_neighbours_i_in_graph_k(int k, int i);
  int neighbour_count_of_neighbour_j_of_i_in_graph_k(int k, int i, int j);
};

#endif
