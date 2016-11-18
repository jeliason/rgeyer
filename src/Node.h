#include <vector>

#ifndef NODE_H_
#define NODE_H_

class Node {
  std::vector<Node *> neighbours;
  public:
    Node();

  virtual ~Node();
  void add_neighbour(Node *);
  void delete_neighbour(Node *);
  void clear_neighbourhood();
  int neighbour_count();
  int neighbour_count_of_neighbhour(int);
  void kill_me();
};

#endif
