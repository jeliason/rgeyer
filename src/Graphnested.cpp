#include <math.h>
#include <vector>
#include "Graphnested.h"

/********************************************************************************************/
Graphnested::~Graphnested()
{
}
/********************************************************************************************/
Graphnested::Graphnested(Pplite *x, std::vector<double> r) {
  this->x = x;
  this->r = r;
  nodes.resize(r.size());
  is_empty = (r.at(0) == 0) & (r.size() == 1);
  int i,j;
  for(i=0; i < nodes.size(); i++) {
    nodes.at(i) = std::vector<Node *> ( x->size() );
    for(j=0; j < x->size(); j++) nodes.at(i).at(j) = new Node();
  }
  // default difference:
  set_pair_suitable_method(&Graphnested::pair_suitable_no_marks);
}

/********************************************************************************************/

void Graphnested::set_pars_i(std::vector<int > n){
  pars_i = n;
};
void Graphnested::set_pars_d(std::vector<double > n){
  pars_d = n;
}

bool Graphnested::suitable_mark(int i){
  return end(pars_i) != std::find(begin(pars_i), end(pars_i), x->getType(&i));
}
/********************************************************************************************/
void Graphnested::compute_edges() {
  if(is_empty) return;
  int i,j,k;
  double d;
  // pairwise, no other way
  for(i = 0; i < x->size()-1; i++) {
    for(j = i+1; j < x->size(); j++)
      if(pair_suitable(i,j))
      {
        d = x->getDist(&i, &j);
        // in any interval
        if(d < r.at(r.size()-1)) {
          //go forwards in ranges
          for(k=0; k < r.size(); k++) {
            if(d < r.at(k)) {
              break; // found the interval.
            }
          }
          nodes.at(k).at(i)->add_neighbour( nodes.at(k).at(j) );
          nodes.at(k).at(j)->add_neighbour( nodes.at(k).at(i) );
        }
      }
  }
}

/********************************************************************************************/
void Graphnested::update_edges_after_delete(int *i){
  if(is_empty) return;
  // point i was removed. It is still in matched location
  for(int k = 0; k < r.size(); k++) {
    nodes.at(k).at(*i)->kill_me();
    nodes.at(k).erase(nodes.at(k).begin() + *i);
  }
}
/********************************************************************************************/
void Graphnested::update_edges_after_addition(){
  if(is_empty) return;
  // assume the last point of the pattern is the point
  int j,k, i = x->size()-1;
  double d;

  // add to nodelists
  for(k = 0; k < r.size(); k++) nodes.at(k).push_back(new Node());

  // check the relations
  for(j = 0; j < i; j++) if(pair_suitable(i,j)) { // go through the old points
    d = x->getDist(&i, &j);
    if(d < r.at(r.size()-1)) {// distance in one of the intervals
      //go forwards in ranges
      for(k=0; k < r.size(); k++) {
        if(d < r.at(k)) {
          break; // found the interval.
        }
      }
      nodes.at(k).at(i)->add_neighbour( nodes.at(k).at(j) );
      nodes.at(k).at(j)->add_neighbour( nodes.at(k).at(i) );
    }
  }
}

/********************************************************************************************/
void Graphnested::update_edges_after_move(int *i){
  if(is_empty) return;
  // this is trickier
  int j,k;
  double d;
  // remove old relations of node i
  for(k=0; k < r.size();k++) nodes.at(k).at(*i)->clear_neighbourhood();

  for(j = 0; j < x->size(); j++) { // go through the old points
    if(pair_suitable(*i,j)){ // others

      // remove possible relations with i and j
      for(k=0; k < r.size();k++) nodes.at(k).at(j)->delete_neighbour(nodes.at(k).at(*i));

      // then add new relations i~j if applicable
      d = x->getDist(i, &j);
      if(d < r.at(r.size()-1)) {// distance in one of the intervals
        //go forwards in ranges
        for(k=0; k < r.size(); k++) {
          if(d < r.at(k)) {
            break; // found the interval.
          }
        }
        nodes.at(k).at(*i)->add_neighbour( nodes.at(k).at(j) );
        nodes.at(k).at(j)->add_neighbour( nodes.at(k).at(*i) );
      }
    }
  }
}
/***********************************************************/
int Graphnested::neighbour_count_of_i_in_graph_k(int k, int i){
  if(is_empty) return 0;
  return nodes.at(k).at(i)->neighbour_count();
};

int Graphnested::has_neighbours_i_in_graph_k(int k, int i){
  if(is_empty) return 0;
  if(neighbour_count_of_i_in_graph_k(k,i) > 0) return 1;
  return 0;
}

int Graphnested::neighbour_count_of_neighbour_j_of_i_in_graph_k(int k, int i, int j){
  if(is_empty) return 0;
  return nodes.at(k).at(i)-> neighbour_count_of_neighbhour(j);
};
/********************************************************************************************/
void Graphnested::set_pair_suitable_method(bool (Graphnested::*p)(int, int) ){
  pair_suitable_pt  = p;
}


bool Graphnested::pair_suitable(int i, int j){
  return (*this.*pair_suitable_pt)(i, j);
}

bool Graphnested::pair_suitable_no_marks(int i, int j){
  return i != j; // basic
}

bool Graphnested::pair_suitable_same_mark(int i, int j){
  if (i != j) return (x->getType(&i) == x->getType(&j));
  return false;
};

bool Graphnested::pair_suitable_different_mark(int i, int j){
  if (i != j) return (x->getType(&i) != x->getType(&j));
  return false;
}


bool Graphnested::pair_suitable_same_mark_but_of_given_set(int i, int j){
  if ( pair_suitable_same_mark(i, j) ){
    if(suitable_mark(i)) return suitable_mark(j);
  }
  return false;
}



bool Graphnested::pair_suitable_different_mark_but_of_given_set(int i, int j){
  if ( pair_suitable_different_mark(i, j) ){
    if(suitable_mark(i)) return suitable_mark(j);
  }
  return false;
}

