#include <math.h>
#include <vector>
#include "Node.h"

/********************************************************************************************/
Node::~Node()
{
}
/********************************************************************************************/
Node::Node() {
  neighbours.clear();
}
/********************************************************************************************/
void Node::add_neighbour(Node *n){
  // check new neighbour
  bool old = false;
  for(int i=0; i < neighbours.size(); i++) {
    if(neighbours.at(i) == n) {
      old = true;
      break;
    }
  }
  if(!old) neighbours.push_back(n);
};

/********************************************************************************************/
void Node::clear_neighbourhood(){
  neighbours.clear();
}
void Node::delete_neighbour(Node *n){
  bool found = false;
  int i;
  for(i=0; i < neighbours.size(); i++)
  if( neighbours.at(i) == n ) {
    found = true;
    break;
  }
  if(found) neighbours.erase(neighbours.begin() + i);
};

void Node::kill_me(){
  // gotta tell neighbours to remove me
  int i;
  for(i=0; i < neighbours.size(); i++) neighbours.at(i)->delete_neighbour(this);
  delete this;
};


int Node::neighbour_count(){
  return neighbours.size();
};

int Node::neighbour_count_of_neighbhour(int j){
  return neighbours.at(j)->neighbour_count();
}

