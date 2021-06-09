#include <iostream>

#include "Vertex.h"

// Constructor for Vertex.
// Initialization cannot be specified
Vertex::Vertex()
{
  // label specifying no vertex
  label_ = -1;
  // white
  color_ = 0;
  // -1 means infinity
  distance_ = -1;
  // no predecessor
  parent_ = NULL;
}