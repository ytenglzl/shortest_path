#pragma once

// Vertex class
class Vertex
{
public:
  // label
  int label_;
  //0 white, 1 grey, 2 black
  int color_;
  //distance
  int distance_;
  //pointer to predecessor
  Vertex *parent_;

  Vertex();
};