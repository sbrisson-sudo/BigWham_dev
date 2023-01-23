///=============================================================================
///
/// \file        griffithcrack-2d.cpp
///
/// \author      Ankit
///
/// \copyright   Copyright (©) 2018 EPFL (Ecole Polytechnique Fédérale
///              de Lausanne)\n
///              Geo-Energy lab
///
/// \brief       Test for 2D Griffith crack
///
///=============================================================================


#include <iostream>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

using Int2D = Array<int, 2, Dynamic>;
using Real2D = Array<double, 2, Dynamic>;
using Real = double;

int main(int argc, char *argv[]) {

  int nelems = 70;
  int dim = 2;

  Real E = 1.0;
  Real nu = 0.2;


  Real2D pts;
  Int2D conn;


  pts.setZero(nelems+1, dim);
  conn.setZero(nelems, dim);

  for (int i = 0; i < nelems+1; ++i) {
    pts(i,0) = -1.0 + (2.0 * i) / nelems;
    pts(i,1) = 0.;
  }

  for (int i = 0; i < nelems; ++i) {
    conn(i,0) = i;
    conn(i,1) = i+1;
  }

  std::cout << "Hello" << std::endl;

  return 0;
}
