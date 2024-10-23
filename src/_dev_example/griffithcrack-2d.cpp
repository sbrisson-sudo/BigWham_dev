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

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "BigWhamIO.h"
#include <il/Array2D.h>

using namespace Eigen;
using Int2D = Array<int, Dynamic, Dynamic, RowMajor>;
using Real2D = Array<double, Dynamic, Dynamic, RowMajor>;
using Real = double;

int main(int argc, char *argv[]) {

  // Using 2DP1 element

  int nelems = 5;
  int dim = 2;

  il::Array2D<double> Coor{nelems + 1, dim, 0.};
  il::Array2D<il::int_t> Conn{nelems, 2, 0};

  Real E = 1.0;
  Real nu = 0.2;
  std::vector<Real> properties = {E, nu};

  std::vector<Real> pts;
  std::vector<int> conn;

  std::string kernel = "2DP1";

  int nleaf = 16;
  Real eta = 30.0;
  Real eps_aca = 0.0001;

  pts.resize((nelems + 1) * dim);
  conn.resize(nelems * dim);

  for (int i = 0; i < nelems + 1; ++i) {
    pts[i * dim + 0] = -1.0 + (2.0 * i) / nelems;
    pts[i * dim + 1] = 0.;
  }

  for (int i = 0; i < nelems; ++i) {
    conn[i * dim + 0] = i;
    conn[i * dim + 1] = i + 1;
  }

  // auto hmat = std::unique_ptr<Bigwhamio>();

  Bigwhamio hmat;

  std::cout << Map<Real2D>(pts.data(), nelems + 1, dim) << std::endl;
  std::cout << Map<Int2D>(conn.data(), nelems, dim) << std::endl;

  hmat.set(pts, conn, kernel, properties, nleaf, eta, eps_aca);

  return 0;
}
