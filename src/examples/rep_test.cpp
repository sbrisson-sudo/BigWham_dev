///=============================================================================
///
/// \file        rep_test.cpp
///
/// \author      Ankit
///
/// \copyright   Copyright (©) 2018 EPFL (Ecole Polytechnique Fédérale
///              de Lausanne)\n
///              Geo-Energy lab
///
/// \brief       Test for Penny shape crack for profiling bigwham
///
///=============================================================================

#include "cnpy.h"
#include "core/BEMesh.h"
#include "core/ElasticProperties.h"
#include "src/core/SquareMatrixGenerator.h"
#include "src/core/elements/Triangle.h"
#include "src/core/hierarchical_representation.h"
#include <src/elasticity/BIE_elastostatic.h>
#include "src/elasticity/3d/BIE_elastostatic_triangle_0_impls.h"
#include "src/hmat/hmatrix/Hmat.h"
#include <algorithm>
#include <cmath>
#include <il/Array2D.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <string>
#include <vector>

using Real1D = il::Array<double>;
using Real2D = il::Array2D<double>;
using Int1D = il::Array<il::int_t>;
using Int2D = il::Array2D<il::int_t>;
using Tri0 = bie::Triangle<0>;
using Mesh = bie::BEMesh<Tri0>;
using MatProp = bie::ElasticProperties;
using KernelH = bie::BIE_elastostatic<Tri0, Tri0, bie::ElasticKernelType::H>;
using MatixGenerator = bie::SquareMatrixGenerator<double, Tri0, KernelH>;

template <typename T> std::string print_array2D(const il::Array2D<T> &);
template <typename T> std::string print_array1D(const il::Array<T> &);
template <typename T>
void copy_array2D(il::Array2D<T> &, const cnpy::NpyArray &);

int main(int argc, char *argv[]) {

  std::string f_coord = "../../../src/examples/mesh_coords.npy";
  std::string f_conn = "../../../src/examples/mesh_conn.npy";

  auto coord_npy = cnpy::npy_load(f_coord);
  auto conn_npy = cnpy::npy_load(f_conn);

  Real2D coord;
  Int2D conn;

  // Penby shape geometery
  uint dim = 3;
  double radius = 1.0;
  double pressure = 1.0;

  // Elastic properties
  double E = 1.0;
  double nu = 0.0;
  double pi = M_PI;

  // H-mat parameters
  il::int_t max_leaf_size = 32;
  double eta = 2.0;
  double eps_aca = 1.e-3;

  copy_array2D(coord, coord_npy);
  copy_array2D(conn, conn_npy);

  std::cout << print_array2D(coord) << std::endl;
  std::cout << print_array2D(conn) << std::endl;

  Mesh my_mesh(coord, conn);
  Real2D xcol = my_mesh.getCollocationPoints();
  MatProp elas(E, nu);
  KernelH ker(elas, coord.size(1));
  Real1D prop{1, 1000.};
  ker.setKernelProperties(prop);

  bie::HRepresentation hr =
      bie::h_representation_square_matrix(my_mesh, max_leaf_size, eta);

  MatixGenerator M(my_mesh, ker, hr.permutation_0_);
  bie::Hmat<double> h_(M, hr, eps_aca);

  Real1D cos{M.size(1), 0.0}, trac_perm{M.size(0), 0.0};
  double pre_fac = (8 * (1 - nu * nu)) / (pi * E);

  double rsq;
  for (il::int_t i = 0; i < M.sizeAsBlocks(1); i++) {
    rsq = xcol(i, 0) * xcol(i, 0) + xcol(i, 1) * xcol(i, 1);
    cos[hr.permutation_1_[i] + 2] = pre_fac * std::sqrt(radius * radius - rsq);
  }
  trac_perm = h_.matvec(cos);
  Real1D trac = trac_perm;

  for (il::int_t i = 0; i < M.sizeAsBlocks(0); i++) {
    for (uint d = 0; d < dim; d++) {
      trac[hr.permutation_0_[i] + d] = trac_perm[i + d];
    }
  }
  std::cout << "Traction \n " << print_array1D(trac) << std::endl;

  Real1D rel_err{M.sizeAsBlocks(0), 0.};
  for (il::int_t i = 0; i < M.sizeAsBlocks(0); i++) {
    rel_err[i] += trac[i + 0] * trac[i + 0];
    rel_err[i] += trac[i + 1] * trac[i + 1] - pressure * pressure;
    rel_err[i] += trac[i + 2] * trac[i + 2];
    rel_err[i] = std::sqrt(rel_err[i]);
    // std::cout << "rel x: " << rel_err[i] << "\n";
  }
  // std::cout << "Linf rel error " << il::norm(rel_err,il::Norm::Linf)
  // <<"\n"; std::cout << "L2 rel error " << il::norm(rel_err,il::Norm::L2)
  // <<"\n";
  std::cout << "Mean rel error " << il::mean(rel_err) << std::endl;

  return 0;
}

template <typename T>
std::string print_array2D(const il::Array2D<T> &instance) {
  std::string t1 = "Size : " + std::to_string(instance.size(0)) + " X " +
                   std::to_string(instance.size(1)) + "  \n\n";
  std::string t2 = "";
  for (il::int_t i = 0; i < instance.size(0); ++i) {
    for (il::int_t j = 0; j < instance.size(1); ++j) {
      t2 += std::to_string(instance(i, j)) + " ";
    }
    t2 += "\n";
  }
  return t1 + t2;
}

template <typename T>
void copy_array2D(il::Array2D<T> &A, const cnpy::NpyArray &B) {
  A.Resize(B.shape[0], B.shape[1]);
  for (uint i = 0; i < B.num_vals; i++) {
    A.Data()[i] = B.data<T>()[i];
  }
}

template <typename T> std::string print_array1D(const il::Array<T> &instance) {
  std::string t1 = "Size : " + std::to_string(instance.size());
  std::string t2 = "";
  for (il::int_t i = 0; i < instance.size(); ++i) {
    t2 += std::to_string(instance[i]) + " ";
    t2 += "\n";
  }
  return t1 + t2;
}
