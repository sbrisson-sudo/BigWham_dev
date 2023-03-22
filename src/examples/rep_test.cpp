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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <math.h>
#include <memory>
#include <string>
#include <vector>

#include <il/Array2D.h>

#include "elasticity/bie_elastostatic.h"
#include "hmat/hierarchical_representation.h"
#include "hmat/hmatrix/Hmat.h"
#include "hmat/square_matrix_generator.h"

#include "core/be_mesh.h"
#include "core/bie_kernel.h"
#include "core/elastic_properties.h"

#include "elements/triangle.h"

#include "cnpy.h"

using namespace bie;
template <typename T> std::string print_array2D(const il::Array2D<T> &);
template <typename T> std::string print_array1D(const il::Array<T> &);
template <typename T>
void copy_array2D(il::Array2D<T> &, const cnpy::NpyArray &);

int main(int argc, char *argv[]) {

  std::string f_coord = "mesh_coords.npy";
  std::string f_conn = "mesh_conn.npy";

  auto coord_npy = cnpy::npy_load(f_coord);
  auto conn_npy = cnpy::npy_load(f_conn);

  il::Array2D<double> coord;
  il::Array2D<il::int_t> conn;

  // Penby shape geometery
  uint dim = 3;
  double radius = 1.0;
  double pressure = 1.0;

  // Elastic properties
  double G = 1.0;
  double nu = 0.25;
  double E = (2 * G) * (1 + nu);
  double pi = M_PI;

  // H-mat parameters
  il::int_t max_leaf_size = 4;
  double eta = 3.0;
  double eps_aca = 1.e-1;

  copy_array2D(coord, coord_npy);
  copy_array2D(conn, conn_npy);

  // std::cout << print_array2D(coord) << std::endl;
  // std::cout << print_array2D(conn) << std::endl;

  // BEMesh<Triangle<0>> my_mesh(coord, conn);
  auto my_mesh = std::make_shared<BEMesh<Triangle<0>>>(coord, conn);
  my_mesh->construct_mesh();

  il::Array2D<double> xcol = my_mesh->collocation_points();

  ElasticProperties elas(E, nu);

  // BieElastostatic<Triangle<0>, Triangle<0>, ElasticKernelType::H> ker(
  //     elas, coord.size(1));

  auto ker = std::make_shared<
      BieElastostatic<Triangle<0>, Triangle<0>, ElasticKernelType::H>>(
      elas, coord.size(1));

  // TODO: Should be default
  il::Array<double> prop{1, 1000.};
  ker->set_kernel_properties(prop);

  std::cout << "Number of Collocation points  = " << xcol.size(0) << " X "
            << xcol.size(1) << std::endl;

  auto hr = h_representation_square_matrix(my_mesh, max_leaf_size, eta);

  SquareMatrixGenerator<double> M(my_mesh, ker, hr);
  bie::Hmat<double> h_(M, hr, eps_aca);

  il::Array<double> dd{M.size(1), 0.0};
  il::Array<double> dd_perm{M.size(1), 0.0};
  il::Array<double> trac{M.size(0), 0.0};
  il::Array<double> trac_perm{M.size(0), 0.0};

  double pre_fac = (8 * (1 - nu * nu)) / (pi * E);

  // std::cout << print_array1D(hr.permutation_1_) << std::endl;
  double rsq;
  for (il::int_t i = 0; i < M.sizeAsBlocks(1); i++) {
    rsq = xcol(i, 0) * xcol(i, 0) + xcol(i, 1) * xcol(i, 1);
    dd[dim * i + 2] = pre_fac * std::sqrt(radius * radius - rsq);
  }

  for (il::int_t i = 0; i < M.sizeAsBlocks(1); i++) {
    il::int_t j = hr.permutation_1_[i];
    for (uint d = 0; d < dim; d++) {
      dd_perm[dim * i + d] = dd[dim * j + d];
    }
  }

  cnpy::npy_save("perm1.npy", hr.permutation_1_.Data(),
                 {static_cast<unsigned long>(hr.permutation_1_.size())});
  cnpy::npy_save("perm0.npy", hr.permutation_0_.Data(),
                 {static_cast<unsigned long>(hr.permutation_0_.size())});
  cnpy::npy_save("collocation.npy", xcol.Data(),
                 {static_cast<unsigned long>(xcol.size(0) * xcol.size(1))});
  cnpy::npy_save("dd.npy", dd.Data(), {static_cast<unsigned long>(dd.size())});
  cnpy::npy_save("dd_perm.npy", dd_perm.Data(),
                 {static_cast<unsigned long>(dd_perm.size())});

  // std::cout << "COS \n " << cos.size() << std::endl;
  // std::cout << print_array1D(cos) << std::endl;

  trac_perm = h_.matvec(dd_perm);

  // std::cout << "Traction Perm \n " << print_array1D(trac_perm) << std::endl;
  for (il::int_t i = 0; i < M.sizeAsBlocks(0); i++) {
    il::int_t j = hr.permutation_0_[i];
    for (uint d = 0; d < dim; d++) {
      trac[dim * j + d] = trac_perm[dim * i + d];
    }
  }
  // std::cout << "Traction \n " << print_array1D(trac) << std::endl;
  cnpy::npy_save("trac.npy", trac.Data(),
                 {static_cast<unsigned long>(trac.size())});
  cnpy::npy_save("trac_perm.npy", trac_perm.Data(),
                 {static_cast<unsigned long>(trac_perm.size())});

  il::Array<double> rel_err{M.sizeAsBlocks(0), 0.};
  for (il::int_t i = 0; i < M.sizeAsBlocks(0); i++) {
    rel_err[i] += trac[dim * i + 0] * trac[dim * i + 0];
    rel_err[i] += trac[dim * i + 1] * trac[dim * i + 1];
    rel_err[i] +=
        (trac[dim * i + 2] - pressure) * (trac[dim * i + 2] - pressure);
    rel_err[i] = std::sqrt(rel_err[i]);
    // std::cout << "rel x: " << rel_err[i] << "\n";
  }
  // std::cout << "Linf rel error " << il::norm(rel_err,il::Norm::Linf)
  // <<"\n"; std::cout << "L2 rel error " << il::norm(rel_err,il::Norm::L2)
  // <<"\n";
  std::cout << "Mean rel error " << il::mean(rel_err) << std::endl;
  std::cout << "L2 rel error " << il::norm(rel_err, il::Norm::L2) << std::endl;

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
  std::string t1 = "Size : " + std::to_string(instance.size()) + "\n\n";
  std::string t2 = "";
  for (il::int_t i = 0; i < instance.size(); ++i) {
    t2 += std::to_string(instance[i]) + " ";
    t2 += "\n";
  }
  return t1 + t2;
}
