///=============================================================================
///
/// \file        rep_test_rowhmat.cpp
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
#include "hmat/bie_matrix_generator.h"
#include "hmat/hierarchical_representation.h"
#include "hmat/hmatrix/hmat.h"
#include "hmat/square_matrix_generator.h"

#include "core/be_mesh.h"
#include "core/bie_kernel.h"
#include "core/elastic_properties.h"

#include "elements/triangle.h"

#include "cnpy.h"

using namespace bigwham;
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
  il::int_t max_leaf_size = 20;
  double eta = 3.0;
  double eps_aca = 1.e-3;

  copy_array2D(coord, coord_npy);
  copy_array2D(conn, conn_npy);

  // std::cout << print_array2D(coord) << std::endl;
  // std::cout << print_array2D(conn) << std::endl;

  // BEMesh<Triangle<0>> my_mesh(coord, conn);
  auto my_mesh = std::make_shared<BEMesh<Triangle<0>>>(coord, conn);
  my_mesh->ConstructMesh();

  ElasticProperties elas(E, nu);
  auto ker = std::make_shared<
      BieElastostatic<Triangle<0>, Triangle<0>, ElasticKernelType::H>>(elas,
                                                                       dim);

  il::Array2D<double> xcol = my_mesh->collocation_points();

  std::cout << "Number of Collocation points  = " << xcol.size(0) << " X "
            << xcol.size(1) << std::endl;

  int nelem = xcol.size(0);
  int nvert = xcol.size(1);
  int nelem0 = nelem / 3;
  int nelem1 = nelem - nelem0;

  std::vector<double> dd(nelem * dim, 0.0);
  std::vector<double> trac(nelem * dim, 0.0);

  // DD boundary conditions
  double pre_fac = (8 * (1 - nu * nu)) / (pi * E);
  double rsq;
  for (il::int_t i = 0; i < nelem; i++) {
    rsq = xcol(i, 0) * xcol(i, 0) + xcol(i, 1) * xcol(i, 1);
    dd[dim * i + 2] = pre_fac * std::sqrt(radius * radius - rsq);
  }

  il::Array2D<il::int_t> conn0{nelem0, nvert, 0};
  il::Array2D<il::int_t> conn1{nelem1, nvert, 0};

  for (int j = 0; j < nvert; ++j) {
    for (int i = 0; i < nelem0; ++i) {
      conn0(i, j) = conn(i, j);
    }
  }
  auto mesh_rec0 = std::make_shared<BEMesh<Triangle<0>>>(coord, conn0);
  mesh_rec0->ConstructMesh();
  auto hr0 =
      HRepresentationRectangularMatrix(my_mesh, mesh_rec0, max_leaf_size, eta);
  BieMatrixGenerator<double> M0(my_mesh, mesh_rec0, ker, hr0);
  bigwham::Hmat<double> hmat0(M0, eps_aca);
  auto trac0 = hmat0.matvecOriginal(dd);

  for (int j = 0; j < nvert; ++j) {
    for (int i = 0; i < nelem1; ++i) {
      conn1(i, j) = conn(i + nelem0, j);
    }
  }
  auto mesh_rec1 = std::make_shared<BEMesh<Triangle<0>>>(coord, conn1);
  mesh_rec1->ConstructMesh();
  auto hr1 =
      HRepresentationRectangularMatrix(my_mesh, mesh_rec1, max_leaf_size, eta);
  BieMatrixGenerator<double> M1(my_mesh, mesh_rec1, ker, hr1);
  bigwham::Hmat<double> hmat1(M1, eps_aca);
  auto trac1 = hmat1.matvecOriginal(dd);

  // std::cout << "Traction \n " << print_array1D(trac) << std::endl;
  cnpy::npy_save("trac.npy", trac.data(),
                 {static_cast<unsigned long>(trac.size())});

  il::Array<double> rel_err{nelem, 0.};
  for (il::int_t i = 0; i < nelem0; i++) {
    rel_err[i] += trac0[dim * i + 0] * trac0[dim * i + 0];
    rel_err[i] += trac0[dim * i + 1] * trac0[dim * i + 1];
    rel_err[i] +=
        (trac0[dim * i + 2] - pressure) * (trac0[dim * i + 2] - pressure);
    rel_err[i] = std::sqrt(rel_err[i]);
  }
  for (il::int_t i = 0; i < nelem1; i++) {
    rel_err[i + nelem0] += trac1[dim * i + 0] * trac1[dim * i + 0];
    rel_err[i + nelem0] += trac1[dim * i + 1] * trac1[dim * i + 1];
    rel_err[i + nelem0] +=
        (trac1[dim * i + 2] - pressure) * (trac1[dim * i + 2] - pressure);
    rel_err[i + nelem0] = std::sqrt(rel_err[i + nelem0]);
  }

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
