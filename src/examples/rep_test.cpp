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
#include <il/Array2D.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using Real2D = il::Array2D<double>;
using Int2D = il::Array2D<il::int_t>;

template <typename T> std::string print_array(const il::Array2D<T> &);

int main(int argc, char *argv[]) {

  std::string f_coord = "mesh_coords.npy";
  std::string f_conn = "mesh_conn.npy";

  // il::Status load_status{};
  // auto load_helper_real = il::LoadHelper<Real2D>();
  // Real2D coord = load_helper_real.load(f_coord, il::io, load_status);

  // auto load_helper_int = il::LoadHelper<Int2D>();
  // Int2D conn = load_helper_int.load(f_conn, il::io, load_status);

  auto coord_npy = cnpy::npy_load(f_coord);
  auto conn_npy = cnpy::npy_load(f_conn);

  Real2D coord{coord_npy.shape[0], coord_npy.shape[1], 0.0};
  for (uint i = 0; i < coord_npy.num_vals; i++) {
    coord.Data()[i] = coord_npy.data<double>()[i];
  }

  Int2D conn{conn_npy.shape[0], conn_npy.shape[1], 0};
  for (uint i = 0; i < conn_npy.num_vals; i++) {
    conn.Data()[i] = conn_npy.data<il::int_t>()[i];
  }

  // std::cout << coord.shape[0] << coord.shape[1] << std::endl;

  std::cout << print_array(coord) << std::endl;
  std::cout << print_array(conn) << std::endl;

  return 0;
}

template <typename T> std::string print_array(const il::Array2D<T> &instance) {
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
