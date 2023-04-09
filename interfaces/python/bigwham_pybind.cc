//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 10.01.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Jan. 12 2021

#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "io/bigwham_io_gen.h"

#include "py_bigwham_helper.h"

namespace py = pybind11;
using namespace bie;
/* -------------------------------------------------------------------------- */

class BigWhamIORect : public BigWhamIOGen {
public:
  BigWhamIORect() {}
  ~BigWhamIORect() {}
};
/* -------------------------------------------------------------------------- */

class PyGetFullBlocks {
private:
  il::Array<double> val_list;
  il::Array<int> rowN;
  il::Array<int> columN;

public:
  PyGetFullBlocks() = default;
  ~PyGetFullBlocks() = default;

  void set(const BigWhamIOGen &BigwhamioObj) {
    il::Array<int> pos_list;
    int nbfentry;

    std::cout << " calling getFullBlocks \n";
    BigwhamioObj.GetFullBlocks(this->val_list, pos_list);
    std::cout << " n entries: " << (this->val_list.size()) << "\n";
    std::cout << " Preparing the vectors \n";

    nbfentry = this->val_list.size();
    this->rowN.Resize(nbfentry);
    this->columN.Resize(nbfentry);

    for (int i = 0; i < nbfentry; i++) {
      this->rowN[i] = pos_list[2 * i];
      this->columN[i] = pos_list[2 * i + 1];
    }
    std::cout << " --- set pyGetFullBlocks completed ---- \n";
  };
  /* --------------------------------------------------------------------------
   */

  pbarray<double> getValList() {
    return as_pyarray<double>(std::move(this->val_list));
  };
  pbarray<int> getColumnN() {
    return as_pyarray<int>(std::move(this->columN));
  };
  pbarray<int> getRowN() { return as_pyarray<int>(std::move(this->rowN)); };
};
/* -------------------------------------------------------------------------- */

PYBIND11_MODULE(py_bigwham, m) {

  //    // Binding the mother class Bigwhamio
  //    // option py::dynamic_attr() added to allow new members to be created
  //    dynamically);
  py::class_<BigWhamIOGen>(m, "BigWhamIOSelf", py::dynamic_attr(),
                           py::module_local())
      .def(py::init<>()) // constructor
      .def("hmat_destructor", &BigWhamIOGen::HmatrixDestructor)
      .def("set", &BigWhamIOGen::SetSelf)
      .def("get_collocation_points", &BigWhamIOGen::GetCollocationPoints)
      .def("get_permutation", &BigWhamIOGen::GetPermutation)
      .def("get_compression_ratio", &BigWhamIOGen::GetCompressionRatio)
      .def("get_kernel_name", &BigWhamIOGen::kernel_name)
      .def("get_spatial_dimension", &BigWhamIOGen::spatial_dimension)
      .def("matrix_size", &BigWhamIOGen::MatrixSize)
      .def("get_hpattern", &BigWhamIOGen::GetHPattern)
      .def("convert_to_global",
           [](BigWhamIOGen &self, const pbarray<double> &x) -> decltype(auto) {
             auto tx = as_array_view<double>(x);
             auto v = self.ConvertToGlobal(tx);
             return as_pyarray<double>(std::move(v));
           })
      .def("convert_to_local",
           [](BigWhamIOGen &self, const pbarray<double> &x) -> decltype(auto) {
             auto tx = as_array_view<double>(x);
             auto v = self.ConvertToLocal(tx);
             return as_pyarray<double>(std::move(v));
           })
      //.def("hdotProductInPermutted", &BigWhamIOGen::hdotProductInPermutted)
      // I change the previous binding of hdotProductInPermutted to return a
      // numpy array!!
      .def(
          "matvec_permute",
          [](BigWhamIOGen &self, const pbarray<double> &x) -> decltype(auto) {
            auto tx = as_array_view<double>(x);
            auto v = self.MatVecPerm(tx);
            return as_pyarray<double>(std::move(v));
          },
          " dot product between hmat and a vector x", py::arg("x"))

      //.def("hdotProduct",            &BigWhamIOGen::matvect, " dot product
      // between hmat and a vector x",py::arg("x"))
      // I change the previous binding of matvect to return a numpy array!!
      // todo: is it possible to move the result of the dot product to an
      // std::array? the array is contiguous in memory but not the vector!!!!!!
      // CP
      .def(
          "matvec",
          [](BigWhamIOGen &self, const pbarray<double> &x) -> decltype(auto) {
            auto tx = as_array_view<double>(x);
            auto v = self.MatVec(tx);
            return as_pyarray<double>(std::move(v));
          },
          " dot product between hmat and a vector x in original ordering",
          py::arg("x"))
      .def("get_hmat_time", &BigWhamIOGen::hmat_time)
      .def("get_omp_threads", &BigWhamIOGen::GetOmpThreads);

  /* --------------------------------------------------------------------------
   */
  py::class_<PyGetFullBlocks>(m, "PyGetFullBlocks")
      .def(py::init<>())
      .def("set", &PyGetFullBlocks::set)
      .def("get_val_list", &PyGetFullBlocks::getValList)
      .def("get_col", &PyGetFullBlocks::getColumnN)
      .def("get_row", &PyGetFullBlocks::getRowN);
  /* --------------------------------------------------------------------------
   */

  declare_array<double>(m, "Real2D");
  declare_array<il::int_t>(m, "Int2D");

  // py::class_<bie::Segment<0>>(m, "Segment0").def(py::init<>());

  /* --------------------------------------------------------------------------
   */

  py::class_<BigWhamIORect>(m, "BigWhamIORect", py::dynamic_attr())
      .def(py::init<>())
      .def("hmat_destructor", &BigWhamIORect::HmatrixDestructor)
      .def("set", &BigWhamIORect::Set)
      .def("get_collocation_points", &BigWhamIORect::GetCollocationPoints)
      .def("get_permutation", &BigWhamIORect::GetPermutation)
      .def("get_compression_ratio", &BigWhamIORect::GetCompressionRatio)
      .def("get_kernel_name", &BigWhamIORect::kernel_name)
      .def("get_spatial_dimension", &BigWhamIORect::spatial_dimension)
      .def("matrix_size", &BigWhamIORect::MatrixSize)
      .def("get_hpattern", &BigWhamIORect::GetHPattern)
      .def("convert_to_global",
           [](BigWhamIORect &self, const pbarray<double> &x) -> decltype(auto) {
             auto tx = as_array_view<double>(x);
             auto v = self.ConvertToGlobal(tx);
             return as_pyarray<double>(std::move(v));
           })
      .def("convert_to_local",
           [](BigWhamIORect &self, const pbarray<double> &x) -> decltype(auto) {
             auto tx = as_array_view<double>(x);
             auto v = self.ConvertToLocal(tx);
             return as_pyarray<double>(std::move(v));
           })
      .def(
          "matvec",
          [](BigWhamIORect &self, const pbarray<double> &x) -> decltype(auto) {
            auto tx = as_array_view<double>(x);
            auto v = self.MatVec(tx);
            return as_pyarray<double>(std::move(v));
          },
          " dot product between hmat and a vector x in original ordering")
      .def("get_hmat_time", &BigWhamIORect::hmat_time)
      .def("get_omp_threads", &BigWhamIORect::GetOmpThreads);
  /* --------------------------------------------------------------------------
   */

  m.def("py_get_collocation_points", &PyGetCollocationPoints);
  m.def("py_get_permutation", &PyGetPermutation);
}
