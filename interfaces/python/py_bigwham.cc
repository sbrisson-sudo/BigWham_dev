//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 10.01.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications ::July 4, 2024 - Ankit Gupta

#include <memory>
#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "io/bigwham_io.h"
#include "io/bigwham_io_helper.h"

#include "elements/boundary_element.h"
#include "elements/segment.h"
#include "elements/point.h"
#include "elements/triangle.h"

#include "py_bigwham_helper.h"

namespace py = pybind11;
using namespace bigwham;
/* -------------------------------------------------------------------------- */

class BigWhamIORect : public BigWhamIO
{
public:
  BigWhamIORect(const std::vector<double> &coor_src,
                const std::vector<int> &conn_src,
                const std::vector<double> &coor_rec,
                const std::vector<int> &conn_rec, const std::string &kernel,
                const std::vector<double> &properties,const int n_openMP_threads,
                const bool verbose)
      : BigWhamIO(coor_src,
                  conn_src,
                  coor_rec,
                  conn_rec, kernel,
                  properties, n_openMP_threads,
                  verbose) {}

  ~BigWhamIORect() {}
};
/* -------------------------------------------------------------------------- */

class PyGetFullBlocks
{
private:
  il::Array<double> val_list;
  il::Array<int> rowN;
  il::Array<int> columN;

public:
  PyGetFullBlocks() = default;
  ~PyGetFullBlocks() = default;

  void set(const BigWhamIO &BigwhamioObj, bool original_order)
  {
    il::Array<int> pos_list;
    int nbfentry;

    if (original_order){
      BigwhamioObj.GetFullBlocks(this->val_list, pos_list);
    } else {
      BigwhamioObj.GetFullBlocksPerm(this->val_list, pos_list);
    }
    
    nbfentry = this->val_list.size();
    this->rowN.Resize(nbfentry);
    this->columN.Resize(nbfentry);

    for (int i = 0; i < nbfentry; i++)
    {
      this->rowN[i] = pos_list[2 * i];
      this->columN[i] = pos_list[2 * i + 1];
    }
  };

  void setRect(const BigWhamIORect &BigwhamioObj, bool original_order)
  {
    il::Array<int> pos_list;
    int nbfentry;

    if (original_order){
      BigwhamioObj.GetFullBlocks(this->val_list, pos_list);
    } else {
      BigwhamioObj.GetFullBlocksPerm(this->val_list, pos_list);
    }

    nbfentry = this->val_list.size();
    this->rowN.Resize(nbfentry);
    this->columN.Resize(nbfentry);

    for (int i = 0; i < nbfentry; i++)
    {
      this->rowN[i] = pos_list[2 * i];
      this->columN[i] = pos_list[2 * i + 1];
    }
  };

  
  /* --------------------------------------------------------------------------
   */

  pbarray<double> getValList()
  {
    return as_pyarray<double>(std::move(this->val_list));
  };
  pbarray<int> getColumnN()
  {
    return as_pyarray<int>(std::move(this->columN));
  };
  pbarray<int> getRowN() { 
    return as_pyarray<int>(std::move(this->rowN)); 
  };
};
/* -------------------------------------------------------------------------- */

PYBIND11_MODULE(py_bigwham, m)
{
  //    // Binding the mother class Bigwhamio_old
  //    // option py::dynamic_attr() added to allow new members to be created
  //    dynamically);
  // Square Self Interaction matrices
  py::class_<BigWhamIO>(m, "BigWhamIOSelf", py::dynamic_attr(),
                        py::module_local())
      .def(py::init<const std::vector<double> &,
                    const std::vector<int> &, const std::string &,
                    const std::vector<double> &,const int &,const bool>()) // constructor
      .def("hmat_destructor", &BigWhamIO::HmatrixDestructor)
      .def("load_from_file", &BigWhamIO::LoadFromFile)
      .def("build_hierarchical_matrix", &BigWhamIO::BuildHierarchicalMatrix)
      .def("build_pattern", &BigWhamIO::BuildPattern)
      .def("get_collocation_points", &BigWhamIO::GetCollocationPoints)
      .def("get_permutation", &BigWhamIO::GetPermutation)
      .def("get_compression_ratio", &BigWhamIO::GetCompressionRatio)
      .def("get_kernel_name", &BigWhamIO::kernel_name)
      .def("get_spatial_dimension", &BigWhamIO::spatial_dimension)
      .def("matrix_size", &BigWhamIO::MatrixSize)
      .def("get_hpattern", &BigWhamIO::GetHPattern)
      .def("write_hmatrix", &BigWhamIO::WriteHmatrix)
      .def("get_hmat_time", &BigWhamIO::hmat_time)
      .def("get_omp_threads", &BigWhamIO::GetOmpThreads)
      .def("get_element_normals",&BigWhamIO::GetElementNormals)
      .def("get_rotation_matrix",&BigWhamIO::GetRotationMatrix)
      .def(py::pickle(
          [](const BigWhamIO &self) { // __getstate__
            /* Return a tuple that fully encodes the state of the object */
            return py::make_tuple(self.kernel_name());
          },
          [](py::tuple t) { // __setstate__
            if (t.size() != 1)
              throw std::runtime_error("Invalid state!");

            /* Create a new C++ instance */
            BigWhamIO p;
            return p;
          }))
      .def("convert_to_global",
           [](BigWhamIO &self, const pbarray<double> &x) -> decltype(auto)
           {
             auto tx = as_array_view<double>(x);
             auto v = self.ConvertToGlobal(tx);
             return as_pyarray<double>(std::move(v));
           })
      .def("convert_to_local",
           [](BigWhamIO &self, const pbarray<double> &x) -> decltype(auto)
           {
             auto tx = as_array_view<double>(x);
             auto v = self.ConvertToLocal(tx);
             return as_pyarray<double>(std::move(v));
           })
      //return a numpy array!!
      .def(
          "matvec_permute",
          [](BigWhamIO &self, const pbarray<double> &x) -> decltype(auto)
          {
            auto tx = as_array_view<double>(x);
            auto v = self.MatVecPerm(tx);
            return as_pyarray<double>(std::move(v));
          },
          " dot product between hmat and a vector x", py::arg("x"))
      // return a numpy array!!
      .def(
          "matvec",
          [](BigWhamIO &self, const pbarray<double> &x) -> decltype(auto)
          {
            auto tx = as_array_view<double>(x);
            il::Array<double> v = self.MatVec(tx);
            return as_pyarray<double>(std::move(v));
          },
          " dot product between hmat and a vector x in original ordering",
          py::arg("x"))
      .def(
           "matvecVoid",
           [](BigWhamIO &self, const pbarray<double> &x, pbarray<double> &y) -> decltype(auto)
           {
               auto tx = as_array_view<double>(x);
               auto ty = as_array_edit<double>(y);
               self.MatVecVoid(tx,ty);
               return;
               },
               " dot product between hmat and a vector x in original ordering - void function",
               py::arg("x"),py::arg("y"))
      .def(
          "compute_displacements",
          [](BigWhamIO &self, const std::vector<double> &coor, const pbarray<double> &x) -> decltype(auto)
          {
            auto tx = as_array_view<double>(x);
            auto v = self.ComputeDisplacements(coor, tx);
            return as_pyarray<double>(std::move(v));
          },
          " compute displacement at set of points",
          py::arg("x"), py::arg("coor"))
      .def(
          "compute_stresses",
          [](BigWhamIO &self, const std::vector<double> &coor, const pbarray<double> &x) -> decltype(auto)
          {
            auto tx = as_array_view<double>(x);
            auto v = self.ComputeStresses(coor, tx);
            return as_pyarray<double>(std::move(v));
          },
          " compute stresses at set of points",
          py::arg("x"), py::arg("coor"))
      .def(
          "get_diagonal", 
          [](const BigWhamIO &self) {
            std::vector<double> val_list;
            self.GetDiagonal(val_list);
            return py::array(val_list.size(), val_list.data());
        });;

  /* --------------------------------------------------------------------------
   */
  py::class_<PyGetFullBlocks>(m, "PyGetFullBlocks")
      .def(py::init<>())
      .def("set", &PyGetFullBlocks::set)
      .def("setRect", &PyGetFullBlocks::setRect)
      .def("get_val_list", &PyGetFullBlocks::getValList)
      .def("get_col", &PyGetFullBlocks::getColumnN)
      .def("get_row", &PyGetFullBlocks::getRowN);
  /* --------------------------------------------------------------------------
   */

//  declare_array<double>(m, "Real2D");
//  declare_array<il::int_t>(m, "Int2D");

  /* --------------------------------------------------------------------------
   */

  py::class_<BigWhamIORect>(m, "BigWhamIORect", py::dynamic_attr())
      .def(py::init<const std::vector<double> &, const std::vector<int> &,
                    const std::vector<double> &, const std::vector<int> &,
                    const std::string &, const std::vector<double> &, const int &, const bool>())
      .def("hmat_destructor", &BigWhamIORect::HmatrixDestructor)
      .def("build_hierarchical_matrix", &BigWhamIORect::BuildHierarchicalMatrix)
      .def("build_pattern", &BigWhamIORect::BuildPattern)
      .def("load_from_file", &BigWhamIORect::LoadFromFile)
      .def("get_collocation_points", &BigWhamIORect::GetCollocationPoints)
      .def("get_permutation", &BigWhamIORect::GetPermutation)
      .def("get_permutation_receivers", &BigWhamIORect::GetPermutationReceivers)
      .def("get_compression_ratio", &BigWhamIORect::GetCompressionRatio)
      .def("get_kernel_name", &BigWhamIORect::kernel_name)
      .def("get_spatial_dimension", &BigWhamIORect::spatial_dimension)
      .def("matrix_size", &BigWhamIORect::MatrixSize)
      .def("get_hpattern", &BigWhamIORect::GetHPattern)
      .def("write_hmatrix", &BigWhamIO::WriteHmatrix) // check here why BigWhamIO and not BigWhamIORect ?
      .def("get_hmat_time", &BigWhamIORect::hmat_time)
      .def("get_omp_threads", &BigWhamIORect::GetOmpThreads)
      .def("convert_to_global",
           [](BigWhamIORect &self, const pbarray<double> &x) -> decltype(auto)
           {
             auto tx = as_array_view<double>(x);
             auto v = self.ConvertToGlobal(tx);
             return as_pyarray<double>(std::move(v));
           })
      .def("convert_to_local",
           [](BigWhamIORect &self, const pbarray<double> &x) -> decltype(auto)
           {
             auto tx = as_array_view<double>(x);
             auto v = self.ConvertToLocal(tx);
             return as_pyarray<double>(std::move(v));
           })
      .def(
          "matvec",
          [](BigWhamIORect &self, const pbarray<double> &x) -> decltype(auto)
          {
            auto tx = as_array_view<double>(x);
            auto v = self.MatVec(tx);
            return as_pyarray<double>(std::move(v));
          },
          " dot product between hmat and a vector x in original ordering")
      .def(
        "get_diagonal", 
        [](const BigWhamIORect &self) {
          std::vector<double> val_list;
          self.GetDiagonal(val_list);
          return py::array(val_list.size(), val_list.data());
      });


/* --------------------------------------------------------------------------*/
/* Misc functions can be used without constructing BigWhamIO object*/

  m.def("py_get_collocation_points", &PyGetCollocationPoints);
  m.def("py_get_permutation", &PyGetPermutation);

}
