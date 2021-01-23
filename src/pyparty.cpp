//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 10.01.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Jan. 12 2021

#include <pybind11-master/include/pybind11/pybind11.h>
#include <pybind11-master/include/pybind11/numpy.h>
#include <pybind11-master/include/pybind11/stl.h>
#include <pybind11-master/include/pybind11/complex.h>
#include <pybind11-master/include/pybind11/functional.h>
#include <pybind11-master/include/pybind11/chrono.h>

#include "BigWham.h"

namespace py = pybind11;

PYBIND11_MODULE(pyparty, m) {

//    // Binding the mother class Bigwhamio
//    // option py::dynamic_attr() added to allow new members to be created dynamically);
    py::class_<Bigwhamio>(m, "Bigwhamio", py::dynamic_attr())
      .def(py::init<>())       // constructor
      .def("set",                    &Bigwhamio::set)
      .def("getCollocationPoints",   &Bigwhamio::getCollocationPoints)
      .def("getPermutation",         &Bigwhamio::getPermutation)
      .def("getCompressionRatio",    &Bigwhamio::getCompressionRatio)
      .def("getKernel",              &Bigwhamio::getKernel)
      .def("getSpatialDimension",    &Bigwhamio::getSpatialDimension)
      .def("matrixSize",             &Bigwhamio::matrixSize)
      .def("getHpattern",            &Bigwhamio::getHpattern)
      .def("getFullBlocks",          &Bigwhamio::getFullBlocks)
      .def("hdotProductInPermutted", &Bigwhamio::hdotProductInPermutted)
      .def("hdotProduct",            &Bigwhamio::hdotProduct, " dot product between hmat and a vector x",py::arg("x"))
      .def("computeStresses", &Bigwhamio::computeStresses, "function to compute the stress at a given set of points")
      .def("getInfluenceCoe", &Bigwhamio::getInfluenceCoe);
    }