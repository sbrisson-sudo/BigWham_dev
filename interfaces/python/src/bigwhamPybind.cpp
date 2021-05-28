//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 10.01.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Jan. 12 2021

#include <interfaces/python/pybind11-master/include/pybind11/pybind11.h>
#include <interfaces/python/pybind11-master/include/pybind11/numpy.h>
#include <interfaces/python/pybind11-master/include/pybind11/stl.h>
#include <interfaces/python/pybind11-master/include/pybind11/complex.h>
#include <interfaces/python/pybind11-master/include/pybind11/functional.h>
#include <interfaces/python/pybind11-master/include/pybind11/chrono.h>

#include "BigWham.h"

namespace py = pybind11;

// get fullBlocks
class pyGetFullBlocks {
private:
    std::vector<double> val_list;
    std::vector<int> rowN;
    std::vector<int> columN;

public:
    pyGetFullBlocks() = default;
    ~pyGetFullBlocks() = default;

    void set(Bigwhamio & BigwhamioObj) {
        std::vector<int> pos_list;
        long nbfentry;

        std::cout << " calling getFullBlocks \n";
        BigwhamioObj.getFullBlocks(this->val_list,pos_list);
        std::cout << " n entries: " <<  (this->val_list.size()) << "\n";
        std::cout << " preparing the Mtensors \n";

        nbfentry = this->val_list.size();
        this->rowN.resize(nbfentry);
        this->columN.resize(nbfentry);

        // be careful makeSparseMatrix require 1-indexing for position !
        for (int i=0;i<pos_list.size()/2;i++){
            this->rowN[i]=pos_list[2*i];
            this->columN[i]=pos_list[2*i+1];
        }

        std::cout << " --- set pyGetFullBlocks completed ---- \n";

    };

    std::vector<double> & getgetValList(){return this->val_list;};
    std::vector<int> & getgetColumnN(){return this->columN;};
    std::vector<int> & getgetRowN(){return this->rowN;};

    // the following lines are
    // taken from https://github.com/pybind/pybind11/issues/1042
    py::array getRowN(){
        auto v = new std::vector<int>(getgetRowN());
        this->rowN = std::vector<int>();
        auto capsule = py::capsule(v, [](void *v) { delete reinterpret_cast<std::vector<int>*>(v); });
        return py::array(v->size(), v->data(), capsule);
    };

    py::array getColumnN(){
        auto v = new std::vector<int>(getgetColumnN());
        this->columN = std::vector<int>();
        auto capsule = py::capsule(v, [](void *(v)) { delete reinterpret_cast<std::vector<int>*>(v); });
        return py::array(v->size(), v->data(), capsule);
    };

    py::array getValList(){
        auto v = new std::vector<double>(getgetValList());
        this->val_list = std::vector<double>();
        auto capsule = py::capsule(v, [](void *v) { delete reinterpret_cast<std::vector<double>*>(v); });
        return py::array(v->size(), v->data(), capsule);
    };

};



PYBIND11_MODULE(bigwhamPybind, m) {

//    // Binding the mother class Bigwhamio
//    // option py::dynamic_attr() added to allow new members to be created dynamically);
    py::class_<Bigwhamio>(m, "Bigwhamio", py::dynamic_attr(), py::module_local())
      .def(py::init<>())       // constructor
      .def("set",                    &Bigwhamio::set)
      .def("getCollocationPoints",   &Bigwhamio::getCollocationPoints)
      .def("getPermutation",         &Bigwhamio::getPermutation)
      .def("getCompressionRatio",    &Bigwhamio::getCompressionRatio)
      .def("getKernel",              &Bigwhamio::getKernel)
      .def("getSpatialDimension",    &Bigwhamio::getSpatialDimension)
      .def("matrixSize",             &Bigwhamio::matrixSize)
      .def("getHpattern",            &Bigwhamio::getHpattern)
      .def("hdotProductInPermutted", &Bigwhamio::hdotProductInPermutted)
      .def("hdotProduct",            &Bigwhamio::hdotProduct, " dot product between hmat and a vector x",py::arg("x"))
      .def("computeStresses", &Bigwhamio::computeStresses, "function to compute the stress at a given set of points")
      .def("getInfluenceCoe", &Bigwhamio::getInfluenceCoe)
      .def("computeDisplacements", &Bigwhamio::computeDisplacements);


    py::class_<pyGetFullBlocks>(m, "pyGetFullBlocks")
            .def(py::init<>())
            .def("set",               &pyGetFullBlocks::set)
            .def("getValList",        &pyGetFullBlocks::getValList )
            .def("getColumnN",        &pyGetFullBlocks::getColumnN)
            .def("getRowN",     &pyGetFullBlocks::getRowN);

    }

