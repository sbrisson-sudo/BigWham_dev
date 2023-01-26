//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 12.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#ifndef BIGWHAM_BIE_KERNEL_H
#define BIGWHAM_BIE_KERNEL_H
#pragma once

#include <string>

#include <src/core/BoundaryElement.h>

namespace bie{


template<typename T,class Es,class Er>
class BIE_Kernel {
protected:
    bool local_unknowns_ {true};
    bool local_co_variables_{true};
    il::int_t dof_dim_{1};
    il::int_t dim_{2};

public:

    //  constructor
    BIE_Kernel(bool local_unknowns,bool local_co_variables)  {
        local_unknowns_=local_unknowns;local_co_variables_=local_co_variables;
    };

    BIE_Kernel()  {};

    // methods
    virtual std::vector<T>  influence(Es source_elt,il::int_t i_s,Er receiver_elt, il::int_t i_r) const =0 ;

    bool isLocalUnknowns() const {return local_unknowns_;};
    bool isLocalCoVariables() const {return local_co_variables_;};
    il::int_t getDofDim() const  {return dof_dim_;};
    il::int_t getSpatialDimension() const {return dim_;};

};
//
//    template<typename T, class Es, class Er>
//    BIE_Kernel<T, Es, Er>::BIE_Kernel() = default;
//
//    template<typename T, class Es, class Er>
//    BIE_Kernel<T, Es, Er>::~BIE_Kernel() =default;


}
#endif //BIGWHAM_BIE_KERNEL_H
