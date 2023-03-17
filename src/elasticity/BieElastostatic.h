//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 23.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_BIEELASTOSTATIC_H
#define BIGWHAM_BIEELASTOSTATIC_H

#include <vector> 
#include "core/BieKernel.h"
#include "core/ElasticProperties.h"
#include "core/BoundaryElement.h"

namespace bie{

    enum ElasticKernelType {U,T,H,S,V};


/*

 Class Square Matrix generator for BIE - note that the element type of the source and receiver elements should be the same!

 Es, Er: Element Type, Triangle<0>
 ElasticKernelType: ElasticKernelType {U,T,H,S,V};
*/
    template<class Es,class Er,ElasticKernelType k>
    class BieElastostatic : public BieKernel<double,Es,Er> {
// generic class for uniform  elasticity - isotopric - full-space.
//note that for pure mode 1 (scalar problem), another class must be written or / derived from this one.
    using BieKernel<double,Es,Er>::BieKernel;

    protected:
        il::Array<double> kernel_properties_{};
        bie::ElasticProperties elas_;
        bool local_unknowns_ {true};
        bool local_co_variables_{true};

    public :

        BieElastostatic() : BieKernel<double, Es, Er>()  {};

        BieElastostatic(bie::ElasticProperties &elas,il::int_t dim) : BieKernel<double, Es, Er>() {
            elas_=elas;
            this->dof_dimension_=dim;
            this->dim_=dim;
        };

        BieElastostatic(bie::ElasticProperties &elas ,il::int_t dim,bool local_unknowns,bool local_co_variables)  : BieKernel<double, Es, Er>()  {
            elas_=elas;
            this->dof_dimension_=dim;
            this->dim_=dim;
            local_unknowns_=local_unknowns;
            local_co_variables_=local_co_variables;
        };

        void setKernelProperties( il::Array<double> &prop) {
            kernel_properties_=prop;
        }

        bool isLocalUnknowns() const {return local_unknowns_;};
        bool isLocalCoVariables() const {return local_co_variables_;};

        virtual std::vector<double>  influence(Es source_elt,il::int_t i_s,Er receiver_elt, il::int_t i_r) const;
    };



}

#endif //BIGWHAM_BIEELASTOSTATIC_H
