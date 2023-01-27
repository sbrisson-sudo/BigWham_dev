//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 23.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_BIE_ELASTOSTATIC_H
#define BIGWHAM_BIE_ELASTOSTATIC_H

#include <src/core/BIE_Kernel.h>
#include <src/core/ElasticProperties.h>
#include <src/core/BoundaryElement.h>

#include <src/elasticity/2d/ElasticS3DP0_element.h>

namespace bie{

    enum ElasticKernelType {U,T,H,S,V};

    template<class Es,class Er,ElasticKernelType k>
    class BIE_elastostatic : public BIE_Kernel<double,Es,Er> {
// generic class for uniform  elasticity - isotopric - full-space.
//note that for pure mode 1 (scalar problem), another class must be written or / derived from this one.
    using BIE_Kernel<double,Es,Er>::BIE_Kernel;

    protected:
        il::Array<double> kernel_properties_{};
        bie::ElasticProperties elas_;

    public :

        BIE_elastostatic() : BIE_Kernel<double, Es, Er>()  {};

        BIE_elastostatic(bie::ElasticProperties &elas,il::int_t dim) : BIE_Kernel<double, Es, Er>() {
            elas_=elas;
            this->dof_dim_=dim;
            this->dim_=dim;
        };

        BIE_elastostatic(bie::ElasticProperties &elas ,il::int_t dim,bool local_unknowns,bool local_co_variables)  : BIE_Kernel<double, Es, Er>(local_unknowns,local_co_variables)  {
            elas_=elas;
            this->dof_dim_=dim;
            this->dim_=dim;
        };

        void setKernelProperties( il::Array<double> &prop) {
            kernel_properties_=prop;
        }

         std::vector<double>  influence(Es source_elt,il::int_t i_s,Er receiver_elt, il::int_t i_r) const { } ;

    };



}

#endif //BIGWHAM_BIE_ELASTOSTATIC_H
