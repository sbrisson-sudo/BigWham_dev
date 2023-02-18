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
        bool local_unknowns_ {true};
        bool local_co_variables_{true};

    public :

        BIE_elastostatic() : BIE_Kernel<double, Es, Er>()  {};

        BIE_elastostatic(bie::ElasticProperties &elas,il::int_t dim) : BIE_Kernel<double, Es, Er>() {
            elas_=elas;
            this->dof_dimension_=dim;
            this->dim_=dim;
        };

        BIE_elastostatic(bie::ElasticProperties &elas ,il::int_t dim,bool local_unknowns,bool local_co_variables)  : BIE_Kernel<double, Es, Er>()  {
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

        virtual std::vector<double>  influence(Es source_elt,il::int_t i_s,Er receiver_elt, il::int_t i_r) const {} ;

    };

    // a dummy derived class for simplified 3D P0 kernel....
    template<class Es,class Er,ElasticKernelType k>
    class BIE_elastostatic_sp3d : public BIE_elastostatic<Es,Er,k> {
        using BIE_elastostatic<Es,Er,k>::BIE_elastostatic;

        public:
            BIE_elastostatic_sp3d() : BIE_elastostatic< Es, Er,k>()  {};

            BIE_elastostatic_sp3d(bie::ElasticProperties &elas,il::int_t dim) : BIE_elastostatic<Es, Er,k>() {
                IL_EXPECT_FAST(dim==2);
                this->elas_=elas;
                this->dof_dimension_=dim;
                this->dim_=dim;
        };

        BIE_elastostatic_sp3d(bie::ElasticProperties &elas ,il::int_t dim,bool local_unknowns,bool local_co_variables)  : BIE_elastostatic<Es, Er,k>()  {
            IL_EXPECT_FAST(dim==2);
            this->elas_=elas;
            this->dof_dimension_=dim;
            this->dim_=dim;
            this->local_unknowns_=local_unknowns;
            this->local_co_variables_=local_co_variables;
        };

        virtual std::vector<double>  influence(Es source_elt,il::int_t i_s,Er receiver_elt, il::int_t i_r) const {} ;

    };



}

#endif //BIGWHAM_BIE_ELASTOSTATIC_H
