//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 23.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE
// file for more details.
//

#pragma once
#ifndef BIGWHAM_BIEELASTOSTATIC_H
#define BIGWHAM_BIEELASTOSTATIC_H

#include <vector>

#include "core/bie_kernel.h"
#include "core/elastic_properties.h"
#include "elements/boundary_element.h"
#include "elements/triangle.h"
#include "elements/point.h"

namespace bigwham {

enum class ElasticKernelType { U, T, H, S, V, W }; // should have T and Ttilde as Tt

/*

 Class  for quasi-static elastic BIE

 Es, Er: Element Type, e.g. Triangle<0>
 ElasticKernelType: ElasticKernelType {U,T,H,S,V};

 (for H-matrix - the 2 element types must be the same usually)
 */
template <class Es, class Er, ElasticKernelType k>
class BieElastostatic : public BieKernel<double> {
  // generic class for uniform  elasticity - isotopric - full-space.
  // note that for pure mode 1 (scalar problem), another class must be written
  // or / derived from this one.

protected:
  bigwham::ElasticProperties elas_;
  bool local_unknowns_;
  bool local_co_variables_;

public:
  BieElastostatic() : BieKernel<double>(){};

  BieElastostatic(const bigwham::ElasticProperties &elas, const il::int_t dim) : BieKernel<double>() {
    elas_ = elas;
    this->dof_dimension_ = dim;
    this->spatial_dimension_ = dim;
    this->local_unknowns_ = true;
    this->local_co_variables_ = true;
  }
  BieElastostatic(bigwham::ElasticProperties &elas, il::int_t dim,bool local_unknowns, bool local_co_variables) : BieKernel<double>() {
    elas_ = elas;
    this->dof_dimension_ = dim;
    this->spatial_dimension_ = dim;
    local_unknowns_ = local_unknowns;
    local_co_variables_ = local_co_variables;
  };

  ~BieElastostatic() {};

  bool local_unknowns() const { return local_unknowns_; };
  bool local_covariables() const { return local_co_variables_; };

  virtual std::vector<double> influence(const BoundaryElement &source_elt,
                                        il::int_t i_s,
                                        const BoundaryElement &receiver_elt,
                                        il::int_t i_r) const override;


};


} // namespace bigwham

#endif
