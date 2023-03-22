//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 23.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_BIEELASTOSTATIC_H
#define BIGWHAM_BIEELASTOSTATIC_H

#include "core/bie_kernel.h"
#include "core/elastic_properties.h"
#include "elements/boundary_element.h"
#include <vector>

namespace bie {

enum class ElasticKernelType { U, T, H, S, V };

/*

 Class Square Matrix generator for BIE - note that the element type of the
 source and receiver elements should be the same!

 Es, Er: Element Type, Triangle<0>
 ElasticKernelType: ElasticKernelType {U,T,H,S,V};
*/
template <class Es, class Er, ElasticKernelType k>
class BieElastostatic : public BieKernel<double> {
  // generic class for uniform  elasticity - isotopric - full-space.
  // note that for pure mode 1 (scalar problem), another class must be written
  // or / derived from this one.

protected:
  il::Array<double> kernel_properties_;
  bie::ElasticProperties elas_;
  bool local_unknowns_;
  bool local_co_variables_;

public:
  BieElastostatic() : BieKernel<double>(){};

  BieElastostatic(const bie::ElasticProperties &elas, const il::int_t dim)
      : BieKernel<double>() {
    elas_ = elas;
    this->dof_dimension_ = dim;
    this->dim_ = dim;
    this->local_unknowns_ = true;
    this->local_co_variables_ = true;
  }
  BieElastostatic(bie::ElasticProperties &elas, il::int_t dim,
                  bool local_unknowns, bool local_co_variables)
      : BieKernel<double>() {
    elas_ = elas;
    this->dof_dimension_ = dim;
    this->dim_ = dim;
    local_unknowns_ = local_unknowns;
    local_co_variables_ = local_co_variables;
  };

  void set_kernel_properties(const il::Array<double> &prop) {
    kernel_properties_ = prop;
  }

  bool is_local_unknowns() const { return local_unknowns_; };
  bool is_local_covariables() const { return local_co_variables_; };

  virtual std::vector<double> influence(const BoundaryElement &source_elt,
                                        il::int_t i_s,
                                        const BoundaryElement &receiver_elt,
                                        il::int_t i_r) const override;
};

} // namespace bie

#endif // BIGWHAM_BIEELASTOSTA
