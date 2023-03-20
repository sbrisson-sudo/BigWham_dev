//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 12.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_BIEKERNEL_H
#define BIGWHAM_BIEKERNEL_H

#include <string>
#include <vector>

#include "core/boundary_element.h"

namespace bie {

template <typename T> class BieKernel {
protected:
  il::int_t dof_dimension_;
  il::int_t dim_;

public:
  //  constructor
  BieKernel(){};
  ~BieKernel(){};

  // methods
  virtual std::vector<T> influence(BoundaryElement &src_element, il::int_t colloc_id_src,
                                   BoundaryElement &rec_element,
                                   il::int_t colloc_id_rec) const = 0;

  il::int_t getDofDimension() const { return dof_dimension_; };
  il::int_t getSpatialDimension() const { return dim_; };
};

} // namespace bie
#endif // BIGWHAM_BIEKERNEL_H
