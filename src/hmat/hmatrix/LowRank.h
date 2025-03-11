//
// This file is part of BigWham.
//
// Created by Francois Fayard - 2018
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#pragma once

#ifndef BIGWHAM_LOWRANK_H
#define BIGWHAM_LOWRANK_H

#include <il/Array2D.h>

namespace bigwham {

template <typename T>
struct LowRank {
  il::Array2D<T> A;
  il::Array2D<T> B;

  // // Constructor to initialize A and B with specific alignment
  // LowRank() 
  //   : A(0, 0, il::align, 64),
  //     B(0, 0, il::align, 64) 
  // {}
};

}  // namespace bigwham
#endif