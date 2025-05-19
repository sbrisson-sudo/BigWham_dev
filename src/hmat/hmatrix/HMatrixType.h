//
// This file is part of BigWham.
//
// Created by Francois Fayard - 2018
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#pragma once

#ifndef BIGWHAM_HMATRIXTYPE_H
#define BIGWHAM_HMATRIXTYPE_H

namespace bigwham {

enum class HMatrixType { FullRank, LowRank, Hierarchical, FullLu };

struct SubHMatrix {
  il::Range range0;
  il::Range range1;
  bigwham::HMatrixType type;
};

}  // namespace bigwham

#endif