#pragma once

namespace il {

enum class HMatrixType { FullRank, LowRank, Hierarchical, FullLu };

struct SubHMatrix {
  il::Range range0;
  il::Range range1;
  il::HMatrixType type;
};

}  // namespace il
