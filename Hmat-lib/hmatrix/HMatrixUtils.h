#pragma once

#include <fstream>
#include <iostream>

#include <Hmat-lib/hmatrix/HMatrix.h>
#include <il/linearAlgebra.h>

namespace il {

template <typename T>
void toArray2D(const il::HMatrix<T> &H, il::spot_t s, il::io_t,
               il::Array2DEdit<T> E) {
  if (H.isFullRank(s)) {
    il::Array2DView<T> A = H.asFullRank(s);
    IL_EXPECT_MEDIUM(A.size(0) == E.size(0));
    IL_EXPECT_MEDIUM(A.size(1) == E.size(1));
    for (il::int_t i1 = 0; i1 < A.size(1); ++i1) {
      for (il::int_t i0 = 0; i0 < A.size(0); ++i0) {
        E(i0, i1) = A(i0, i1);
      }
    }
  } else if (H.isLowRank(s)) {
    il::Array2DView<T> A = H.asLowRankA(s);
    il::Array2DView<T> B = H.asLowRankB(s);
    IL_EXPECT_MEDIUM(A.size(0) == E.size(0));
    IL_EXPECT_MEDIUM(B.size(0) == E.size(1));
    IL_EXPECT_MEDIUM(A.size(1) == B.size(1));
    il::blas(static_cast<T>(1), A, B, il::Dot::Transpose, static_cast<T>(0),
             il::io, E);
  } else if (H.isHierarchical(s)) {
    const il::spot_t s00 = H.child(s, 0, 0);
    const il::spot_t s10 = H.child(s, 1, 0);
    const il::spot_t s01 = H.child(s, 0, 1);
    const il::spot_t s11 = H.child(s, 1, 1);
    il::toArray2D(H, s00, il::io, E.Edit(il::Range{0, H.size(0, s00)},
                                         il::Range{0, H.size(1, s00)}));
    il::toArray2D(
        H, s10, il::io,
        E.Edit(il::Range{H.size(0, s00), H.size(0, s00) + H.size(0, s10)},
               il::Range{0, H.size(1, s00)}));
    il::toArray2D(
        H, s01, il::io,
        E.Edit(il::Range{0, H.size(0, s00)},
               il::Range{H.size(1, s00), H.size(1, s00) + H.size(1, s01)}));
    il::toArray2D(
        H, s11, il::io,
        E.Edit(il::Range{H.size(0, s00), H.size(0, s00) + H.size(0, s10)},
               il::Range{H.size(1, s00), H.size(1, s00) + H.size(1, s01)}));
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
il::Array2D<double> toArray2D(const il::HMatrix<T> &H) {
  const il::int_t n0 = H.size(0);
  const il::int_t n1 = H.size(1);
  il::Array2D<T> A{n0, n1};
  il::Array2DEdit<T> E = A.Edit();
  il::spot_t s = H.root();
  toArray2D(H, s, il::io, E);
  return A;
}

template <typename T>
il::int_t nbElements(const il::HMatrix<T> &H, il::spot_t s) {
  if (H.isFullRank(s)) {
    il::Array2DView<T> A = H.asFullRank(s);
    return A.size(0) * A.size(1);
  } else if (H.isFullLu(s)) {
    il::Array2DView<T> A = H.asFullLu(s);
    return A.size(0) * A.size(1);
  } else if (H.isLowRank(s)) {
    il::Array2DView<T> A = H.asLowRankA(s);
    il::Array2DView<T> B = H.asLowRankB(s);
    return A.size(0) * A.size(1) + B.size(0) * B.size(1);
  } else if (H.isHierarchical(s)) {
    const il::spot_t s00 = H.child(s, 0, 0);
    const il::spot_t s10 = H.child(s, 1, 0);
    const il::spot_t s01 = H.child(s, 0, 1);
    const il::spot_t s11 = H.child(s, 1, 1);
    return il::nbElements(H, s00) + nbElements(H, s10) + nbElements(H, s01) +
        nbElements(H, s11);
  } else {
    IL_UNREACHABLE;
  }
  IL_UNREACHABLE;
  return 0;
}

template <typename T>
double compressionRatio(const il::HMatrix<T> &H) {
  const il::int_t n0 = H.size(0);
  const il::int_t n1 = H.size(1);
  return nbElements(H, H.root()) / static_cast<double>(n0 * n1);
}


/////////
template <typename T>
void hmatPatternF(const il::HMatrix<T> &H, il::spot_t s, il::ArrayView<int> x,
                  il::ArrayView<int> y, il::io_t, std::ostream &output) {
  // write to output filestream a HMatrix block pattern
  // the call of this function should be with x and y as vector containing -
  // 0-Ndof-1
  //
  // output : one row per block with the following format
  //   istart,jstart,iend,jend,Bool(0 for full, 1 for LR), number of memory entry
  IL_EXPECT_FAST(H.size(0, s) == y.size());
  IL_EXPECT_FAST(H.size(1, s) == x.size());

  // std::ofstream  output(filename);
  if (H.isFullRank(s)) {
    il::Array2DView<T> A = H.asFullRank(s);

    output << y[0] + 1 << "," << x[0] + 1 << "," << (y[0] + y.size()) << ","
           << (x[0] + x.size()) << ", " << 0 << ","
           <<  (A.size(0) * A.size(1)) << std::endl;
    return;
  } else if (H.isLowRank(s)) {
    il::Array2DView<T> A = H.asLowRankA(s);
    il::Array2DView<T> B = H.asLowRankB(s);
    output << y[0] + 1 << "," << x[0] + 1 << "," << (y[0] + y.size()) << ","
           << (x[0] + x.size()) << ", " << 1 << ","
           << (A.size(0) * A.size(1) + B.size(0) * B.size(1)) << std::endl;
    return;
  } else if (H.isHierarchical(s)) {
    const il::spot_t s00 = H.child(s, 0, 0);
    const il::spot_t s10 = H.child(s, 1, 0);
    const il::spot_t s01 = H.child(s, 0, 1);
    const il::spot_t s11 = H.child(s, 1, 1);
    const il::int_t n00 = H.size(0, s00);
    const il::int_t n10 = H.size(0, s10);
    const il::int_t n01 = H.size(1, s00);
    const il::int_t n11 = H.size(1, s01);
    il::ArrayView<int> x0 = x.view(il::Range{0, n01});
    il::ArrayView<int> x1 = x.view(il::Range{n01, n01 + n11});
    il::ArrayView<int> y0 = y.view(il::Range{0, n00});
    il::ArrayView<int> y1 = y.view(il::Range{n00, n00 + n10});
    hmatPatternF(H, s00, x0, y0, il::io_t(), output);
    hmatPatternF(H, s10, x0, y1, il::io_t(), output);
    hmatPatternF(H, s01, x1, y0, il::io_t(), output);
    hmatPatternF(H, s11, x1, y1, il::io_t(), output);
    return;
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void output_hmatPatternF(const il::HMatrix<T> &A,const std::string& filename) {
  // write the Hierarchical matrix A block pattern to a file
  // A :: Hierarchical matrix
  // filename :: string containing the path of the output file

  il::spot_t s = A.root();
  il::Array<int> xl{A.size(0), 0};
  for (il::int_t i = 0; i < A.size(0); i++) {
    xl[i] = static_cast<int>(i);
  }
  std::ofstream osh(filename);
  hmatPatternF(A, s, xl.view(), xl.view(), il::io, osh);
  osh.close();
}
//
//
/////////// OUTPUT H-PATTERN - matrix.....

template <typename T>
void hmatPattern(const il::HMatrix<T> &H, il::spot_t s,il::ArrayView<int> x,
                      il::ArrayView<int> y, il::io_t,
                      Array2D<il::int_t> &pattern) { //  //
  // write to output filestream a HMatrix block pattern
  // the call of this function should be with x and y as vector containing -
  // 0-Ndof-1
  //
  // output : one row per block with the following format
  //   spot index,istart,jstart

  il::int_t nr=pattern.size(0);
  il::int_t nc=pattern.size(1);

  IL_EXPECT_FAST(nr==3);

  if (H.isFullRank(s)) {
///    il::Array2DView<T> A = H.asFullRank(s);
  //  il::Array2DView<T> A = H.asFullRank(s);
    pattern.Resize(nr,nc+1);

    pattern(0,nc)=s.index;
    pattern(1,nc)=y[0] + 1;
    pattern(2,nc)=x[0] + 1;

//    output << y[0] + 1 << "," << x[0] + 1 << "," << (y[0] + y.size()) << ","
//           << (x[0] + x.size()) << ", " << 0 << ","
//           <<  (A.size(0) * A.size(1)) << std::endl;
    return;
  } else if (H.isLowRank(s)) {
    pattern.Resize(nr,nc+1);
    pattern(0,nc)=s.index;
    pattern(1,nc)=y[0] + 1;  // row.
    pattern(2,nc)=x[0] + 1;  // col.

//    output << y[0] + 1 << "," << x[0] + 1 << "," << (y[0] + y.size()) << ","
//           << (x[0] + x.size()) << ", " << 1 << ","
//           << (A.size(0) * A.size(1) + B.size(0) * B.size(1)) << std::endl;
    return;
  } else if (H.isHierarchical(s)) {
    const il::spot_t s00 = H.child(s, 0, 0);
    const il::spot_t s10 = H.child(s, 1, 0);
    const il::spot_t s01 = H.child(s, 0, 1);
    const il::spot_t s11 = H.child(s, 1, 1);
    const il::int_t n00 = H.size(0, s00);
    const il::int_t n10 = H.size(0, s10);
    const il::int_t n01 = H.size(1, s00);
    const il::int_t n11 = H.size(1, s01);
    il::ArrayView<int> x0 = x.view(il::Range{0, n01});
    il::ArrayView<int> x1 = x.view(il::Range{n01, n01 + n11});
    il::ArrayView<int> y0 = y.view(il::Range{0, n00});
    il::ArrayView<int> y1 = y.view(il::Range{n00, n00 + n10});
    hmatPattern(H, s00, x0, y0, il::io, pattern);
    hmatPattern(H, s10, x0, y1, il::io, pattern);
    hmatPattern(H, s01, x1, y0, il::io, pattern);
    hmatPattern(H, s11, x1, y1, il::io, pattern);
    return;
  } else {
    IL_UNREACHABLE;
  }
}
template <typename T>
il::Array2D<il::int_t> output_hmatPattern(const il::HMatrix<T> &A) {
  // write the Hierarchical matrix A block pattern to a matrix
  // A :: Hierarchical matrix
  // output:
  // and Array  with the spot index of all blocks

  il::spot_t s = A.root();
  il::Array<int> xl{A.size(0), 0};
  for (il::int_t i = 0; i < A.size(0); i++) {
    xl[i] = static_cast<int>(i);
  }

  il::Array2D<il::int_t> pattern{3,0};
  pattern.Reserve(3,A.size(0)/2); // here careful, reserve may be insufficient in some cases...

  hmatPattern(A, s, xl.view(), xl.view(), il::io, pattern);

  return pattern;

}

// output full rank entry only
template <typename T>
void hmatFullBlocks(const il::HMatrix<T> &H, il::spot_t s,
                      il::ArrayView<int> x, il::ArrayView<int> y, il::io_t,
                      il::Array<double> &val, il::Array2C<il::int_t> &pos) {
  // extract the full blocks entry
  // the call of this function should be with x and y as vector containing -
  // 0-Ndof-1
  //
  // output : pos the i,j of the full block entries
  //          val the corresponding values
  //
  IL_EXPECT_FAST(H.size(0, s) == y.size());
  IL_EXPECT_FAST(H.size(1, s) == x.size());

  il::int_t nr=pos.size(0);
  il::int_t nc=pos.size(1);

  IL_EXPECT_FAST(nc==2); // pos must be a 2 colummns array

  if (H.isFullRank(s)) {
    // std::cout << " full rank \n";
    il::Array2DView<T> A = H.asFullRank(s);
    il::int_t entb= A.size(0) * A.size(1);
    pos.Resize(nr+entb,2);
    val.Resize(nr+entb);
    il::int_t index=0;
    for (il::int_t j=0;j<A.size(1);j++){
      for (il::int_t i=0;i<A.size(0);i++){
        pos(nr+index,0)=i+y[0];
        pos(nr+index,1)=j+x[0];
        val[nr+index]=A(i,j);
        index++;
      }
    }
    return;
  } else if (H.isLowRank(s)) {
    //std::cout << " low rank \n";

    return;
  } else if (H.isHierarchical(s)) {
//    std::cout << " hierarchical \n";
    const il::spot_t s00 = H.child(s, 0, 0);
    const il::spot_t s10 = H.child(s, 1, 0);
    const il::spot_t s01 = H.child(s, 0, 1);
    const il::spot_t s11 = H.child(s, 1, 1);
    const il::int_t n00 = H.size(0, s00);
    const il::int_t n10 = H.size(0, s10);
    const il::int_t n01 = H.size(1, s00);
    const il::int_t n11 = H.size(1, s01);
    il::ArrayView<int> x0 = x.view(il::Range{0, n01});
    il::ArrayView<int> x1 = x.view(il::Range{n01, n01 + n11});
    il::ArrayView<int> y0 = y.view(il::Range{0, n00});
    il::ArrayView<int> y1 = y.view(il::Range{n00, n00 + n10});
    hmatFullBlocks(H, s00, x0, y0, il::io, val, pos);
    hmatFullBlocks(H, s10, x0, y1, il::io, val, pos);
    hmatFullBlocks(H, s01, x1, y0, il::io, val, pos);
    hmatFullBlocks(H, s11, x1, y1, il::io, val, pos);
    return;
  } else {
    IL_UNREACHABLE;
  }
}


template <typename T>
void output_hmatFullBlocks(const il::HMatrix<T> &A, il::Array<double>  &val,il::Array2C<il::int_t> &pos) {
  // write the Hierarchical matrix A block pattern to a matrix
  // A :: Hierarchical matrix
  // in/output:
  //   pos Array2C  (ne,2) containing the position of the full block entries
  //   val corresponding entry

  il::spot_t s = A.root();
  il::Array<int> xl{A.size(0), 0};
  for (il::int_t i = 0; i < A.size(0); i++) {
    xl[i] = static_cast<int>(i);
  }

  hmatFullBlocks(A, s, xl.view(), xl.view(), il::io, val, pos);
  std::cout <<" end of get full blocks - n FB entries " << val.size()  <<" \n";

}





}  // namespace il
