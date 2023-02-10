//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 15.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications ::January 31. 2023 - cleaning up and using the new code interface.

#pragma once

#include <iostream>
#include <string>
#include <string_view>

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/Dynamic.h>
#include <il/Map.h>

#include <hmat/cluster/cluster.h>
#include <hmat/hmatrix/Hmat.h>

#include <src/core/BEMesh.h>
#include <src/core/BoundaryElement.h>
#include <src/core/elements/Segment.h>
#include <src/core/elements/Triangle.h>
#include <src/core/elements/Rectangle.h>

#include <src/core/BIE_Kernel.h>
#include <src/elasticity/BIE_elastostatic.h>
#include <src/core/SquareMatrixGenerator.h>

#include <src/elasticity/2d/BIE_elastostatic_segment_0_impls.h>
#include <src/elasticity/2d/BIE_elastostatic_segment_1_impls.h>
#include <src/elasticity/3d/BIE_elastostatic_triangle_0_impls.h>

#include <src/core/ElasticProperties.h>


// utilities for switch with string in C++17
// https://learnmoderncpp.com/2020/06/01/strings-as-switch-case-labels/
inline constexpr auto hash_djb2a(const std::string_view sv) {
  unsigned long hash{5381};
  for (unsigned char c : sv) {
    hash = ((hash << 5) + hash) ^ c;
  }
  return hash;
}

inline constexpr auto operator"" _sh(const char *str, size_t len) {
  return hash_djb2a(std::string_view{str, len});
}
//////////////////////////// utility for mesh object creation from std::vector
template<class El>
bie::BEMesh<El> createMeshFromVect(int spatial_dimension,int n_vertex_elt,const
 std::vector<double>& coor, const std::vector<int>& conn){
    il::int_t npoints = coor.size() / spatial_dimension;
    il::int_t nelts = conn.size() / spatial_dimension;
    il::Array2D<double> Coor{npoints, spatial_dimension, 0.};  //
    il::Array2D<il::int_t> Conn{nelts, n_vertex_elt, 0};
    // populate mesh  ...
    int index = 0;
    for (il::int_t i = 0; i < Coor.size(0); i++) {
        for (il::int_t j = 0; j < Coor.size(1); j++) {
            Coor(i, j) = coor[index];
            index++;
        }
    }
    index = 0;
    for (il::int_t i = 0; i < Conn.size(0); i++) {
        for (il::int_t j = 0; j < Conn.size(1); j++) {
            Conn(i, j) = conn[index];
            index++;
        }
    }
    bie::BEMesh<El> mesh(Coor,Conn);
    std::cout << "in create mesh - done "  << "\n";
    return mesh;
}
//////////////////////////// the 'infamous' Bigwhamio class
class Bigwhamio {
private:
  bie::Hmat<double> h_{}; // the  Hmat object

  il::Array<il::int_t>   permutation_; // permutation list of the collocation points
  il::Array2D<double> collocationPoints_; //  collocation points coordinates ?

  int dimension_;     // spatial dimension
  int dof_dimension_; // number of dof per nodes / collocation points

  bool isBuilt_; // True if the class instance is built

  // H-matrix parameters
  int max_leaf_size_;
  double eta_;
  double epsilon_aca_;
  // kernel
  std::string kernel_;

  // statistics
  double h_representation_time_;
  double hmat_time_;

public:
  //---------------------------------------------------------------------------
  Bigwhamio() { // default constructor
    dimension_ = 0;
    dof_dimension_ = 0;
    isBuilt_ = false;
    eta_ = 0.0;
    epsilon_aca_ = 0.001;
    max_leaf_size_ = 100;
    kernel_ = "none";
    h_representation_time_ = 0.;
    hmat_time_ = 0.;
  };

  ~Bigwhamio() = default;

  void set(const std::vector<double> &coor, const std::vector<int> &conn,
           const std::string &kernel, const std::vector<double> &properties,
           const int max_leaf_size, const double eta, const double eps_aca) {
    // coor and conn are assumed to be passed in row-major storage format
    kernel_ = kernel;
    // switch depending on Kernels for mesh building
    max_leaf_size_ = max_leaf_size;
    eta_ = eta;
    epsilon_aca_ = eps_aca;
    bie::ElasticProperties elas(properties[0], properties[1]);
    std::cout << " Now setting things for kernel ... " << kernel_ <<" prop si" << properties.size()<<  "\n";
    il::Timer tt;
    if (kernel_ == "S3DP0") {
          IL_ASSERT(properties.size() == 3);
      } else  {
          IL_ASSERT(properties.size() == 2);
    }

   switch(hash_djb2a(kernel_)) {
          case "S3DP0"_sh:{
              dimension_ = 2;
              using EltType = bie::Segment<0>;
              int  nvertices_per_elt_=dimension_;
              std::cout << " in S3DP0 " << coor.size() << "-" << conn.size();
              bie::BEMesh<EltType> mesh=createMeshFromVect<EltType>(dimension_,nvertices_per_elt_,coor, conn);
              auto tes=mesh.getVertices(0);
              mesh.setCurrentElement(0);
              std::cout << " ELement vertices " << tes(0,0) << "- " <<tes(0,1) << tes(1,0) << "- " <<tes(1,1)<<"\n";
              collocationPoints_=mesh.getCollocationPoints(); // be careful returning it in original ordering.  note this is only for the output function getCollocationPoints ... could be deleted possibly
              std::cout << " mesh - done - n_elts" << mesh.numberOfElts() << "\n";
              //tt.Start();
              bie::HRepresentation hr=bie::h_representation_square_matrix(mesh,max_leaf_size_,eta_);
              std::cout << " pattern created \n";
              //tt.Stop();
              h_representation_time_=tt.time();tt.Reset();tt.Start();
              const auto ker_type = bie::ElasticKernelType::H;
              bie::BIE_elastostatic_sp3d<EltType,EltType,ker_type>  ker(elas,dimension_);
              il::Array<double> prop{1,properties[2]}; // for the SP3D0
              ker.setKernelProperties(prop);
              bie::SquareMatrixGenerator<double,EltType,bie::BIE_elastostatic_sp3d<EltType,EltType,ker_type>> M(mesh,ker,hr.permutation_0_);
              h_.toHmat(M,hr,epsilon_aca_);
              tt.Stop();
              permutation_=hr.permutation_0_;
              hmat_time_=tt.time();
              break;
          }
          case "2DP1"_sh :{
              dimension_ = 2;
              using EltType = bie::Segment<1>;
              int  nvertices_per_elt_=dimension_;
              bie::BEMesh<EltType> mesh=createMeshFromVect<EltType>(dimension_,nvertices_per_elt_,coor, conn);
              tt.Start();
              bie::HRepresentation hr=bie::h_representation_square_matrix(mesh,max_leaf_size,eta);
              tt.Stop();
              std::cout << " pattern created \n";
              collocationPoints_=mesh.getCollocationPoints(); // be careful returning it in original ordering.  note this is only for the output function getCollocationPoints ... could be deleted possibly
              h_representation_time_=tt.time();tt.Reset();tt.Start();
              const auto ker_type = bie::ElasticKernelType::H;
              bie::BIE_elastostatic<EltType,EltType,ker_type> ker(elas,dimension_);
              bie::SquareMatrixGenerator<double,EltType,bie::BIE_elastostatic<EltType,EltType,ker_type>> M(mesh,ker,hr.permutation_0_);
              h_.toHmat(M,hr,epsilon_aca_);
              tt.Stop();
              permutation_=hr.permutation_0_;
              hmat_time_=tt.time();
              break;
          }
          case "3DT0"_sh :{
              dimension_ = 3;
              int  nvertices_per_elt_=3;
              using EltType = bie::Triangle<0>;
              bie::BEMesh<EltType> mesh=createMeshFromVect<EltType>(dimension_,nvertices_per_elt_,coor, conn);
              tt.Start();
              bie::HRepresentation hr=bie::h_representation_square_matrix(mesh,max_leaf_size,eta);
              tt.Stop();
              collocationPoints_=mesh.getCollocationPoints(); // be careful returning it in original ordering.  note this is only for the output function getCollocationPoints ... could be deleted possibly
              h_representation_time_=tt.time();tt.Reset();tt.Start();
              const auto ker_type = bie::ElasticKernelType::H;
              bie::BIE_elastostatic<EltType,EltType,ker_type>  ker(elas,dimension_);
              bie::SquareMatrixGenerator<double,EltType,bie::BIE_elastostatic<EltType,EltType,ker_type>> M(mesh,ker,hr.permutation_0_);
              h_.toHmat(M,hr,epsilon_aca_);
              tt.Stop();
              permutation_=hr.permutation_0_;
              hmat_time_=tt.time();
              break;
          }
          default:
          {
              std::cout << "wrong inputs -abort \n";
              il::abort();
          }
   }
    if (h_.isBuilt()) {
      isBuilt_ = true;
      dof_dimension_=h_.dofDimension();
      std::cout << "HMAT --> built \n";
      double test_cr=h_.compressionRatio();
      std::cout << "H mat set " << " CR = " <<  test_cr << " eps_aca "     << epsilon_aca_ << " eta " << eta_ << "\n";
      //std::cout << "H-mat construction time = :  " << hmat_time_ << "\n";
    } else {
      isBuilt_ = false;
    }

    std::cout << "end of set() of bigwhamio object set \n";
  }

  bool isBuilt() { return isBuilt_; };

  void hmatDestructor() {
    // this function will free the memory and set the hmat obj to its initial
    // status prior to initialization this will avoid ownership specifications
    // at binding time
    this->h_.hmatMemFree();
  };

  //---------------------------------------------------------------------------
  //  get and other methods below
  double getHmatTime() const { return hmat_time_; };
  double getPatternTime() const { return h_representation_time_; };

  std::vector<double> getCollocationPoints() {
    IL_EXPECT_FAST(isBuilt_);
    std::cout << "beginning of getCollocationPoints bigwham \n";
    std::cout << " Spatial dim :" << dimension_ << " collocation dim size :" << collocationPoints_.size(1) << "\n";
    std::cout << " collocation npoints :" << collocationPoints_.size(0) << "\n";
    std::cout << "coll points dim " << collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";

    IL_EXPECT_FAST(collocationPoints_.size(1) == dimension_);

    il::int_t npoints = collocationPoints_.size(0);
    std::vector<double> flat_col;
    flat_col.assign(npoints * dimension_, 0.);
    int index = 0;
    for (il::int_t i = 0; i < collocationPoints_.size(0); i++) {
      for (il::int_t j = 0; j < collocationPoints_.size(1); j++) {
        flat_col[index] = collocationPoints_(i, j);
        index++;
      }
    }
    std::cout << "end of getCollocationPoints bigwham \n";
    return flat_col;
  };

  //---------------------------------------------------------------------------
  std::vector<long> getPermutation() {
    IL_EXPECT_FAST(isBuilt_);
    std::cout << "coll points dim " << collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";
    std::vector<long> permut;
    permut.assign(permutation_.size(), 0);
    for (il::int_t i = 0; i < permutation_.size(); i++) {
      permut[i] = permutation_[i];
    }
    return permut;
  }
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  double getCompressionRatio() {
    IL_EXPECT_FAST(isBuilt_); return h_.compressionRatio();
  }

  std::string getKernel() { return kernel_; }

  int getSpatialDimension() const { return dimension_; }

  int getProblemDimension() const { return dof_dimension_; }

  long matrixSize(int k) { return h_.size(k); };

  //---------------------------------------------------------------------------
  std::vector<long> getHpattern() {
    // API function to output the hmatrix pattern
    //  as flattened list via a pointer
    //  the numberofblocks is also returned (by reference)
    //
    //  the pattern matrix is formatted as
    // row = 1 block : i_begin,j_begin, i_end,j_end,FLAG,entry_size
    // with FLAG=0 for full rank and FLAG=1 for low rank
    //
    // we output a flatten row-major order std::vector

    IL_EXPECT_FAST(isBuilt_);

    bie::HPattern pattern = h_.pattern();

    long numberofblocks = pattern.n_B;
    long len = 6 * numberofblocks;
    std::cout << "number of blocks " << numberofblocks << "\n";

    std::vector<long> patternlist(len, 0);

    int index = 0;
    //  starts with full rank
    for (il::int_t j = 0; j < pattern.n_FRB; j++) {
      // check is low rank or not
      //  il::Array2DView<double> A = h_.asFullRank(s);
      patternlist[index++] = pattern.FRB_pattern(1, j);
      patternlist[index++] = pattern.FRB_pattern(2, j);
      patternlist[index++] = pattern.FRB_pattern(3, j);
      patternlist[index++] = pattern.FRB_pattern(4, j);
      patternlist[index++] = 0;
      patternlist[index++] =
          (pattern.FRB_pattern(4, j) - pattern.FRB_pattern(2, j)) *
          (pattern.FRB_pattern(3, j) -
           pattern.FRB_pattern(1, j)); // size of that sub blocks
    }
    // then low ranks
    for (il::int_t j = 0; j < pattern.n_LRB; j++) {

      patternlist[index++] = pattern.LRB_pattern(1, j);
      patternlist[index++] = pattern.LRB_pattern(2, j);
      patternlist[index++] = pattern.LRB_pattern(3, j);
      patternlist[index++] = pattern.LRB_pattern(4, j);
      patternlist[index++] = 1;
      patternlist[index++] = pattern.LRB_pattern(5, j); // the rank
    }
    // return a row major flatten vector
    return patternlist;
  }
  //---------------------------------------------------------------------------
  void getFullBlocks(std::vector<double> &val_list,std::vector<int> &pos_list) {
    // return the full dense block entries of the hmat as
    // flattened lists
    // val_list(i) = H(pos_list(2*i),pos_list(2*i+1));
    // output in the original dof state (accounting for the permutation)

    IL_EXPECT_FAST(isBuilt_);

    il::Array<double> values{};
    il::Array<int> positions{};
    h_.fullBlocksOriginal(permutation_, il::io, values, positions);
    //    std::cout << " checking values size" << values.size() <<  "\n";
    val_list.reserve(values.size());
    for (il::int_t i = 0; i < values.size(); i++) {
      val_list.push_back(values[i]);
    }
    pos_list.reserve(positions.size());
    for (il::int_t i = 0; i < positions.size(); i++) {
      pos_list.push_back(positions[i]);
    }
    std::cout << "number of entries " << val_list.size() << " - "     << pos_list.size() << "\n";
    std::cout << " End of Bigwhamio getFullBlocks \n";
  }

//---------------------------------------------------------------------------
    void getDiagonal(std::vector<double> &val_list) {
        // return the diagonal of the h-matrix
        // output in the original dof state (accounting for the permutation)

        IL_EXPECT_FAST(isBuilt_);
        val_list=h_.diagonalOriginal(permutation_);

        std::cout << " End of Bigwhamio getDiagonal() \n";
    }

  // ---------------------------------------------------------------------------
  std::vector<double> matvect(const std::vector<double> &x) {
      // in the original / natural ordering
    IL_EXPECT_FAST(this->isBuilt_);
    IL_EXPECT_FAST(h_.size(0) == h_.size(1));
    IL_EXPECT_FAST(h_.size(1) == x.size());
    std::vector<double> y = h_.matvecOriginal(permutation_, x);
    return y;
  }

  // ---------------------------------------------------------------------------
  std::vector<double> hdotProductInPermutted(const std::vector<double> &x) {
      // in the permutted state.
    IL_EXPECT_FAST(this->isBuilt_);
    IL_EXPECT_FAST(h_.size(0) == h_.size(1));
    IL_EXPECT_FAST(h_.size(1) == x.size());
    std::vector<double> y = h_.matvec(x);
    return y;
  }




}; // end class bigwhamio
