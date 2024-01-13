//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 15.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Dec. 2023 - new interface improvements


#include "bigwham_io_gen.h"
#include "bigwham_io_helper.h"

#include "core/be_mesh.h"
#include "core/bie_kernel.h"
#include "core/elastic_properties.h"


#include "elements/point.h"
#include "elements/rectangle.h"
#include "elements/segment.h"
#include "elements/triangle.h"

#include "elasticity/bie_elastostatic.h"
#include "elasticity/fullspace_iso_axisymmetry_flat_unidirectional/bie_elastostatic_axi3d0.h"
#include "elasticity/fullspace_iso_sp3d_segment/bie_elastostatic_sp3d.h"
#include "elasticity/fullspace_iso_2d_segment/bie_elastostatic_segment_0_influence.h"
#include "elasticity/fullspace_iso_2d_segment/bie_elastostatic_segment_1_influence.h"
#include "elasticity/fullspace_iso_3d_triangle/bie_elastostatic_triangle_0_influence.h"

/* -------------------------------------------------------------------------- */
using namespace bie;

// square matrix case
BigWhamIOGen::BigWhamIOGen(const std::vector<double> &coor, const std::vector<int> &conn, const std::string &kernel, const std::vector<double> &properties) {

    this->kernel_name_ = kernel;

    std::cout << " Now setting things for kernel ... " << kernel_name_
              << " with properties size " << properties.size() << "\n";

    switch (hash_djb2a(kernel_name_)) {
        case "2DP0"_sh: {
            IL_ASSERT(properties.size() == 2);
            ElasticProperties elas(properties[0], properties[1]);
            spatial_dimension_ = 2;
            dof_dimension_ = 2;
            flux_dimension_ = 3;
            int nvertices_per_elt_ = 2;
            using EltType = Segment<0>;
            mesh_ = bie::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                coor, conn);
            ker_obj_ = std::make_shared<bie::BieElastostatic<EltType, EltType, bie::ElasticKernelType::H>>(
                    elas, spatial_dimension_);
            using ObsType = Point<2>;
            ker_obs_u_=std::make_shared<bie::BieElastostatic<EltType, ObsType, bie::ElasticKernelType::T>>(
                    elas, spatial_dimension_);
            ker_obs_q_=std::make_shared<bie::BieElastostatic<EltType, ObsType, bie::ElasticKernelType::W>>(
                    elas, spatial_dimension_);
            break;
        }
       case "2DP1"_sh: {
            IL_ASSERT(properties.size() == 2);
            ElasticProperties elas(properties[0], properties[1]);

            spatial_dimension_ = 2;
            dof_dimension_ = 2;
            flux_dimension_ = 3;
            int nvertices_per_elt_ = 2;
            using EltType = bie::Segment<1>;
            mesh_ = bie::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                coor, conn);
            ker_obj_ = std::make_shared<bie::BieElastostatic<EltType, EltType, bie::ElasticKernelType::H>>(
                    elas, spatial_dimension_);
            using ObsType = Point<2>;
            ker_obs_q_=std::make_shared<bie::BieElastostatic<EltType,ObsType, bie::ElasticKernelType::W>>(
                   elas, spatial_dimension_);
            // still missing displacement for that element
            break;
        }
        case "S3DP0"_sh: {
            IL_ASSERT(properties.size() == 3);
            ElasticProperties elas(properties[0], properties[1]);
            spatial_dimension_ = 2;
            dof_dimension_ = 2;
            flux_dimension_ = 3;
            int nvertices_per_elt_ = 2;
            using EltType = bie::Segment<0>;
            mesh_ = bie::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                     coor, conn);
            ker_obj_ = std::make_shared<
                    bie::BieElastostaticSp3d<EltType, EltType, bie::ElasticKernelType::H>>(
                    elas, spatial_dimension_);

            il::Array<double> prop{1, properties[2]};
            ker_obj_->set_kernel_properties(prop);
            break;
        }
        case "3DT0"_sh: {
            IL_ASSERT(properties.size() == 2);
            ElasticProperties elas(properties[0], properties[1]);
            spatial_dimension_ = 3;
            dof_dimension_ = 3;
            flux_dimension_ = 6; // 6 stress components
            int nvertices_per_elt_ = 3;
            using EltType = bie::Triangle<0>;
            mesh_ = bie::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                coor, conn);
            ker_obj_ = std::make_shared<
                    bie::BieElastostatic<EltType, EltType, bie::ElasticKernelType::H>>(
                    elas, spatial_dimension_);
            using ObsType = Point<3>;
            ker_obs_u_=std::make_shared<bie::BieElastostatic<EltType, ObsType, bie::ElasticKernelType::T>>(
                    elas, spatial_dimension_);
            ker_obs_q_=std::make_shared<bie::BieElastostatic<EltType,ObsType, bie::ElasticKernelType::W>>(
                    elas, spatial_dimension_);
            break;
        }
        case "Axi3DP0"_sh: {
            IL_ASSERT(properties.size() == 2);
            ElasticProperties elas(properties[0], properties[1]);
            spatial_dimension_ = 2;
            dof_dimension_ = 2;
            int nvertices_per_elt_ = 2;
            using EltType = bie::Segment<0>;
            mesh_ = bie::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                coor, conn);
            ker_obj_ = std::make_shared<bie::ElasticAxiSymmRingKernel>(
                    elas, spatial_dimension_); // todo - change this API to the new one.... should be std::make_shared<
           // bie::BieElastostatic<EltType, EltType, bie::ElasticKernelType::H>>
            break;
        }
        default: {
            std::cout << "wrong inputs -abort \n";
            il::abort();
        }
    }
    mesh_src_ = mesh_;
    mesh_rec_ = mesh_;

};

// rectangular Hmatrix
BigWhamIOGen::BigWhamIOGen(const std::vector<double> &coor_src,
             const std::vector<int> &conn_src,
             const std::vector<double> &coor_rec,
             const std::vector<int> &conn_rec, const std::string &kernel,
             const std::vector<double> &properties){
    this->kernel_name_ = kernel;
    std::cout << " Now setting things for kernel ... " << kernel_name_
              << " with properties size " << properties.size() << "\n";

    switch (hash_djb2a(kernel_name_)) {
        case "3DT0-3DT0-H"_sh: {
            ElasticProperties elas(properties[0], properties[1]);
            dof_dimension_ = 3;
            spatial_dimension_ = 3;
            flux_dimension_ = 6 ;
            using src_elem = Triangle<0>;
            using rec_elem = Triangle<0>;
            mesh_src_ = bie::CreateMeshFromVect<src_elem>(
                    spatial_dimension_, /* num vertices */ 3, coor_src, conn_src);
            mesh_rec_ = bie::CreateMeshFromVect<rec_elem>(
                    spatial_dimension_, /* num vertices */ 3, coor_rec, conn_rec);
            ker_obj_ = std::make_shared<
                    BieElastostatic<src_elem, rec_elem, ElasticKernelType::H>>(
                    elas, spatial_dimension_);
            break;
        }
        case "2DS0-2DP-T"_sh: {
            // 2D Segment0 and 2D Point
            ElasticProperties elas(properties[0], properties[1]);
            dof_dimension_ = 2;
            spatial_dimension_ = 2;
            flux_dimension_ = 3;
            using src_elem = Segment<0>;
            using rec_elem = Point<2>;
            mesh_src_ = bie::CreateMeshFromVect<src_elem>(
                    spatial_dimension_, /* num vertices */ 2, coor_src, conn_src);
            mesh_rec_ = bie::CreateMeshFromVect<rec_elem>(
                    spatial_dimension_, /* num vertices */ 1, coor_rec, conn_rec);
            ker_obj_ = std::make_shared<
                    BieElastostatic<src_elem, rec_elem, ElasticKernelType::T>>(
                    elas, spatial_dimension_);
            break;
        }
        default: {
            std::cout << "wrong inputs -abort \n";
            il::abort();
        }
    }
    mesh_ = mesh_src_;
}

// Method for pattern Construction
// API to the Hierarchical representation function
void BigWhamIOGen::BuildPattern(const int max_leaf_size, const double eta) {

    max_leaf_size_ = max_leaf_size;
    eta_ = eta;
    il::Timer tt;
    std::cout << "--------------------\n";
    std::cout << "Hierarchical representation creation ...\n";
    tt.Start();
    this->hr_ = HRepresentationRectangularMatrix(mesh_src_, mesh_rec_,max_leaf_size_, eta_);
    is_pattern_built_= true;
    tt.Stop();
    std::cout << "Hierarchical representation complete.\n";
}

// construct Hierarchical matrix
void BigWhamIOGen::BuildHierarchicalMatrix(const int max_leaf_size, const double eta, const double eps_aca) {

  epsilon_aca_ = eps_aca;
  if (! is_pattern_built_){
      BigWhamIOGen::BuildPattern(max_leaf_size, eta);
  }

  il::Timer tt;
  std::cout << "--------------------\n";
  std::cout << "Populating Hierarchical matrix ...\n";
  tt.Start();
  bie::BieMatrixGenerator<double> M(mesh_src_, mesh_rec_, ker_obj_, hr_);
  hmat_ = std::make_shared<Hmat<double>>(M, epsilon_aca_);
  tt.Stop();
  hmat_time_ = tt.time();
  if (hmat_->isBuilt()) {
    is_built_ = true;
    IL_EXPECT_FAST(dof_dimension_ == hmat_->dofDimension());
    std::cout << "Hierarchical matrix construction complete.\n";
    double test_cr = hmat_->compressionRatio();
    std::cout <<  "Compression Ratio = " << test_cr << ", eps_aca = " << epsilon_aca_
              << ", eta = " << eta_ << "\n";
    std::cout << "Hierarchical matrix  construction time = :  " << hmat_time_ << "\n";
  } else {
    is_built_ = false;
  }
  std::cout << "--------------------\n";

}
/* -------------------------------------------------------------------------- */

std::vector<long> BigWhamIOGen::GetHPattern() const {
  // API function to output the hmatrix pattern
  //  as flattened list via a pointer
  //  the numberofblocks is also returned (by reference)
  //
  //  the pattern matrix is formatted as
  // row = 1 block : i_begin,j_begin, i_end,j_end,FLAG,entry_size
  // with FLAG=0 for full rank and FLAG=1 for low rank
  //
  // we output a flatten row-major order std::vector

  IL_EXPECT_FAST(is_pattern_built_ || is_built_);

  bie::HPattern pattern =hr_->pattern_;
  if (is_built_) {
      pattern =hmat_->pattern();
  }

  long numberofblocks = pattern.n_B;
  long len = 6 * numberofblocks;
  std::cout << "number of blocks " << numberofblocks << "\n";

  std::vector<long> patternlist(len, 0);

  int index = 0;
  //  starts with full rank
  for (il::int_t j = 0; j < pattern.n_FRB; j++) {
    // check is low rank or not
    //  il::Array2DView<double> A = hmat_.asFullRank(s);
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
  // then low ranks (only valid values if the Hmat has been built !)
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
/* --------------------------------------------------------------------------*/

void BigWhamIOGen::GetFullBlocks(il::Array<double> &val_list,
                                 il::Array<int> &pos_list) const {
  // return the full dense block entries of the hmat as
  // flattened lists
  // val_list(i) = H(pos_list(2*i),pos_list(2*i+1));
  // output in the original dof state (accounting for the permutation)

  IL_EXPECT_FAST(is_built_);

  il::Array<double> values{};
  il::Array<int> positions{};
  hmat_->fullBlocksOriginal(il::io, values, positions);
  //    std::cout << " checking values size" << values.size() <<  "\n";
  val_list.Resize(values.size());
  for (il::int_t i = 0; i < values.size(); i++) {
    val_list[i] = values[i];
  }
  pos_list.Resize(positions.size());
  for (il::int_t i = 0; i < positions.size(); i++) {
    pos_list[i] = positions[i];
  }
 // std::cout << "number of entries " << val_list.size() << " - "
 //           << pos_list.size() << "\n";
 // std::cout << " End of Bigwhamio getFullBlocks \n";
}
/* --------------------------------------------------------------------------*/

void BigWhamIOGen::GetFullBlocks(std::vector<double> &val_list,
                                 std::vector<int> &pos_list) const {
  // return the full dense block entries of the hmat as
  // flattened lists
  // val_list(i) = H(pos_list(2*i),pos_list(2*i+1));
  // output in the original dof state (accounting for the permutation)

  IL_EXPECT_FAST(is_built_);

  il::Array<double> values{};
  il::Array<int> positions{};
  hmat_->fullBlocksOriginal(il::io, values, positions);
  //    std::cout << " checking values size" << values.size() <<  "\n";
  val_list.reserve(values.size());
  for (il::int_t i = 0; i < values.size(); i++) {
    val_list.push_back(values[i]);
  }
  pos_list.reserve(positions.size());
  for (il::int_t i = 0; i < positions.size(); i++) {
    pos_list.push_back(positions[i]);
  }
  //std::cout << "number of entries " << val_list.size() << " - "
  //          << pos_list.size() << "\n";
  //std::cout << " End of Bigwhamio getFullBlocks \n";
}
/* --------------------------------------------------------------------------*/

void BigWhamIOGen::GetDiagonal(std::vector<double> &val_list) const {
  // return the diagonal of the h-matrix
  // output in the original dof state (accounting for the permutation)

  IL_EXPECT_FAST(is_built_);
  val_list = hmat_->diagonalOriginal();

  std::cout << " End of Bigwhamio getDiagonal() \n";
}
/* -------------------------------------------------------------------------- */

std::vector<double> BigWhamIOGen::MatVec(const std::vector<double> &x) const {
  // in the original / natural ordering
  IL_EXPECT_FAST(this->is_built_);
  IL_EXPECT_FAST(hmat_->size(1) == x.size());
  std::vector<double> y = hmat_->matvecOriginal(x);
  return y;
}
/* -------------------------------------------------------------------------- */

il::Array<double> BigWhamIOGen::MatVec(il::ArrayView<double> x) const {
  // in the original / natural ordering
  IL_EXPECT_FAST(this->is_built_);
  IL_EXPECT_FAST(hmat_->size(1) == x.size());
  auto y = hmat_->matvecOriginal(x);
  return y;
}
/* -------------------------------------------------------------------------- */
il::Array<double> BigWhamIOGen::MatVecPerm(il::ArrayView<double> x) const {
  // in the permutted state.
  IL_EXPECT_FAST(this->is_built_);
  IL_EXPECT_FAST(hmat_->size(1) == x.size());
  auto y = hmat_->matvec(x);
  return y;
}
/* -------------------------------------------------------------------------- */

std::vector<double>
BigWhamIOGen::MatVecPerm(const std::vector<double> &x) const {
  // in the permutted state.
  IL_EXPECT_FAST(this->is_built_);
  IL_EXPECT_FAST(hmat_->size(1) == x.size());
  std::vector<double> y = hmat_->matvec(x);
  return y;
}
/* -------------------------------------------------------------------------- */

il::Array<double>
BigWhamIOGen::ConvertToGlobal(il::ArrayView<double> x_local) const {
  // Input: x in original state (not permutted)
  // Output: in original state (not permutted)
  return mesh_rec_->ConvertToGlobal(x_local);
}
/* -------------------------------------------------------------------------- */
il::Array<double>
BigWhamIOGen::ConvertToLocal(il::ArrayView<double> x_global) const {
  // Input: x in original state (not permutted)
  // Output: in original state (not permutted)
  return mesh_src_->ConvertToLocal(x_global);
}
/* -------------------------------------------------------------------------- */

std::vector<double>
BigWhamIOGen::ConvertToGlobal(const std::vector<double> &x_local) const {
  // Input: x in original state (not permutted)
  // Output: in original state (not permutted)
  return mesh_rec_->ConvertToGlobal(x_local);
}
/* -------------------------------------------------------------------------- */

std::vector<double>
BigWhamIOGen::ConvertToLocal(const std::vector<double> &x_global) const {
  // Input: x in original state (not permutted)
  // Output: in original state (not permutted)
  return mesh_src_->ConvertToLocal(x_global);
}
/* -------------------------------------------------------------------------- */

std::vector<double> BigWhamIOGen::GetCollocationPoints() const {
  IL_EXPECT_FAST(is_built_);

  auto col_pts = mesh_rec_->collocation_points();  // this should be the receiver mesh for  generality
  IL_EXPECT_FAST(col_pts.size(1) == spatial_dimension_);
  il::int_t npoints = col_pts.size(0);
  std::vector<double> flat_col;
  flat_col.assign(npoints * spatial_dimension_, 0.);
  int index = 0;
  for (il::int_t i = 0; i < col_pts.size(0); i++) {
    for (il::int_t j = 0; j < col_pts.size(1); j++) {
      flat_col[index] = col_pts(i, j);
      index++;
    }
  }
  return flat_col;
}
/* -------------------------------------------------------------------------- */

std::vector<long> BigWhamIOGen::GetPermutation() const {
  IL_EXPECT_FAST(is_built_);
  std::vector<long> permut;
  permut.assign(hr_->permutation_1_.size(), 0);
  for (il::int_t i = 0; i < hr_->permutation_1_.size(); i++) {
    permut[i] = hr_->permutation_1_[i];
  }
  return permut;
}
/* -------------------------------------------------------------------------- */

il::Array<double> BigWhamIOGen::ComputePotentials(const std::vector<double> &coor_obs , const il::ArrayView<double> sol_local) const
{ // return the potential in the global system of coordinates, the solution on the source mesh being defined in the local eletment system of the source mesh
    IL_EXPECT_FAST(sol_local.size() == mesh_src_->num_collocation_points() * dof_dimension_);
    IL_EXPECT_FAST(coor_obs.size() % spatial_dimension_ == 0);
    if (this->ker_obs_u_ == nullptr) {
        std::cout << "Kernel Not Implemented !\n";
    }
    IL_EXPECT_FAST(this->ker_obs_u_!= nullptr);

    il::int_t npts = coor_obs.size() / spatial_dimension_;
    std::shared_ptr<bie::Mesh> mesh_obs;
    std::vector<int> conn_obs(npts,0);
    for (int i =0 ; i<npts;i++){
        // dummy connectivity for a mesh of points.
        conn_obs[i]=i;
    }
    IL_EXPECT_FAST(spatial_dimension_==2 || spatial_dimension_==3);

    switch (spatial_dimension_) {
        case 2 : {
            mesh_obs = bie::CreateUniqueMeshFromVect<bie::Point<2>>(spatial_dimension_, 1, coor_obs, conn_obs);
        }
        case 3: {
            mesh_obs = bie::CreateUniqueMeshFromVect<bie::Point<3>>(spatial_dimension_, 1, coor_obs, conn_obs);
        }
    }

    il::Array<double> obs_potential{npts*dof_dimension_,0.};

// loop on source collocation points
#pragma omp for
    for (il::int_t i=0;i<mesh_src_->num_collocation_points();i++){

        il::int_t e_i = this->mesh_src_->GetElementId(i);
        il::int_t is_l = this->mesh_src_->GetElementCollocationId(i);
        auto source_element = this->mesh_src_->GetElement(e_i);

        // local solution (assumed in original ordering....)
        il::Array<double> elt_solu{dof_dimension_};
        for (il::int_t k=0;k<dof_dimension_;k++){
            elt_solu[k]=sol_local[i*dof_dimension_+k];
        }

        // loop on obs points mesh
        for (il::int_t j_obs =0;j_obs<mesh_obs->num_collocation_points();j_obs++){
            il::int_t e_j_r = mesh_obs->GetElementId(j_obs);
            auto receiver_element = mesh_obs->GetElement(e_j_r);
            il::int_t ir_l = mesh_obs->GetElementCollocationId(j_obs);
            std::vector<double> st = this->ker_obs_u_->influence(*source_element, is_l, *receiver_element,ir_l);
            for (il::int_t j=0;j<dof_dimension_;j++){
                for(il::int_t k=0;k<dof_dimension_;k++){
                    obs_potential[j_obs*dof_dimension_+j]+=st[k*dof_dimension_+j]*elt_solu[k];
                }
            }
        }
    }

return obs_potential;
}


il::Array<double> BigWhamIOGen::ComputeFluxes(const std::vector<double> &coor_obs , const il::ArrayView<double> sol_local) const
{
// return the potential in the global system of coordinates, the solution on the source mesh being defined in the local eletment system of the source mesh
    IL_EXPECT_FAST(sol_local.size() == mesh_src_->num_collocation_points() * dof_dimension_);
    IL_EXPECT_FAST(coor_obs.size() % spatial_dimension_ == 0);
    if (this->ker_obs_q_ == nullptr) {
        std::cout << "Kernel Not Implemented !\n";
    }
    IL_EXPECT_FAST(this->ker_obs_q_!= nullptr);

    il::int_t npts = coor_obs.size() / spatial_dimension_;

    std::shared_ptr<bie::Mesh> mesh_obs;
    std::vector<int> conn_obs(npts,0);
    for (int i =0 ; i<npts;i++){
        // dummy connectivity for a mesh of points.
        conn_obs[i]=i;
    }
    IL_EXPECT_FAST(spatial_dimension_==2 || spatial_dimension_==3);

    switch (spatial_dimension_) {
        case 2 : {
            mesh_obs = bie::CreateUniqueMeshFromVect<bie::Point<2>>(spatial_dimension_, 1, coor_obs, conn_obs);
        }
        case 3: {
            mesh_obs = bie::CreateUniqueMeshFromVect<bie::Point<3>>(spatial_dimension_, 1, coor_obs, conn_obs);
        }
    }

    il::Array<double> obs_flux{npts*flux_dimension_,0.};

// loop on source collocation points
#pragma omp for
    for (il::int_t i = 0; i < mesh_src_->num_collocation_points(); i++) {

        il::int_t e_i = this->mesh_src_->GetElementId(i);
        il::int_t is_l = this->mesh_src_->GetElementCollocationId(i);

        auto source_element = this->mesh_src_->GetElement(e_i);
        // local solution (assumed in original ordering....)
        il::Array<double> elt_solu{dof_dimension_}, e_potential{dof_dimension_, 0.};
        for (il::int_t k = 0; k < dof_dimension_; k++) {
            elt_solu[k] = sol_local[i * dof_dimension_ + k];
        }

        // loop on obs points mesh
        for (il::int_t j_obs = 0; j_obs < mesh_obs->num_collocation_points(); j_obs++) {
            il::int_t e_i_r  = mesh_obs->GetElementId(j_obs);
            auto receiver_element = mesh_obs->GetElement(e_i_r);
            il::int_t ir_l = mesh_obs->GetElementCollocationId(j_obs);
            std::vector<double> st = this->ker_obs_q_->influence(*source_element, is_l, *receiver_element, ir_l);
            for (il::int_t j = 0; j < flux_dimension_; j++) {
                for (il::int_t k = 0; k < dof_dimension_; k++) {
                    obs_flux[j_obs * flux_dimension_ + j] += st[flux_dimension_ * k + j] * elt_solu[k];
                }
            }
        }
    }
    return obs_flux;
}
