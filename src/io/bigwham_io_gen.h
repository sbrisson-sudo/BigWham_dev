//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 15.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Dec. 2023 - new interface improvements

#pragma once

#ifndef BIGWHAM_IO_GEN_H
#define BIGWHAM_IO_GEN_H

#ifdef IL_OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <memory>

#include "hmat/hierarchical_representation.h"
#include "hmat/bie_matrix_generator.h"
#include "hmat/hmatrix/hmat.h"

class BigWhamIOGen {
private:
  std::string kernel_name_;
  il::int_t spatial_dimension_;
  il::int_t dof_dimension_;
  il::int_t flux_dimension_;

  il::int_t max_leaf_size_;
  double eta_;
  double epsilon_aca_;

  bool is_built_=false;
  bool is_pattern_built_=false;
  double h_representation_time_;
  double hmat_time_;

  std::shared_ptr<bigwham::Hmat<double>> hmat_;
  std::shared_ptr<bigwham::HRepresentation> hr_;
  std::shared_ptr<bigwham::Mesh> mesh_src_;
  std::shared_ptr<bigwham::Mesh> mesh_rec_;
  std::shared_ptr<bigwham::Mesh> mesh_; // for square matrix
  std::shared_ptr<bigwham::BieKernel<double>> ker_obj_; // BieKernel description for the matrix operator
  std::shared_ptr<bigwham::BieKernel<double>> ker_obs_u_; // BieKernel description for computing observation of potential/displacement (u) at points
  std::shared_ptr<bigwham::BieKernel<double>> ker_obs_q_; // BieKernel description for computing observation of 'flux'/'stress' (q) at points

public:
    BigWhamIOGen() {};
    // square matrices
    BigWhamIOGen(const std::vector<double> &coor, const std::vector<int> &conn,
               const std::string &kernel, const std::vector<double> &properties) ;

    // rectangular Hmat
    BigWhamIOGen(const std::vector<double> &coor_src,
                 const std::vector<int> &conn_src,
                 const std::vector<double> &coor_rec,
                 const std::vector<int> &conn_rec, const std::string &kernel,
                 const std::vector<double> &properties);
  ~BigWhamIOGen() {};

  void BuildPattern(const int max_leaf_size, const double eta);

  void BuildHierarchicalMatrix(const int max_leaf_size, const double eta, const double eps_aca); // construct Hierarchical matrix

  // load from file
  void LoadFromFile(const std::string &filename) {
    this->hmat_ = std::make_shared<bigwham::Hmat<double>>(filename);
  }

  void WriteHmatrix(const std::string &filename) {
    this->hmat_->writeToFile(filename);
  }

  il::Array<double> m_yout_; // output vector of matvec
  void HmatDestructor();
  [[nodiscard]] std::vector<double> GetCollocationPoints() const;
  [[nodiscard]] std::vector<long> GetPermutation() const;
  [[nodiscard]] std::vector<long> GetHPattern() const;
  void GetFullBlocks(std::vector<double> &val_list,
                     std::vector<int> &pos_list) const;
  void GetFullBlocks(il::Array<double> &val_list,
                     il::Array<int> &pos_list) const;
  void GetDiagonal(std::vector<double> &val_list) const;
  [[nodiscard]] std::vector<double> MatVec(const std::vector<double> &x) const;
  [[nodiscard]] il::Array<double> MatVec(il::ArrayView<double> x) const;
  [[nodiscard]] std::vector<double> MatVecPerm(const std::vector<double> &x) const;
  [[nodiscard]] il::Array<double> MatVecPerm(il::ArrayView<double> x) const;
  [[nodiscard]] std::vector<double> ConvertToGlobal(const std::vector<double> &x_local) const;
  [[nodiscard]] std::vector<double> ConvertToLocal(const std::vector<double> &x_global) const;
  [[nodiscard]] il::Array<double> ConvertToGlobal(il::ArrayView<double> x_local) const;
  [[nodiscard]] il::Array<double> ConvertToLocal(il::ArrayView<double> x_global) const;
  [[nodiscard]] il::Array<double> ComputePotentials(const std::vector<double> &coor_obs , const il::ArrayView<double> sol_local) const;
  [[nodiscard]] il::Array<double> ComputeFluxes(const std::vector<double> &coor_obs , const il::ArrayView<double> sol_local) const;
  [[nodiscard]] il::Array<double> ComputeDisplacements(const std::vector<double> &coor_obs , const il::ArrayView<double> sol_local) const
  {return ComputePotentials(coor_obs, sol_local);};
  [[nodiscard]] il::Array<double> ComputeStresses(const std::vector<double> &coor_obs , const il::ArrayView<double> sol_local) const
  {return ComputeFluxes(coor_obs,sol_local); };

  void MatVec(il::ArrayView<double> x);

  int MatrixSize(const int k) { return hmat_->size(k); };
  [[nodiscard]] double GetCompressionRatio() const {
    IL_EXPECT_FAST(is_built_);
    return hmat_->compressionRatio();
  };
  [[nodiscard]] double hmat_time() const { return hmat_time_; };
  [[nodiscard]] double pattern_time() const { return h_representation_time_; };
  [[nodiscard]] std::string kernel_name() const { return kernel_name_; };
  [[nodiscard]] int spatial_dimension() const { return spatial_dimension_; };
  [[nodiscard]] int dof_dimension() const { return dof_dimension_; };
  [[nodiscard]] bool is_built() const { return is_built_; };
  void HmatrixDestructor() {
    // this function will free the memory and set the hmat obj to its initial
    // status prior to initialization this will avoid ownership specifications
    // at binding time
    this->hmat_->hmatMemFree();
  }

  int GetOmpThreads() {
    int threads = 1;
#ifdef IL_OPENMP
#pragma omp parallel
    {
#pragma omp single
       std::cout << "NUM OF OMP THREADS: " << omp_get_num_threads() <<
       std::endl;
      threads = omp_get_num_threads();
    }
#endif
    return threads;
  }
};

#endif // BIGWHAM_IO_GEN_H
