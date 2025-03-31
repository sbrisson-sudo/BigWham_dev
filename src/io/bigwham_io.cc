//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 15.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Dec. 2023 - new interface improvements

#include "bigwham_io.h"
#include "bigwham_io_helper.h"

#include "core/be_mesh.h"
#include "core/bie_kernel.h"
#include "core/elastic_properties.h"

#include "elements/point.h"
#include "elements/rectangle.h"
#include "elements/segment.h"
#include "elements/triangle.h"

#include "elasticity/bie_elastostatic.h"
#include "elasticity/fullspace_iso_axisymmetry_flat_unidirectional/bie_elastostatic_axi3d_uni.h"
#include "elasticity/fullspace_iso_sp3d_segment/bie_elastostatic_sp3d.h"
#include "elasticity/fullspace_iso_2d_segment/bie_elastostatic_segment_0_influence.h"
#include "elasticity/fullspace_iso_2d_segment/bie_elastostatic_segment_1_influence.h"
#include "elasticity/fullspace_iso_3d_triangle/bie_elastostatic_triangle_0_influence.h"
#include "elasticity/fullspace_iso_3d_rectangle/bie_elastostatic_rectangle_0_influence.h"
#include "elasticity/fullspace_iso_3d_rectangle/bie_elastostatic_rectangle_0_mode1_influence.h"

#include "elasticity/fullspace_iso_3d_triangle/bie_elastostatic_triangle_2_influence.h"

/* -------------------------------------------------------------------------- */
using namespace bigwham;

// SQUARE matrix case
BigWhamIO::BigWhamIO(const std::vector<double> &coor,
                     const std::vector<int> &conn,
                     const std::string &kernel,
                     const std::vector<double> &properties,
                     const int n_openMP_threads,
                     const bool verbose,
                     const bool homogeneous_size_pattern,
                     const bool useCuda,
                     const int fixed_rank
                    )
{
    this->verbose_ = verbose;
    this->homogeneous_size_pattern_ = homogeneous_size_pattern;
    this->kernel_name_ = kernel;
    this->fixed_rank_ = fixed_rank;

    // Ensuring CUDA compatible arguments
    this->use_cuda_hmat = useCuda;

    if( this->use_cuda_hmat  && (!GetCudaAvailable())){
        std::cerr << "WARNING: Bigwham hasn't been compiled with CUDA support, falling back to CPU matvec" << std::endl;
        this->use_cuda_hmat = false;
    }
    
    if (this->use_cuda_hmat  && (fixed_rank <= 0)) {
        std::cerr << "WARNING: CUDA matvec can't be used without using fixed rank, falling back to CPU matvec" << std::endl;
        this->use_cuda_hmat = false;
    }
    if (this->use_cuda_hmat  && !homogeneous_size_pattern) {
        std::cerr << "WARNING: CUDA matvec can't be used without using homogeneous sized pattern, falling back to CPU matvec" << std::endl;
        this->use_cuda_hmat = false;
    }
    
    // This should be cleaned properly
    this->n_openMP_threads_ = n_openMP_threads;
    int n_available = this->GetAvailableOmpThreads();
    if (n_openMP_threads > n_available)
    {
        this->n_openMP_threads_ = n_available;
    };

#ifdef BIGWHAM_OPENMP
    if (this->verbose_)
            std::cout << "Forcing the number of OpenMP threads to " << this->n_openMP_threads_ << std::endl;
    omp_set_max_active_levels(1);  // Limit parallel region depth
    // omp_set_nested(0);  // Disable nested parallelism

    omp_set_num_threads(this->n_openMP_threads_);
#endif

    if (this->verbose_){
        std::cout << "BigWham using " << this->n_openMP_threads_ << " OpenMP threads\n";
        std::cout << " Now setting things for kernel ... " << kernel_name_
                  << " with properties size " << properties.size() << "\n";
    }
    

    switch (hash_djb2a(kernel_name_))
    {
    case "2DS0-H"_sh:
    { // 2D segment piece-wise ct 0 element, H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 2;
        dof_dimension_ = 2;
        flux_dimension_ = 3;
        int nvertices_per_elt_ = 2;
        using EltType = Segment<0>;
        mesh_ = bigwham::CreateMeshFromVect<Segment<0>>(spatial_dimension_, nvertices_per_elt_,
                                                        coor, conn);
        ker_obj_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, Segment<0>, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        using ObsType = Point<2>;
        ker_obs_u_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, ObsType, bigwham::ElasticKernelType::T>>(
            elas, spatial_dimension_);
        ker_obs_q_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, ObsType, bigwham::ElasticKernelType::W>>(
            elas, spatial_dimension_);
        break;
    }
    case "2DS1-H"_sh:
    { // 2D linear segment , H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 2;
        dof_dimension_ = 2;
        flux_dimension_ = 3;
        int nvertices_per_elt_ = 2;
        using EltType = bigwham::Segment<1>;
        mesh_ = bigwham::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                     coor, conn);
        ker_obj_ = std::make_shared<bigwham::BieElastostatic<EltType, EltType, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        using ObsType = Point<2>;
        ker_obs_q_ = std::make_shared<bigwham::BieElastostatic<EltType, ObsType, bigwham::ElasticKernelType::W>>(
            elas, spatial_dimension_);
        // still missing displacement for that kernel + element
        break;
    }
    case "S3DS0-H"_sh:
    { // simplified 3D segment piece-wise constant, H-kernel
        IL_ASSERT(properties.size() == 3);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 2;
        dof_dimension_ = 2;
        flux_dimension_ = 3;
        int nvertices_per_elt_ = 2;
        using EltType = bigwham::Segment<0>;
        mesh_ = bigwham::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                     coor, conn);
        ker_obj_ = std::make_shared<
            bigwham::BieElastostaticSp3d<EltType, EltType, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        il::Array<double> prop{1, properties[2]};
        ker_obj_->set_kernel_properties(prop);
        // still missing displacement & stresses representation for that kernel + element
        break;
    }
    case "Axi3DS0-H"_sh:
    { // flat axisymmetry unidirectional shear+tensile, 2D ring  piece-wise constant (0) element, H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 2;
        dof_dimension_ = 2;
        int nvertices_per_elt_ = 2;
        using EltType = bigwham::Segment<0>;
        mesh_ = bigwham::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                     coor, conn);
        ker_obj_ = std::make_shared<bigwham::BieElastostaticAxi3D<EltType, EltType, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        break;
    }
    case "3DT0-H"_sh:
    { // 3D triangle P0 element, H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 3;
        dof_dimension_ = 3;
        flux_dimension_ = 6; // 6 stress components
        int nvertices_per_elt_ = 3;
        using EltType = bigwham::Triangle<0>;
        mesh_ = bigwham::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                     coor, conn);
        ker_obj_ = std::make_shared<
            bigwham::BieElastostatic<EltType, EltType, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        using ObsType = Point<3>;
        ker_obs_u_ = std::make_shared<bigwham::BieElastostatic<EltType, ObsType, bigwham::ElasticKernelType::T>>(
            elas, spatial_dimension_);
        ker_obs_q_ = std::make_shared<bigwham::BieElastostatic<EltType, ObsType, bigwham::ElasticKernelType::W>>(
            elas, spatial_dimension_);
        break;
    }
    case "3DT6-H"_sh:
    { // 3D triangle P2 element, H-kernel (currently, there is an error in this element (some chges of coordinates or others)
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 3;
        dof_dimension_ = 3;
        flux_dimension_ = 6; // 6 stress components
        int nvertices_per_elt_ = 3;
        using EltType = bigwham::Triangle<2>;
        mesh_ = bigwham::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                     coor, conn);
        ker_obj_ = std::make_shared<
            bigwham::BieElastostatic<EltType, EltType, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        // still missing displacement & stresses representation for that kernel + element
        // using ObsType = Point<3>;
        break;
    }
    case "3DR0-H"_sh:
    { // 3D rectangular P0 element, H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 3;
        dof_dimension_ = 3;
        flux_dimension_ = 6; // 6 stress components
        int nvertices_per_elt_ = 4;
        using EltType = bigwham::Rectangle<0>;
        mesh_ = bigwham::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                     coor, conn);
        ker_obj_ = std::make_shared<
            bigwham::BieElastostatic<EltType, EltType, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        using ObsType = Point<3>;
        ker_obs_u_ = std::make_shared<bigwham::BieElastostatic<EltType, ObsType, bigwham::ElasticKernelType::T>>(
            elas, spatial_dimension_);
        ker_obs_q_ = std::make_shared<bigwham::BieElastostatic<EltType, ObsType, bigwham::ElasticKernelType::W>>(
            elas, spatial_dimension_);
        break;
    }
    case "3DR0-H-mode1"_sh:
    { // 3D rectangular P0 element - tensile only, H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 3;
        dof_dimension_ = 1;
        flux_dimension_ = 6; // 6 stress components
        int nvertices_per_elt_ = 4;
        using EltType = bigwham::Rectangle<0>;
        mesh_ = bigwham::CreateMeshFromVect<EltType>(spatial_dimension_, nvertices_per_elt_,
                                                     coor, conn);
        ker_obj_ = std::make_shared<
            bigwham::BieElastostaticModeI<EltType, EltType, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        // observation kernels to implement....
        // still missing displacement & stresses representation for that kernel + element
        break;
    }
    default:
    {
        std::cout << "wrong inputs -abort \n";
        il::abort();
    }
    }
    mesh_src_ = mesh_;
    mesh_rec_ = mesh_;
};

// RECTANGULAR Hmatrix case
BigWhamIO::BigWhamIO(const std::vector<double> &coor_src,
                     const std::vector<int> &conn_src,
                     const std::vector<double> &coor_rec,
                     const std::vector<int> &conn_rec, const std::string &kernel,
                     const std::vector<double> &properties,
                     const int n_openMP_threads,
                    const bool verbose,
                    const bool homogeneous_size_pattern,
                    const bool useCuda,
                    const int fixed_rank)
{
    this->verbose_ = verbose;
    this->homogeneous_size_pattern_ = homogeneous_size_pattern;
    this->fixed_rank_ = fixed_rank;

#ifdef USE_CUDA
    this->use_cuda_hmat = useCuda;
#endif

    this->n_openMP_threads_ = n_openMP_threads;
    auto n_available = this->GetAvailableOmpThreads();
    if (n_openMP_threads > n_available)
    {
        this->n_openMP_threads_ = n_available;
    };

#ifdef BIGWHAM_OPENMP
    if (this->verbose_)
        std::cout << "Forcing the number of OpenMP threads to " << this->n_openMP_threads_ << std::endl;
    omp_set_num_threads(this->n_openMP_threads_);
    // omp_set_nested(0);  // Disable nested parallelism
    omp_set_max_active_levels(1);  // Limit parallel region depth
#endif

    this->kernel_name_ = kernel;

    if (this->verbose_){
        std::cout << "BigWham using " << this->n_openMP_threads_ << " OpenMP threads\n";
        std::cout << "Now setting things for kernel ... " << kernel_name_
                  << " with properties size " << properties.size() << "\n";
    }

    

    switch (hash_djb2a(kernel_name_))
    {
    case "2DS0-2DS0-H"_sh:
    { // 2D segment piece-wise ct 0 element, H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 2;
        dof_dimension_ = 2;
        flux_dimension_ = 3;
        using src_elem = Segment<0>;
        using rec_elem = Segment<0>;
        mesh_src_ = bigwham::CreateMeshFromVect<src_elem>(
            spatial_dimension_, /* num vertices */ 2, coor_src, conn_src);
        mesh_rec_ = bigwham::CreateMeshFromVect<rec_elem>(
            spatial_dimension_, /* num vertices */ 2, coor_rec, conn_rec);
        ker_obj_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, Segment<0>, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        using ObsType = Point<2>;
        ker_obs_u_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, ObsType, bigwham::ElasticKernelType::T>>(
            elas, spatial_dimension_);
        ker_obs_q_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, ObsType, bigwham::ElasticKernelType::W>>(
            elas, spatial_dimension_);
        break;
    }
    case "2DS1-2DS1-H"_sh:
    { // 2D segment piece-wise linear element, H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 2;
        dof_dimension_ = 2;
        flux_dimension_ = 3;
        using src_elem = Segment<1>;
        using rec_elem = Segment<1>;
        mesh_src_ = bigwham::CreateMeshFromVect<src_elem>(
            spatial_dimension_, /* num vertices */ 2, coor_src, conn_src);
        mesh_rec_ = bigwham::CreateMeshFromVect<rec_elem>(
            spatial_dimension_, /* num vertices */ 2, coor_rec, conn_rec);
        ker_obj_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, Segment<0>, bigwham::ElasticKernelType::H>>(
            elas, spatial_dimension_);
        using ObsType = Point<2>;
        ker_obs_u_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, ObsType, bigwham::ElasticKernelType::T>>(
            elas, spatial_dimension_);
        ker_obs_q_ = std::make_shared<bigwham::BieElastostatic<Segment<0>, ObsType, bigwham::ElasticKernelType::W>>(
            elas, spatial_dimension_);
        break;
    }
    case "2DS0-2DP-T"_sh:
    {
        // 2D Segment0 and 2D Point, computation of T kernel (displacement due to dislocation)
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        dof_dimension_ = 2;
        spatial_dimension_ = 2;
        flux_dimension_ = 3;
        using src_elem = Segment<0>;
        using rec_elem = Point<2>;
        mesh_src_ = bigwham::CreateMeshFromVect<src_elem>(
            spatial_dimension_, /* num vertices */ 2, coor_src, conn_src);
        mesh_rec_ = bigwham::CreateMeshFromVect<rec_elem>(
            spatial_dimension_, /* num vertices */ 1, coor_rec, conn_rec);
        ker_obj_ = std::make_shared<
            BieElastostatic<src_elem, rec_elem, ElasticKernelType::T>>(
            elas, spatial_dimension_);
        break;
    }
    case "3DT0-3DT0-H"_sh:
    { // 3D triangle P0 elements
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        dof_dimension_ = 3;
        spatial_dimension_ = 3;
        flux_dimension_ = 6;
        using src_elem = Triangle<0>;
        using rec_elem = Triangle<0>;
        mesh_src_ = bigwham::CreateMeshFromVect<src_elem>(
            spatial_dimension_, /* num vertices */ 3, coor_src, conn_src);
        mesh_rec_ = bigwham::CreateMeshFromVect<rec_elem>(
            spatial_dimension_, /* num vertices */ 3, coor_rec, conn_rec);
        ker_obj_ = std::make_shared<BieElastostatic<src_elem, rec_elem, ElasticKernelType::H>>(
            elas, spatial_dimension_);
        // observations....
        using ObsType = Point<3>;
        ker_obs_u_ = std::make_shared<bigwham::BieElastostatic<src_elem, ObsType, bigwham::ElasticKernelType::T>>(
            elas, spatial_dimension_);
        ker_obs_q_ = std::make_shared<bigwham::BieElastostatic<src_elem, ObsType, bigwham::ElasticKernelType::W>>(
            elas, spatial_dimension_);
        break;
    }
    case "3DR0-3DR0-H"_sh:
    { // 3D Rectangular P0 elements
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        dof_dimension_ = 3;
        spatial_dimension_ = 3;
        flux_dimension_ = 6;
        using src_elem = Rectangle<0>;
        using rec_elem = Rectangle<0>;
        mesh_src_ = bigwham::CreateMeshFromVect<src_elem>(
            spatial_dimension_, /* num vertices */ 4, coor_src, conn_src);
        mesh_rec_ = bigwham::CreateMeshFromVect<rec_elem>(
            spatial_dimension_, /* num vertices */ 4, coor_rec, conn_rec);
        ker_obj_ = std::make_shared<BieElastostatic<src_elem, rec_elem, ElasticKernelType::H>>(
            elas, spatial_dimension_);
        // observations....
        using ObsType = Point<3>;
        ker_obs_u_ = std::make_shared<bigwham::BieElastostatic<src_elem, ObsType, bigwham::ElasticKernelType::T>>(
            elas, spatial_dimension_);
        ker_obs_q_ = std::make_shared<bigwham::BieElastostatic<src_elem, ObsType, bigwham::ElasticKernelType::W>>(
            elas, spatial_dimension_);
        break;
    }
    case "3DR0-3DR0-H-mode1"_sh:
    { // 3D rectangular P0 element - tensile only, H-kernel
        IL_ASSERT(properties.size() == 2);
        ElasticProperties elas(properties[0], properties[1]);
        spatial_dimension_ = 3;
        dof_dimension_ = 1;
        flux_dimension_ = 6; // 6 stress components
        using src_elem = Rectangle<0>;
        using rec_elem = Rectangle<0>;
        mesh_src_ = bigwham::CreateMeshFromVect<src_elem>(
            spatial_dimension_, /* num vertices */ 4, coor_src, conn_src);
        mesh_rec_ = bigwham::CreateMeshFromVect<rec_elem>(
            spatial_dimension_, /* num vertices */ 4, coor_rec, conn_rec);
        ker_obj_ = std::make_shared<BieElastostatic<src_elem, rec_elem, ElasticKernelType::H>>(
            elas, spatial_dimension_);
        // observation kernels to implement....
        // still missing displacement & stresses representation for that kernel + element
        break;
    }
    default:
    {
        std::cout << "wrong inputs -abort \n";
        il::abort();
    }
    }
    mesh_ = mesh_src_;
}

/* -------------------------------------------------------------------------- */

// Method for pattern Construction
// API to the Hierarchical representation function
void BigWhamIO::BuildPattern(const int max_leaf_size, const double eta)
{

    max_leaf_size_ = max_leaf_size;
    eta_ = eta;
    il::Timer tt;
    if (this->verbose_){
        std::cout << "--------------------\n";
        std::cout << "Hierarchical representation creation ...\n";
        if (this->homogeneous_size_pattern_){
            std::cout << "Using homegeneous size clustering\n";
        }
    }
    tt.Start();
    this->hr_ = HRepresentationRectangularMatrix(mesh_src_, mesh_rec_, max_leaf_size_, eta_, verbose_, homogeneous_size_pattern_);
    is_pattern_built_ = true;
    tt.Stop();
    if (this->verbose_)
        std::cout << "Hierarchical representation complete.\n";
}
/* -------------------------------------------------------------------------- */

// method to only get the permutation of the source mesh...

std::vector<long> BigWhamIO::GetPermutation() const
{
    // this function can be call either before or after the Hmat is built.
    std::vector<long> permut;
    if (is_built_ || is_pattern_built_)
    {
        permut.assign(hr_->permutation_1_.size(), 0);
        for (il::int_t i = 0; i < hr_->permutation_1_.size(); i++)
        {
            permut[i] = hr_->permutation_1_[i];
        }
    }
    else
    { // the H-mat is not built, here the use case is to just get permutation
        il::Array2D<double> Xcol_source = mesh_src_->collocation_points();
        Cluster cluster_s = cluster(max_leaf_size_, il::io, Xcol_source);
        il::Array<il::int_t> per_i = cluster_s.permutation;
        // convert to std::vect
        permut.assign(per_i.size(), 0);
        for (il::int_t i = 0; i < per_i.size(); i++)
        {
            permut[i] = per_i[i];
        }
    }
    return permut;
}

std::vector<long> BigWhamIO::GetPermutationReceivers() const
{
    // this function can be call either before or after the Hmat is built.
    std::vector<long> permut;
    if (is_built_ || is_pattern_built_)
    {
        permut.assign(hr_->permutation_0_.size(), 0);
        for (il::int_t i = 0; i < hr_->permutation_0_.size(); i++)
        {
            permut[i] = hr_->permutation_0_[i];
        }
    }
    else
    { // the H-mat is not built, here the use case is to just get permutation
        std::cerr << "Returning sources permutation" << std::endl;
        il::Array2D<double> Xcol_source = mesh_src_->collocation_points();
        Cluster cluster_s = cluster(max_leaf_size_, il::io, Xcol_source);
        il::Array<il::int_t> per_i = cluster_s.permutation;
        // convert to std::vect
        permut.assign(per_i.size(), 0);
        for (il::int_t i = 0; i < per_i.size(); i++)
        {
            permut[i] = per_i[i];
        }
    }
    return permut;
}
/* -------------------------------------------------------------------------- */

// construct Hierarchical matrix
void BigWhamIO::BuildHierarchicalMatrix(const int max_leaf_size, const double eta, const double eps_aca)
{

    epsilon_aca_ = eps_aca;
    if (!is_pattern_built_)
    {
        BigWhamIO::BuildPattern(max_leaf_size, eta);
    }

    il::Timer tt;
    if (this->verbose_){
        std::cout << "--------------------\n";
        std::cout << "Populating Hierarchical matrix ...\n";
    }
    
    tt.Start();
    bigwham::BieMatrixGenerator<double> M(mesh_src_, mesh_rec_, ker_obj_, hr_);

#ifdef USE_CUDA
    if (this->use_cuda_hmat){
        hmat_ = std::make_shared<HmatCuda<double>>(M, epsilon_aca_, this->n_openMP_threads_, this->verbose_, this->fixed_rank_);
    } else {
#endif

    hmat_ = std::make_shared<Hmat<double>>(M, epsilon_aca_, this->n_openMP_threads_, this->verbose_, this->fixed_rank_);

#ifdef USE_CUDA
    }
#endif

    tt.Stop();
    hmat_time_ = tt.time();
    if (hmat_->isBuilt())
    {
        is_built_ = true;
        IL_EXPECT_FAST(dof_dimension_ == hmat_->dofDimension());
        if (this->verbose_)
            std::cout << "Hierarchical matrix construction complete.\n";
        double test_cr = hmat_->compressionRatio();
        if (this->verbose_){
            std::cout << "Compression Ratio = " << test_cr << ", eps_aca = " << epsilon_aca_
                  << ", eta = " << eta_ << "\n";
            std::cout << "Hierarchical matrix  construction time = :  " << hmat_time_ << "\n";
        }
        
    }
    else
    {
        is_built_ = false;
    }
    if (this->verbose_){
        std::cout << "BigWham constructed Hmat of size "
                << "(" << mesh_rec_->num_collocation_points()
                << " x " << mesh_rec_->spatial_dimension() << ")"
                << " X " << "(" << mesh_src_->num_collocation_points()
                << " x " << mesh_rec_->spatial_dimension() << ")" << "\n";
        std::cout << "--------------------\n";
    }
}
/* -------------------------------------------------------------------------- */

std::vector<double> BigWhamIO::GetHPattern() const
{
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

    bigwham::HPattern pattern = hr_->pattern_;
    if (is_built_)
    {
        pattern = hmat_->pattern();
    }

    long numberofblocks = pattern.n_B;
    long len = 6 * numberofblocks;
    std::cout << "number of blocks " << numberofblocks << "\n";

    // std::vector<long> patternlist(len, 0);
    std::vector<double> patternlist(len, 0); // to contain error

    int index = 0;
    //  starts with full rank
    for (il::int_t j = 0; j < pattern.n_FRB; j++)
    {
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
    for (il::int_t j = 0; j < pattern.n_LRB; j++)
    {
        patternlist[index++] = pattern.LRB_pattern(1, j);
        patternlist[index++] = pattern.LRB_pattern(2, j);
        patternlist[index++] = pattern.LRB_pattern(3, j);
        patternlist[index++] = pattern.LRB_pattern(4, j);
        patternlist[index++] = 1;

        // If not using fixde rank we return the rank, otherwise we return the
        // frobenius norm error
        if (this->fixed_rank_ > 0){
            patternlist[index++] = hmat_->getLowRankBlock(j)->error_on_approximation;
        } else {
            patternlist[index++] = pattern.LRB_pattern(5, j); // the rank
        }
        
    }
    // return a row major flatten vector
    return patternlist;
}

double BigWhamIO::GetMaxErrorACA() const {
    bigwham::HPattern pattern = hmat_->pattern();
    double max_error = 0;    
    for (il::int_t j = 0; j < pattern.n_LRB; j++){
        double error_aca = hmat_->getLowRankBlock(j)->error_on_approximation;
        if (error_aca > max_error) max_error = error_aca;
    }
    return max_error;
}

/* --------------------------------------------------------------------------*/

void BigWhamIO::GetFullBlocks(il::Array<double> &val_list,
                              il::Array<int> &pos_list) const
{
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
    for (il::int_t i = 0; i < values.size(); i++)
    {
        val_list[i] = values[i];
    }
    pos_list.Resize(positions.size());
    for (il::int_t i = 0; i < positions.size(); i++)
    {
        pos_list[i] = positions[i];
    }
    // std::cout << "number of entries " << val_list.size() << " - "
    //           << pos_list.size() << "\n";
    // std::cout << " End of Bigwhamio_old getFullBlocks \n";
}
/* --------------------------------------------------------------------------*/

void BigWhamIO::GetFullBlocks(std::vector<double> &val_list,
                              std::vector<int> &pos_list) const
{
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
    for (il::int_t i = 0; i < values.size(); i++)
    {
        val_list.push_back(values[i]);
    }
    pos_list.reserve(positions.size());
    for (il::int_t i = 0; i < positions.size(); i++)
    {
        pos_list.push_back(positions[i]);
    }
    // std::cout << "number of entries " << val_list.size() << " - "
    //           << pos_list.size() << "\n";
    // std::cout << " End of Bigwhamio_old getFullBlocks \n";
}

void BigWhamIO::GetFullBlocksPerm(il::Array<double> &val_list,
    il::Array<int> &pos_list) const
{
// return the full dense block entries of the hmat as
// flattened lists
// val_list(i) = H(pos_list(2*i),pos_list(2*i+1));
// output in the original dof state (accounting for the permutation)

IL_EXPECT_FAST(is_built_);

il::Array<double> values{};
il::Array<int> positions{};
hmat_->fullBlocksPerm(il::io, values, positions);

//    std::cout << " checking values size" << values.size() <<  "\n";
val_list.Resize(values.size());
for (il::int_t i = 0; i < values.size(); i++)
{
val_list[i] = values[i];
}
pos_list.Resize(positions.size());
for (il::int_t i = 0; i < positions.size(); i++)
{
pos_list[i] = positions[i];
}
// std::cout << "number of entries " << val_list.size() << " - "
//           << pos_list.size() << "\n";
// std::cout << " End of Bigwhamio_old getFullBlocks \n";
}
/* --------------------------------------------------------------------------*/

void BigWhamIO::GetFullBlocksPerm(std::vector<double> &val_list,
    std::vector<int> &pos_list) const
{
// return the full dense block entries of the hmat as
// flattened lists
// val_list(i) = H(pos_list(2*i),pos_list(2*i+1));
// output in the original dof state (accounting for the permutation)

IL_EXPECT_FAST(is_built_);

il::Array<double> values{};
il::Array<int> positions{};
hmat_->fullBlocksPerm(il::io, values, positions);
//    std::cout << " checking values size" << values.size() <<  "\n";
val_list.reserve(values.size());
for (il::int_t i = 0; i < values.size(); i++)
{
val_list.push_back(values[i]);
}
pos_list.reserve(positions.size());
for (il::int_t i = 0; i < positions.size(); i++)
{
pos_list.push_back(positions[i]);
}
// std::cout << "number of entries " << val_list.size() << " - "
//           << pos_list.size() << "\n";
// std::cout << " End of Bigwhamio_old getFullBlocks \n";
}
/* --------------------------------------------------------------------------*/

void BigWhamIO::GetDiagonal(std::vector<double> &val_list) const
{
    // return the diagonal of the h-matrix
    // output in the original dof state (accounting for the permutation)

    IL_EXPECT_FAST(is_built_);
    val_list = hmat_->diagonalOriginal();
}
/* -------------------------------------------------------------------------- */

std::vector<double> BigWhamIO::MatVec(const std::vector<double> &x) const
{
    // in the original / natural ordering
    IL_EXPECT_FAST(this->is_built_);
    std::vector<double> y = hmat_->matvecOriginal(x);
    return y;
}
/* -------------------------------------------------------------------------- */
il::Array<double> BigWhamIO::MatVec(il::ArrayView<double> x) const
{
    // in the original / natural ordering
    IL_EXPECT_FAST(this->is_built_);
    auto y = hmat_->matvecOriginal(x);
    return y;
}
/* -------------------------------------------------------------------------- */
void BigWhamIO::MatVecVoid(const il::ArrayView<double> xin)
{
    // in the original / natural ordering
    // this function is used in python/julia to avoid copying the output
    // but use internal variable m_yout_ to store the result
    // it will save creating new memory at each call
    // this->m_yout_ = hmat_* xin
    IL_EXPECT_FAST(this->is_built_);
    hmat_->matvecOriginal(xin, this->m_yout_);
}
/* -------------------------------------------------------------------------- */
void BigWhamIO::MatVecVoid(const il::ArrayView<double> xin, il::ArrayEdit<double> yout)
{
    // in the original / natural ordering
    // this function is used in python/julia to avoid copying the output
    // but use internal variable m_yout_ to store the result
    // it will save creating new memory at each call
    // yout = hmat_* xin
    IL_EXPECT_FAST(this->is_built_);
    hmat_->matvecOriginal(xin, yout);
}
/* -------------------------------------------------------------------------- */
il::Array<double> BigWhamIO::MatVecPerm(il::ArrayView<double> x) const
{
    // in the permutted state.
    IL_EXPECT_FAST(this->is_built_);
    auto y = hmat_->matvec(x);
    return y;
}
/* -------------------------------------------------------------------------- */
std::vector<double> BigWhamIO::MatVecPerm(const std::vector<double> &x) const
{
    // in the permutted state.
    IL_EXPECT_FAST(this->is_built_);
    IL_EXPECT_FAST(hmat_->size(1) == x.size());
    std::vector<double> y = hmat_->matvec(x);
    return y;
}

/* -------------------------------------------------------------------------- */
il::Array<double> BigWhamIO::ConvertToGlobal(il::ArrayView<double> x_local) const
{
    // Input: x in original state (not permutted)
    // Output: in original state (not permutted)
    return mesh_rec_->ConvertToGlobal(x_local);
}
/* -------------------------------------------------------------------------- */
il::Array<double> BigWhamIO::ConvertToLocal(il::ArrayView<double> x_global) const
{
    // Input: x in original state (not permutted)
    // Output: in original state (not permutted)
    return mesh_src_->ConvertToLocal(x_global);
}
/* -------------------------------------------------------------------------- */
std::vector<double> BigWhamIO::ConvertToGlobal(const std::vector<double> &x_local) const
{
    // Input: x in original state (not permutted)
    // Output: in original state (not permutted)
    return mesh_rec_->ConvertToGlobal(x_local);
}
/* -------------------------------------------------------------------------- */

std::vector<double> BigWhamIO::ConvertToLocal(const std::vector<double> &x_global) const
{
    // Input: x in original state (not permutted)
    // Output: in original state (not permutted)
    return mesh_src_->ConvertToLocal(x_global);
}
/* -------------------------------------------------------------------------- */
std::vector<double> BigWhamIO::GetElementNormals() const
{
    /*
  # normals
  Get the normal vectors of the elements
  :return: 1D array of the normal vectors of the elements
  3d case:
  [n0x,n0y,n0z,n1x,n1y,n1z,...]
  2d case:
  [n0x,n0y,n1x,n1y,...]
    */

    auto mesh = this->mesh_src_;
    auto dim = mesh->spatial_dimension();
    std::vector<double> normals;
    normals.assign(mesh->num_elements() * dim, 0.0);
    // normals.Resize(dim * mesh->num_elements(), 0.0);
#pragma omp parallel num_threads(this->n_openMP_threads_)
    {
#pragma omp parallel for
        for (il::int_t i = 0; i < mesh->num_elements(); ++i)
        {
            auto elem = mesh->GetElement(i);
            auto n = elem->normal();
            for (il::int_t j = 0; j < dim; ++j)
            {
                normals[i * dim + j] = n[j];
            }
        }
    }
    return normals;
}

/* -------------------------------------------------------------------------- */
std::vector<double> BigWhamIO::ComputeTractionsFromUniformStress(const std::vector<double> &stress) const
{
    // stress = s_11,s22,s_12 in 2D , s_11,s_22,s_33,s_12,s_13,s_23 in 3D
    // return traction on a element facet in the global system

    auto mesh = this->mesh_src_;
    auto dim = mesh->spatial_dimension();
    if (dim == 2)
    {
        IL_EXPECT_FAST(stress.size() == 3);
    }
    if (dim == 3)
    {
        IL_EXPECT_FAST(stress.size() == 6);
    }

    std::vector<double> global_tractions;
    il::Array2D<double> sij{dim, dim, 0.};
    int k = 0;
    for (int i = 0; i < dim; i++)
    { // diag
        sij(i, i) = stress[k];
        k++;
    }
    for (int i = 1; i < dim; i++)
    { // off diag
        sij(0, i) = stress[k];
        sij(i, 0) = stress[k];
        k++;
    }
    if (dim == 3)
    { // s_23
        sij(1, 2) = stress[k];
        sij(2, 1) = stress[k];
    }

    global_tractions.assign(mesh->num_elements() * dim, 0.0);
#pragma omp parallel num_threads(this->n_openMP_threads_)
    {
#pragma omp parallel for
        for (il::int_t i = 0; i < mesh->num_elements(); ++i)
        {
            auto elem = mesh->GetElement(i);
            auto n = elem->normal();
            auto g_t = il::dot(sij, n);
            for (il::int_t j = 0; j < dim; ++j)
            {
                global_tractions[i * dim + j] = g_t[j];
            }
        }
    }
    return global_tractions;
}
/* -------------------------------------------------------------------------- */
std::vector<double> BigWhamIO::GetRotationMatrix() const
{
    /*
  # rotation_matrix
  Get the rotation matrix of the elements
  :return: 1D array of the rotation matrix of the elements
  3d case:
  [r00,r01,r02,r10,r11,r12,r20,r21,r22,...]
  2d case:
  [r00,r01,r10,r11,...]
     */

    auto mesh = this->mesh_src_;
    auto dim = mesh->spatial_dimension();
    std::vector<double> flattend_rotation_matrix;
    int num_elements = mesh->num_elements();
    int num_total_collocation_points = mesh->num_collocation_points();
    flattend_rotation_matrix.assign(dim * dim * num_total_collocation_points, 0.0);
    //  flattend_rotation_matrix.Resize(dim * dim * num_total_collocation_points, 0.0);
#pragma omp parallel num_threads(this->n_openMP_threads_)
    {
#pragma omp parallel for
        for (il::int_t i = 0; i < mesh->num_elements(); ++i)
        {
            auto elem = mesh->GetElement(i);
            auto mat = elem->rotation_matrix();
            int num_local_collocation_points = elem->num_collocation_points();
            for (il::int_t col_id = 0; col_id < num_local_collocation_points; ++col_id)
            {
                for (il::int_t j = 0; j < dim; ++j)
                {
                    for (il::int_t k = 0; k < dim; ++k)
                    {
                        flattend_rotation_matrix[i * dim * dim * num_local_collocation_points + col_id * dim * dim + j * dim + k] = mat(j, k);
                        // flattend_rotation_matrix.push_back(mat(j, k));
                    }
                }
            }
        }
    }
    return flattend_rotation_matrix;
}

/* -------------------------------------------------------------------------- */
std::vector<double> BigWhamIO::GetCollocationPoints() const
{
    //  IL_EXPECT_FAST(is_built_);

    auto col_pts = mesh_rec_->collocation_points(); // this should be the receiver mesh for  generality
    IL_EXPECT_FAST(col_pts.size(1) == spatial_dimension_);
    il::int_t npoints = col_pts.size(0);
    std::vector<double> flat_col;
    flat_col.assign(npoints * spatial_dimension_, 0.);
    int index = 0;
    for (il::int_t i = 0; i < col_pts.size(0); i++)
    {
        for (il::int_t j = 0; j < col_pts.size(1); j++)
        {
            flat_col[index] = col_pts(i, j);
            index++;
        }
    }
    return flat_col;
}
/* -------------------------------------------------------------------------- */

il::Array<double> BigWhamIO::ComputePotentials(const std::vector<double> &coor_obs, const il::ArrayView<double> sol_local) const
{ // return the potential in the global system of coordinates, the solution on the source mesh being defined in the local eletment system of the source mesh
    IL_EXPECT_FAST(sol_local.size() == mesh_src_->num_collocation_points() * dof_dimension_);
    IL_EXPECT_FAST(coor_obs.size() % spatial_dimension_ == 0);
    if (this->ker_obs_u_ == nullptr)
    {
        std::cout << "Kernel Not Implemented !\n";
    }
    IL_EXPECT_FAST(this->ker_obs_u_ != nullptr);

    il::int_t npts = coor_obs.size() / spatial_dimension_;
    std::shared_ptr<bigwham::Mesh> mesh_obs;
    std::vector<int> conn_obs(npts, 0);
    for (int i = 0; i < npts; i++)
    {
        // dummy connectivity for a mesh of points.
        conn_obs[i] = i;
    }
    IL_EXPECT_FAST(spatial_dimension_ == 2 || spatial_dimension_ == 3);

    switch (spatial_dimension_)
    {
    case 2:
    {
        mesh_obs = bigwham::CreateMeshFromVect<bigwham::Point<2>>(spatial_dimension_, 1, coor_obs, conn_obs);
        break;
    }
    case 3:
    {
        mesh_obs = bigwham::CreateMeshFromVect<bigwham::Point<3>>(spatial_dimension_, 1, coor_obs, conn_obs);
        break;
    }
    default:
    {
        std::cout << "Invalid spatial_dimension ! \n";
        break;
    }
    }

    il::Array<double> obs_potential{npts * dof_dimension_, 0.};

// loop on source collocation points
#pragma omp parallel if (std::sqrt((mesh_src_->num_collocation_points()) * (mesh_obs->num_collocation_points())) > 400)
    {
#pragma omp for
        for (il::int_t i = 0; i < mesh_src_->num_collocation_points(); i++)
        {

            il::int_t e_i = this->mesh_src_->GetElementId(i);
            il::int_t is_l = this->mesh_src_->GetElementCollocationId(i);
            auto source_element = this->mesh_src_->GetElement(e_i);

            // local solution (assumed in original ordering....)
            il::Array<double> elt_solu{dof_dimension_};
            for (il::int_t k = 0; k < dof_dimension_; k++)
            {
                elt_solu[k] = sol_local[i * dof_dimension_ + k];
            }

            // loop on obs points mesh
            for (il::int_t j_obs = 0; j_obs < mesh_obs->num_collocation_points(); j_obs++)
            {
                il::int_t e_j_r = mesh_obs->GetElementId(j_obs);
                auto receiver_element = mesh_obs->GetElement(e_j_r);
                il::int_t ir_l = mesh_obs->GetElementCollocationId(j_obs);
                std::vector<double> st = this->ker_obs_u_->influence(*source_element, is_l, *receiver_element, ir_l);
                for (il::int_t j = 0; j < dof_dimension_; j++)
                {
                    for (il::int_t k = 0; k < dof_dimension_; k++)
                    {
                        obs_potential[j_obs * dof_dimension_ + j] += st[k * dof_dimension_ + j] * elt_solu[k];
                    }
                }
            }
        }
    }

    return obs_potential;
}

il::Array<double> BigWhamIO::ComputeFluxes(const std::vector<double> &coor_obs, const il::ArrayView<double> sol_local) const
{
    // return the potential in the global system of coordinates, the solution on the source mesh being defined in the local eletment system of the source mesh
    IL_EXPECT_FAST(sol_local.size() == mesh_src_->num_collocation_points() * dof_dimension_);
    IL_EXPECT_FAST(coor_obs.size() % spatial_dimension_ == 0);
    if (this->ker_obs_q_ == nullptr)
    {
        std::cout << "Kernel Not Implemented !\n";
    }
    IL_EXPECT_FAST(this->ker_obs_q_ != nullptr);

    il::int_t npts = coor_obs.size() / spatial_dimension_;

    std::shared_ptr<bigwham::Mesh> mesh_obs;
    std::vector<int> conn_obs(npts, 0);
    for (int i = 0; i < npts; i++)
    {
        // dummy connectivity for a mesh of points.
        conn_obs[i] = i;
    }
    IL_EXPECT_FAST(spatial_dimension_ == 2 || spatial_dimension_ == 3);

    switch (spatial_dimension_)
    {
    case 2:
    {
        mesh_obs = bigwham::CreateMeshFromVect<bigwham::Point<2>>(spatial_dimension_, 1, coor_obs, conn_obs);
        break;
    }
    case 3:
    {
        mesh_obs = bigwham::CreateMeshFromVect<bigwham::Point<3>>(spatial_dimension_, 1, coor_obs, conn_obs);
        break;
    }
    default:
    {
        std::cout << "Invalid spatial_dimension ! \n";
        break;
    }
    }

    il::Array<double> obs_flux{npts * flux_dimension_, 0.};

// loop on source collocation points
#pragma omp parallel if (std::sqrt((mesh_src_->num_collocation_points()) * (mesh_obs->num_collocation_points())) > 400)
    {
#pragma omp for
        for (il::int_t i = 0; i < mesh_src_->num_collocation_points(); i++)
        {

            il::int_t e_i = this->mesh_src_->GetElementId(i);
            il::int_t is_l = this->mesh_src_->GetElementCollocationId(i);

            auto source_element = this->mesh_src_->GetElement(e_i);
            // local solution (assumed in original ordering....)
            il::Array<double> elt_solu{dof_dimension_}, e_potential{dof_dimension_, 0.};
            for (il::int_t k = 0; k < dof_dimension_; k++)
            {
                elt_solu[k] = sol_local[i * dof_dimension_ + k];
            }
            // loop on obs points mesh
            for (il::int_t j_obs = 0; j_obs < mesh_obs->num_collocation_points(); j_obs++)
            {
                il::int_t e_i_r = mesh_obs->GetElementId(j_obs);
                auto receiver_element = mesh_obs->GetElement(e_i_r);
                il::int_t ir_l = mesh_obs->GetElementCollocationId(j_obs);
                std::vector<double> st = this->ker_obs_q_->influence(*source_element, is_l, *receiver_element, ir_l);
                for (il::int_t j = 0; j < flux_dimension_; j++)
                {
                    for (il::int_t k = 0; k < dof_dimension_; k++)
                    {
                        obs_flux[j_obs * flux_dimension_ + j] += st[flux_dimension_ * k + j] * elt_solu[k];
                    }
                }
            }
        }
    }
    return obs_flux;
}


bool BigWhamIO::GetCudaAvailable(){
#ifdef USE_CUDA
    return true;
#endif 
    return false;
}

