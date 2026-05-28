#include <memory>
#include <iostream>
#include <iomanip> // For std::setw


#include "cnpy.h"
#include <gsl/gsl_errno.h>


#include <il/core/core.h>


// Bigwham imports
#include "core/be_mesh.h"
#include "io/bigwham_io_helper.h"

#include "core/elastic_properties.h"
#include "elements/segment.h"

#include "core/bie_kernel.h"
#include "elasticity/bie_elastostatic.h"
#include "elasticity/bie_elastostatic_eigenstrain.h"

#include "hmat/hierarchical_representation.h"
#include "hmat/bie_matrix_generator.h"
#include "hmat/bie_matrix_generator_by_dof.h"

#include <hmat/arrayFunctor/FullMatrix.h>
#include "hmat/compression/adaptiveCrossApproximation.h"



using namespace bigwham;

template <typename T>
void prettyPrintArray2D(const il::Array2D<T> &M) {
    // Get the dimensions of the matrix
    il::int_t rows = M.size(0);
    il::int_t cols = M.size(1);

    // Determine the maximum width needed for any element in the matrix
    int maxWidth = 0;
    for (il::int_t i = 0; i < rows; ++i) {
        for (il::int_t j = 0; j < cols; ++j) {
            std::ostringstream oss;
            oss << M(i, j);
            int width = oss.str().length();
            if (width > maxWidth) {
                maxWidth = width;
            }
        }
    }

    // Print the matrix with proper formatting
    for (il::int_t i = 0; i < rows; ++i) {
        for (il::int_t j = 0; j < cols; ++j) {
            std::cout << std::setw(maxWidth + 2) << M(i, j);
        }
        std::cout << std::endl;
    }
}


int main(int argc, char * argv[]) {

    gsl_set_error_handler_off();

    // get mesh file names form command line arguments
    if (argc < 4){
        std::cerr << "Usage : ./mat_vec <coor.npy> <conn_tri.npy> <conn_seg.npy>" << std::endl;
        std::exit(1);
    }

    std::string f_coor = argv[1];
    std::string f_conn_tri = argv[2];
    std::string f_conn_seg = argv[3];

    // Load mesh 

    auto coor_npy = cnpy::npy_load(f_coor);
    auto conn_tri_npy = cnpy::npy_load(f_conn_tri);
    auto conn_seg_npy = cnpy::npy_load(f_conn_seg);

    int num_points = coor_npy.shape[0];
    int dim = coor_npy.shape[1];

    int num_elmt_tri = conn_tri_npy.shape[0];
    int num_elmt_seg = conn_seg_npy.shape[0];

    int num_dof_src = num_elmt_tri;
    int num_dof_rcv = num_elmt_seg * dim;

    // std::cout << "Dimension = " << dim << std::endl; 
    // std::cout << "Number of nodes = " << num_points << std::endl; 
    std::cout << "num_elmt_tri = " << num_elmt_tri << std::endl; 
    std::cout << "num_elmt_seg = " << num_elmt_seg << std::endl; 
    // std::cout << "num_dof = " << num_dof << std::endl; 

    std::vector<double> coor_vec(coor_npy.data<double>(), coor_npy.data<double>() + coor_npy.num_vals);
    std::vector<int> conn_tri_vec(conn_tri_npy.data<int>(), conn_tri_npy.data<int>() + conn_tri_npy.num_vals);
    std::vector<int> conn_seg_vec(conn_seg_npy.data<int>(), conn_seg_npy.data<int>() + conn_seg_npy.num_vals);

    // // Print vectors 
    // std::cout << "conn = [";  
    // for (const auto& elem : conn_vec) {
    //     std::cout << elem << ", ";
    // }
    // std::cout << "]\n";

    int spatial_dim = 2;
    int nvertices_per_elt = 2;
    auto mesh_tri = CreateMeshFromVect<Triangle<0>>(
        /*spatial_dim*/ 2, 
        /*nvertices_per_elt*/ 3,
        coor_vec, conn_tri_vec
    );

    auto mesh_seg = CreateMeshFromVect<Segment<0>>(
        /*spatial_dim*/ 2, 
        /*nvertices_per_elt*/ 2,
        coor_vec, conn_seg_vec
    );

    // 2D segment piece-wise ct 0 element, H-kernel
    ElasticProperties elas(1.0, 0.25);
    using EltType = Segment<0>;

    // Kernel
    auto ker_obj = std::make_shared<BieElastostaticEigenstrain<Triangle<0>, Segment<0>, ElasticKernelType::V>>(elas, spatial_dim);

    // HR pattern 
    int max_leaf_size_ = 32;
    double eta = 4.0; 

    auto hr = HRepresentationRectangularMatrix(
        mesh_tri, // source mesh
        mesh_seg, // receiver mesh
        max_leaf_size_, eta, 
        /*verbose_*/ true, 
        /*homogeneous_size_pattern_*/ true
        );

    // Get matrix generator 
    BieMatrixGenerator<double> matgen(mesh_tri, mesh_seg, ker_obj, hr);

    // Get matrix generator by dof 
    BieMatrixGeneratorByDof<double> matgen_bydof{matgen};

    // // -----------
    // // Testing the BieMatrixGenerator by dof 

    // // Generate subset of matrix 
    // int b0 = 4;
    // int b1 = 4;
    // il::Array2D<double> M{4, 4}; // here the difs is 2x1 

    // matgen.set(
    //     b0, b1, il::io, M.Edit()
    // );

    // // switching to dof idices 
    // b0 *= 2;

    // std::cout << "Base matgen :\n";
    // prettyPrintArray2D(M);

    // il::Array2D<double> M_bydof{4, 4}; // here the difs is 2x2 

    // matgen_bydof.set(
    //     b0, b1, il::io, M_bydof.Edit()
    // );

    // std::cout << "By dof matgen 1 :\n";
    // prettyPrintArray2D(M_bydof);

    // // Only extractting the first row
    // il::Array2D<double> row{1, 4}; // here the difs is 2x2 

    // matgen_bydof.set(
    //     b0, b1, il::io, row.Edit()
    // );

    // std::cout << "By dof matgen row 0 :\n";
    // prettyPrintArray2D(row);

    // matgen_bydof.set(
    //     b0 + 1, b1, il::io, row.Edit()
    // );

    // std::cout << "By dof matgen row 1 :\n";
    // prettyPrintArray2D(row);


    // // Only extractting the first row
    // il::Array2D<double> column{4, 1};

    // matgen_bydof.set(
    //     b0, b1 +1, il::io, column.Edit()
    // );
    // std::cout << "By dof matgen column 1 :\n";
    // prettyPrintArray2D(column);

    // -----------
    // Debugging the ACA

    // int LR_block_i = 0;

    // il::int_t i0 = hr->pattern_.LRB_pattern(1, LR_block_i);
    // il::int_t j0 = hr->pattern_.LRB_pattern(2, LR_block_i);
    // il::int_t iend = hr->pattern_.LRB_pattern(3, LR_block_i);
    // il::int_t jend = hr->pattern_.LRB_pattern(4, LR_block_i);


    // il::Range range0{i0, iend};
    // il::Range range1{j0, jend};

    // double epsilon = 1e-4;

    // auto lra = adaptiveCrossApproximation<1>(matgen_bydof, range0, range1, epsilon, /*fixed_rank*/-1);

    // std::cerr << "ACA results on LR block of dof size (" << iend-i0 << ", " << jend-j0 << ") : A.size = (" << lra->A.size(0) << ", " << lra->A.size(1) << "), B.size = (" << lra->B.size(0) << ", " << lra->B.size(1) << ")\n";

    // // Extract the block data 

    // il::Array2D<double> LR_data{iend-i0, jend-j0};

    // matgen_bydof.set(
    //     i0, j0, il::io, LR_data.Edit()
    // );

    // const il::FullMatrix<double> matgen_from_dense{LR_data};

    // il::Range range0_dense{0, iend-i0};
    // il::Range range1_dense{0, jend-j0};

    // auto lra_from_dense = adaptiveCrossApproximation<1>(matgen_from_dense, range0_dense, range1_dense, epsilon, /*fixed_rank*/-1);

    // std::cerr << "[After export to dense] ACA results on LR block of dof size (" << iend-i0 << ", " << jend-j0 << ") : A.size = (" << lra->A.size(0) << ", " << lra->A.size(1) << "), B.size = (" << lra->B.size(0) << ", " << lra->B.size(1) << ")\n";

    // --------------
    // Filling a dense array

    // Fill dense array 
    il::Array2D<double> M{num_dof_rcv, num_dof_src};

    matgen.set(
        0, 0, il::io, M.Edit()
    );

    // std::cout << "M[0,0] = " << M(0,0) << "\n";

    // Transpose into row major
    std::vector<double> buffer(M.size(0) * M.size(1));
    for (int i = 0; i < M.size(0); i++) {
        for (int j = 0; j < M.size(1); j++) {
            buffer[i * M.size(1) + j] = M(i, j);
        }
    }

    
    // Save result to npy
    std::vector<size_t> shape = {
        static_cast<size_t>(M.size(0)),
        static_cast<size_t>(M.size(1))
    };    
    cnpy::npy_save("2dt0s0_v.npy", buffer.data(), shape, "w");

    std::cout << "Dense matrix saved in 2dt0s0_v.npy\n";

    // Also export permutations
    cnpy::npy_save("2dt0s0_v_perm_tri.npy", hr->permutation_1_.data(), {static_cast<size_t>(num_elmt_tri)}, "w");
    cnpy::npy_save("2dt0s0_v_perm_seg.npy", hr->permutation_0_.data(), {static_cast<size_t>(num_elmt_seg)}, "w");

    std::cout << "Permutations saved in 2dt0s0_v_perm_seg.npy and 2dt0s0_v_perm_tri.npy\n";



    return 0;
}