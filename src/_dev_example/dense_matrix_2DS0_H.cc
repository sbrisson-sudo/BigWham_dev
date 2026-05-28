#include <memory>
#include <iostream>
#include <iomanip> // For std::setw

#include "cnpy.h"

#include <il/core/core.h>


// Bigwham imports
#include "core/be_mesh.h"
#include "io/bigwham_io_helper.h"

#include "core/elastic_properties.h"
#include "elements/segment.h"

#include "core/bie_kernel.h"
#include "elasticity/bie_elastostatic.h"

#include "hmat/hierarchical_representation.h"
#include "hmat/bie_matrix_generator.h"
#include "hmat/bie_matrix_generator_by_dof.h"



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

    // get mesh file names form command line arguments
    if (argc < 3){
        std::cerr << "Usage : ./mat_vec <coor.npy> <conn.npy>" << std::endl;
        std::exit(1);
    }

    std::string f_coor = argv[1];
    std::string f_conn = argv[2];

    // Load mesh 

    auto coor_npy = cnpy::npy_load(f_coor);
    auto conn_npy = cnpy::npy_load(f_conn);

    int num_points = coor_npy.shape[0];
    int dim = coor_npy.shape[1];
    int num_elemts = conn_npy.shape[0];
    int type_elemts = conn_npy.shape[1];
    int num_dof = num_elemts * dim;

    // std::cout << "Dimension = " << dim << std::endl; 
    // std::cout << "Number of nodes = " << num_points << std::endl; 
    std::cout << "num_elemts = " << num_elemts << std::endl; 
    std::cout << "num_dof = " << num_dof << std::endl; 

    std::vector<double> coor_vec(coor_npy.data<double>(), coor_npy.data<double>() + coor_npy.num_vals);
    std::vector<int> conn_vec(conn_npy.data<int>(), conn_npy.data<int>() + conn_npy.num_vals);

    // // Print vectors 
    // std::cout << "conn = [";  
    // for (const auto& elem : conn_vec) {
    //     std::cout << elem << ", ";
    // }
    // std::cout << "]\n";

    int spatial_dim = 2;
    int nvertices_per_elt = 2;
    auto mesh = CreateMeshFromVect<Segment<0>>(
        spatial_dim, nvertices_per_elt,
        coor_vec, conn_vec
    );

    // 2D segment piece-wise ct 0 element, H-kernel
    ElasticProperties elas(1.0, 0.25);
    using EltType = Segment<0>;

    // Kernel
    auto ker_obj = std::make_shared<BieElastostatic<Segment<0>, Segment<0>, ElasticKernelType::H>>(elas, spatial_dim);

    // HR pattern 
    int max_leaf_size_ = 16;
    double eta = 0.0; 

    auto hr = HRepresentationSquareMatrix(
        mesh, max_leaf_size_, eta, 
        /*verbose_*/ true, 
        /*homogeneous_size_pattern_*/ false, 
        /*fixed_rank_*/ -1
    );

    // Get matrix generator 
    BieMatrixGenerator<double> matgen(mesh, mesh, ker_obj, hr);

    // // Generate subset of matrix 
    // int b0 = 0;
    // int b1 = 0;
    // il::Array2D<double> M{4, 4}; // here the difs is 2x2 

    // matgen.set(
    //     b0, b1, il::io, M.Edit()
    // );

    // std::cout << "Base matgen :\n";
    // prettyPrintArray2D(M);

    // // Get matrix generator by dof 
    // BieMatrixGeneratorByDof<double> matgen_bydof{matgen};

    // il::Array2D<double> M_bydof{4, 4}; // here the difs is 2x2 

    // matgen_bydof.set(
    //     b0, b1, il::io, M_bydof.Edit()
    // );

    // std::cout << "By dof matgen :\n";
    // prettyPrintArray2D(M_bydof);

    // // Only extractting the first row
    // il::Array2D<double> row{1, 4}; // here the difs is 2x2 

    // matgen_bydof.set(
    //     b0, b1, il::io, row.Edit()
    // );

    // std::cout << "By dof matgen row :\n";
    // prettyPrintArray2D(row);


    // Fill dense array 
    il::Array2D<double> M{num_dof, num_dof};

    matgen.set(
        0, 0, il::io, M.Edit()
    );

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
    cnpy::npy_save("2ds0_h.npy", buffer.data(), shape, "w");

    std::cout << "Dense matrix saved in 2ds0_h.npy\n";

    return 0;
}