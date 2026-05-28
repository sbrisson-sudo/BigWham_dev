#include <memory>


#include "cnpy.h"
#include <gsl/gsl_errno.h>


#include <il/core/core.h>


// Bigwham imports
#include "core/be_mesh.h"
#include "core/be_mesh.h"
#include "io/bigwham_io_helper.h"
#include "io/bigwham_io.h"


#include "core/elastic_properties.h"
#include "elements/segment.h"

#include "core/bie_kernel.h"
#include "elasticity/bie_elastostatic.h"
#include "elasticity/bie_elastostatic_eigenstrain.h"

#include "hmat/hierarchical_representation.h"
#include "hmat/bie_matrix_generator.h"



using namespace bigwham;


int main(int argc, char * argv[]) {

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

    // Bigwham IO object
    BigWhamIO hmat_io(
        coor_vec,       // coor src
        conn_tri_vec,   // conn src
        coor_vec,       // coor rcv
        conn_seg_vec,   // conn rcv
        "2DT0-2DS0-V", 
        {1.0, 0.25}, 
        /*n_omp_threads*/ 1, 
        /*n_GPUs*/ -1, 
        /*verbose*/ true, 
        /*homogeneous_size*/ false, 
        /*use_Cuda*/ false
    );

    // build pattern 
    int max_leaf_size = 64;
    double eta = 6.0;
    hmat_io.BuildPattern(max_leaf_size, eta);

    // export patternhmat_io
    auto pattern_export = hmat_io.GetHPattern();

    // Export pattern to cnpy
    size_t n_cols = 6;  // 4 bounds + 1 flag + 1 value
    size_t n_rows = pattern_export.size() / n_cols;

    // Define shape (row-major)
    std::vector<size_t> shape = {n_rows, n_cols};

    // Save as NumPy .npy file
    cnpy::npy_save("2dt0s0_pattern.npy", pattern_export.data(), shape, "w");
    std::cout << "Hierachical pattern saved in 2dt0s0_pattern.npy\n";

    // Save the permutations
    cnpy::npy_save("2dt0s0_pattern_perm_tri.npy", hmat_io.GetPermutation().data(), {static_cast<size_t>(num_elmt_tri)}, "w");
    cnpy::npy_save("2dt0s0_pattern_perm_seg.npy",  hmat_io.GetPermutationReceivers().data(), {static_cast<size_t>(num_elmt_seg)}, "w");

    std::cout << "Permutations saved in 2dt0s0_pattern_perm_seg.npy and 2dt0s0_pattern_perm_tri.npy\n";

    return 0;

}