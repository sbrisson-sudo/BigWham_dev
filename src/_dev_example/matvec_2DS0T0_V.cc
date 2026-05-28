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
    BigWhamIO hmat_io_dense(
        coor_vec,       // coor src
        conn_tri_vec,   // conn src
        coor_vec,       // coor rcv
        conn_seg_vec,   // conn rcv
        "2DT0-2DS0-V", 
        {1.0, 0.25}, 
        /*n_omp_threads*/ 10, 
        /*n_GPUs*/ -1, 
        /*verbose*/ true, 
        /*homogeneous_size*/ false, 
        /*use_Cuda*/ false
    );

    // build the hmat without compression
    int max_leaf_size = 64;
    double eta = 0.0;
    double eps_aca = 1e-4;
    hmat_io_dense.BuildHierarchicalMatrix(max_leaf_size, eta, eps_aca);

    // rebuilding the hmat with compression
    BigWhamIO hmat_io(
        coor_vec,       // coor src
        conn_tri_vec,   // conn src
        coor_vec,       // coor rcv
        conn_seg_vec,   // conn rcv
        "2DT0-2DS0-V", 
        {1.0, 0.25}, 
        /*n_omp_threads*/ 10, 
        /*n_GPUs*/ -1, 
        /*verbose*/ true, 
        /*homogeneous_size*/ false, 
        /*use_Cuda*/ false
    );

    eta = 5.0;
    hmat_io.BuildHierarchicalMatrix(max_leaf_size, eta, eps_aca);

    // Compute matvec

    // matvec vectors
    il::Array<double> eps0{num_dof_src, il::align_t(), 64};
    auto eps0_edit = eps0.Edit();
    eps0_edit[0] = 0;
    for (int i(1); i<num_dof_src; i++) eps0_edit[i] = eps0_edit[i-1] + 1/(static_cast<double>(num_dof_src)-1);
    auto eps0_view = eps0.view();


    // With uncompressed hmat    
    il::Array<double> t{num_dof_rcv, il::align_t(), 64};
    t = hmat_io_dense.MatVec(eps0_view);
    
    double l2_norm = 0;
    double res_sum = 0;
    auto t_view = t.view();
    for (int i(0); i<num_dof_rcv; i++){
        l2_norm += t_view[i] * t_view[i];
        res_sum += t_view[i];
    } 
    l2_norm = std::sqrt(l2_norm);
    std::cout << "[Uncompressed] L2 norm of the product of H with [0, 1/dof, ...,  1] = " << l2_norm << std::endl;
    std::cout << "[Uncompressed] sum of the product of H with [0, 1/dof, ...,  1] = " << res_sum  << std::endl;

    // With compressed hmat
    il::Array<double> t_b{num_dof_rcv, il::align_t(), 64};
    t_b = hmat_io.MatVec(eps0_view);
    
    l2_norm = 0;
    res_sum = 0;
    auto t_b_view = t_b.view();
    for (int i(0); i<num_dof_rcv; i++){
        l2_norm += t_b_view[i] * t_b_view[i];
        res_sum += t_b_view[i];
    } 
    l2_norm = std::sqrt(l2_norm);
    std::cout << "[Compressed] L2 norm of the product of H with [0, 1/dof, ...,  1] = " << l2_norm << std::endl;
    std::cout << "[Compressed] sum of the product of H with [0, 1/dof, ...,  1] = " << res_sum  << std::endl;


    return 0;

}