///=============================================================================
///
/// \author      Ankit
///
/// \copyright   Copyright (©) 2018 EPFL (Ecole Polytechnique Fédérale
///              de Lausanne)\n
///              Geo-Energy lab
///
/// \brief       Test for Penny shape crack for profiling bigwham
///
///=============================================================================

/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <cmath>
#include <il/Array2D.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <string>
#include <vector>
// #include <chrono>
#include <ctime>

//#ifdef BIWGHAM_OPENMP
#include <omp.h>
//#endif

#ifdef USE_ITT
#include <ittnotify.h>
#endif

#include "cuda_profiler_api.h"

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "cnpy.h"
#include "npy_tools.h"
/* -------------------------------------------------------------------------- */
#include "core/be_mesh.h"
#include "core/bie_kernel.h"
#include "core/elastic_properties.h"

#include "elasticity/bie_elastostatic.h"
#include "elements/triangle.h"
#include "hmat/hierarchical_representation.h"
#include "hmat/hmatrix/hmat.h"
#include "hmat/square_matrix_generator.h"
#include "io/bigwham_io.h"
/* -------------------------------------------------------------------------- */

using namespace bigwham;
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {

#ifdef USE_ITT
  __itt_pause();
#endif

  if (argc < 3){
    std::cerr << "Usage : ./mat_vec <coor.npy> <conn.npy>" << std::endl;
    std::exit(1);
  }

  std::string f_coor = argv[1];
  std::string f_conn = argv[2];

  // std::cout << "Reading mesh coordinates from " << f_coor << std::endl; 
  // std::cout << "Reading mesh connectivity from " << f_conn << std::endl; 

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

  std::string kernel = "3DT0-H";

  if (dim != 3){
    std::cerr << "We are using the 3DT0-H kernel, mesh has to be 3D" << std::endl;
    std::exit(1);
  }
  if (type_elemts != 3){
    std::cerr << "We are using the 3DT0-H kernel, elements have to be triangles" << std::endl;
    std::exit(1);
  }

  // Elastic properties
  double G = 1.0;
  double nu = 0.25;
  double E = (2 * G) * (1 + nu);
  std::vector<double> properties = {G, nu};

  // H-mat parameters
  il::int_t max_leaf_size = 32;
  double eta = 3.0;
  double eps_aca = 1.e-3;

  std::cout << "max_leaf_size = " << max_leaf_size << std::endl; 
  std::cout << "eta = " << eta << std::endl; 

  // matvec vectors
  il::Array<double> dd{num_dof, il::align_t(), 64};
  il::Array<double> t{num_dof, il::align_t(), 64};
  auto dd_edit = dd.Edit();
  dd_edit[0] = 0;
  for (int i(1); i<num_dof; i++) dd_edit[i] = dd_edit[i-1] + 1/static_cast<double>(num_dof);
  auto dd_view = dd.view();

  // Timing
  const int N_matvec = 10;
  struct timespec start, end;
  double gpu_hmat_construction_time, gpu_matvec_time;
  // std::chrono::duration<double> gpu_hmat_construction_time;
  // std::chrono::duration<double> gpu_matvec_time;

  // LR approximation error
  double max_error_aca;
  
  // GPU matvec
  int n_omp_threads = 999;
  bool verbose = true;
  bool homogeneous_size = true;
  bool use_Cuda = true;
  const int rank = 15;
  std::cout << "[GPU] fixed rank = " << rank << std::endl; 

  // auto start = std::chrono::high_resolution_clock::now();
  clock_gettime(CLOCK_MONOTONIC, &start);

  BigWhamIO hmat_io_2(coor_vec, conn_vec, kernel, properties, n_omp_threads, verbose, homogeneous_size, use_Cuda, rank);
  hmat_io_2.BuildPattern(max_leaf_size, eta);
  hmat_io_2.BuildHierarchicalMatrix(max_leaf_size, eta, eps_aca);

  clock_gettime(CLOCK_MONOTONIC, &end);
  gpu_hmat_construction_time = (end.tv_sec - start.tv_sec) +
                      (end.tv_nsec - start.tv_nsec) / 1e9;
  std::cout << "[GPU] Hmat construction time = " << gpu_hmat_construction_time*1000 << " ms" << std::endl;

  // auto end = std::chrono::high_resolution_clock::now();
  // gpu_hmat_construction_time = end - start;
  // std::cout << "[CPU] Hmat construction time = " << gpu_hmat_construction_time.count() << " ms" << std::endl;


  max_error_aca = hmat_io_2.GetMaxErrorACA();
  std::cout << "[GPU] Max ACA error = " << max_error_aca << std::endl;

  t = hmat_io_2.MatVec(dd_view); 


  // Compute matvec

  cudaProfilerStart();  

  t = hmat_io_2.MatVec(dd_view); 

  cudaProfilerStop();  

  
  // start = std::chrono::high_resolution_clock::now(); 
  clock_gettime(CLOCK_MONOTONIC, &start);

  for (int i(0); i<N_matvec; i++) t = hmat_io_2.MatVec(dd_view);

  clock_gettime(CLOCK_MONOTONIC, &end);
  gpu_matvec_time = (end.tv_sec - start.tv_sec) +
  (end.tv_nsec - start.tv_nsec) / 1e9;
  std::cout << "[GPU] Hmat matvec time = " << gpu_matvec_time/N_matvec*1000 << " ms" << std::endl;
  
  // end = std::chrono::high_resolution_clock::now();
  // gpu_matvec_time = end - start;
  // std::cout << "[GPU] Matvec time = " << gpu_matvec_time.count()/N_matvec << " ms" << std::endl;

  // cnpy::npy_save("result_gpu.npy", t.data(), {static_cast<size_t>(num_dof)}, "w");
  
  // Compute l2 norm
  double l2_norm = 0;
  auto t_view = t.view();
  for (int i(0); i<num_dof; i++) l2_norm += t_view[i] * t_view[i];
  l2_norm = std::sqrt(l2_norm);
  std::cout << "[GPU] L2 norm of the product of H with [0, 1/dof, ...,  1] = " << l2_norm << std::endl;

}
