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
#include <iomanip>
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


std::string formatBytes(size_t bytes) {
  const double KB = 1024.0;
  const double MB = KB * 1024.0;
  const double GB = MB * 1024.0;

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2);

  if (bytes >= GB) {
      oss << (bytes / GB) << " GB";
  } else if (bytes >= MB) {
      oss << (bytes / MB) << " MB";
  } else if (bytes >= KB) {
      oss << (bytes / KB) << " KB";
  } else {
      oss << bytes << " B";
  }

  return oss.str();
}
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

  int max_leaf_size = 64;
  double eta = 5.0;
  int rank = 16;

  //   // Ask for max leaf size
  //   while (true) {
  //     std::cout << "Enter max leaf size (strictly positive integer): ";
  //     std::cin >> max_leaf_size;

  //     if (std::cin.fail() || max_leaf_size <= 0) {
  //         std::cin.clear(); // Clear the error flag
  //         std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore invalid input
  //         std::cout << "Invalid input. Please enter a strictly positive integer.\n";
  //     } else {
  //         break;
  //     }
  // }

  // // Ask for eta
  // while (true) {
  //     std::cout << "Enter eta (strictly positive double): ";
  //     std::cin >> eta;

  //     if (std::cin.fail() || eta <= 0.0) {
  //         std::cin.clear(); // Clear the error flag
  //         std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore invalid input
  //         std::cout << "Invalid input. Please enter a strictly positive double.\n";
  //     } else {
  //         break;
  //     }
  // }

  // // Ask for rank
  // while (true) {
  //     std::cout << "Enter rank (strictly positive integer, lower than max leaf size): ";
  //     std::cin >> rank;

  //     if (std::cin.fail() || rank <= 0 || rank >= max_leaf_size) {
  //         std::cin.clear(); // Clear the error flag
  //         std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore invalid input
  //         std::cout << "Invalid input. Please enter a strictly positive integer lower than max leaf size.\n";
  //     } else {
  //         break;
  //     }
  // }

  std::cout << "max_leaf_size = " << max_leaf_size << std::endl; 
  std::cout << "eta = " << eta << std::endl; 
  std::cout << "rank = " << rank << std::endl; 


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
  int n_GPUs = 999;
  bool verbose = false;
  bool homogeneous_size = true;
  bool use_Cuda = true;

  // auto start = std::chrono::high_resolution_clock::now();
  clock_gettime(CLOCK_MONOTONIC, &start);

  BigWhamIO hmat_io_2(coor_vec, conn_vec, kernel, properties, n_omp_threads, n_GPUs, verbose, homogeneous_size, use_Cuda, rank);
  hmat_io_2.BuildPattern(max_leaf_size, eta);


  size_t memory_required = hmat_io_2.GetGPUStorageRequirement();

  std::cout << "[GPU] memory required = " << formatBytes(memory_required) << std::endl; 

}
