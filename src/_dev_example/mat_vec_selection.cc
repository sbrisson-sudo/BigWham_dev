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
#include <ctime>
//#ifdef BIWGHAM_OPENMP
#include <omp.h>
//#endif

#ifdef USE_ITT
#include <ittnotify.h>
#endif
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
#include "hmat/hmatrix/hmat_selection.h"
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
  for (int i(1); i<num_dof; i++) dd_edit[i] = dd_edit[i-1] + 1/(static_cast<double>(num_dof)-1);
  auto dd_view = dd.view();

  // Timing
  const int N_matvec = 10;
  struct timespec start, end;
  double cpu_hmat_construction_time, cpu_matvec_time;
  // std::chrono::duration<double> cpu_hmat_construction_time;
  // std::chrono::duration<double> cpu_matvec_time;

  // LR approximation error
  double max_error_aca;

  // Normal hmat (CPU)
  int n_omp_threads = 999;
  int n_GPUs = 999;
  bool verbose = false;
  bool homogeneous_size = false;
  // bool homogeneous_size = true;
  bool use_Cuda = false;
  std::cout << "[CPU] eps_aca = " << eps_aca << std::endl; 

  // auto start = std::chrono::high_resolution_clock::now();
  clock_gettime(CLOCK_MONOTONIC, &start);

  BigWhamIO hmat_io(coor_vec, conn_vec, kernel, properties, n_omp_threads, n_GPUs, verbose, homogeneous_size, use_Cuda);
  hmat_io.BuildPattern(max_leaf_size, eta);
  hmat_io.BuildHierarchicalMatrix(max_leaf_size, eta, eps_aca);

  clock_gettime(CLOCK_MONOTONIC, &end);
  cpu_hmat_construction_time = (end.tv_sec - start.tv_sec) +
                      (end.tv_nsec - start.tv_nsec) / 1e9;
  std::cout << "[CPU] Hmat construction time = " << cpu_hmat_construction_time*1000 << " ms" << std::endl;

  // auto end = std::chrono::high_resolution_clock::now();
  // cpu_hmat_construction_time = end - start;
  // std::cout << "[CPU] Hmat construction time = " << cpu_hmat_construction_time.count() << " ms" << std::endl;


  max_error_aca = hmat_io.GetMaxErrorACA();
  std::cout << "[CPU] Max ACA error = " << max_error_aca << std::endl;

  // Compute matvec
  t = hmat_io.MatVec(dd_view); 

  // start = std::chrono::high_resolution_clock::now(); 
  clock_gettime(CLOCK_MONOTONIC, &start);

  for (int i(0); i<N_matvec; i++) t = hmat_io.MatVec(dd_view);

  clock_gettime(CLOCK_MONOTONIC, &end);
  cpu_matvec_time = (end.tv_sec - start.tv_sec) +
  (end.tv_nsec - start.tv_nsec) / 1e9;
  std::cout << "[CPU] Hmat matvec time = " << cpu_matvec_time/N_matvec*1000 << " ms" << std::endl;

  // end = std::chrono::high_resolution_clock::now();
  // cpu_matvec_time = end - start;
  // std::cout << "[CPU] Matvec time = " << cpu_matvec_time.count()/N_matvec << " ms" << std::endl;

  // cnpy::npy_save("result_cpu.npy", t.data(), {static_cast<size_t>(num_dof)}, "w");
  
  // Compute l2 norm
  // Compute l2 norm
  double l2_norm = 0;
  double res_sum = 0;
  auto t_view = t.view();
  for (int i(0); i<num_dof; i++){
    l2_norm += t_view[i] * t_view[i];
    res_sum += t_view[i];
  } 
  l2_norm = std::sqrt(l2_norm);
  std::cout << "[CPU] L2 norm of the product of H with [0, 1/dof, ...,  1] = " << l2_norm << std::endl;
  std::cout << "[CPU] sum of the product of H with [0, 1/dof, ...,  1] = " << res_sum  << std::endl;

  cnpy::npy_save("mat_vec_res.npy", t_view.data(), {static_cast<size_t>(t_view.size())}, "w");

  // Here we perform matrix selection
  bigwham::Hmat<double>& hmat_base = *hmat_io.getHmat();

  // Selection 
  il::Array<int> row_indices(3);
  il::Array<int> col_indices(3);

  auto row_indices_edit = row_indices.Edit();
  row_indices_edit[0] = 1;
  row_indices_edit[1] = 2;
  row_indices_edit[2] = 3;
  auto col_indices_edit = col_indices.Edit();
  col_indices_edit[0] = 4;
  col_indices_edit[1] = 5;
  col_indices_edit[2] = 6;

  HmatSelection<double> hmat_selected(hmat_base, row_indices, col_indices);

  // Now we try a matvec operation 
  il::Array<double> x_selec{static_cast<il::int_t>(col_indices.size()*3), 1.0};
  il::Array<double> y_selec = hmat_selected.matvecOriginal(x_selec.view());

  // WE compare the result to the full matvec
  il::Array<double> x_full{hmat_base.size(0), 0.0};
  auto x_full_edit = x_full.Edit();
  for (auto i : col_indices){
    for (int j=0; j<3; j++){
      x_full_edit[i*3+j] = 1.0;
    }
  }

  il::Array<double> y_full = hmat_base.matvecOriginal(x_full.view());

  // Then we extract the elemnets we need 
  auto y_full_view = y_full.view();
  il::Array<double> y_selec_from_full{static_cast<il::int_t>(row_indices.size()*3)};
  auto y_selec_from_full_edit = y_selec_from_full.Edit();

  for (int i=0; i<row_indices.size(); i++){
    for (int j=0; j<3; j++){
      y_selec_from_full_edit[i*3+j] = y_full_view[row_indices[i]*3+j];
    }
  }
  // COmpute l2 nomr of res
  double sum_full = 0.0;
  double sum_selec = 0.0;

  auto y_selec_view = y_selec.view();
  auto y_selec_from_full_view = y_selec_from_full.view();

  for (int i=0; i<y_selec.size(); i++) sum_selec += y_selec_view[i];

  for (int i=0; i<y_selec_from_full.size(); i++) sum_full += y_selec_from_full_view[i];

  std::cout << "[CPU] sum of result using selection = " << sum_selec << std::endl;
  std::cout << "[CPU] sum of result using full hmat = " << sum_full << std::endl;

  // Now a test using bigwhamIO rather than hmat 
  BigWhamIO hmat_io_selec = hmat_io.hmatSelection(row_indices, col_indices);
  y_selec = hmat_io_selec.MatVec(x_selec.view());

  // Trying to get the FR blocks data
  il::Array<double> val_list;
  il::Array<int> pos_list;
  hmat_io_selec.GetFullBlocks(val_list, pos_list);

  std::cout << "Call to hmat_io_selec.GetFullBlocks successfull " <<  std::endl;


  // std::cout << "[CPU] res = [" ;
  // for (int i(0); i<num_dof; i++) std::cout << t_view[i] << ", ";
  // std::cout << "\n";

  // // Exporting stuff
  // auto hmat = hmat_io.getHmat();
  // auto hmat_cuda = std::dynamic_pointer_cast<bigwham::HmatCuda<double>>(hmat);

  // // Test getting blocks
  // const int block_id = 0;
  // auto block_host = hmat_cuda->getFRBlockDataHost(block_id);
  // size_t block_size = static_cast<size_t>(block_host.size(0));
  // cnpy::npy_save("full_block_host.npy", block_host.data(), {block_size, block_size}, "w");

  // block_host = hmat_cuda->getFRBlockDataDevice(block_id);
  // cnpy::npy_save("full_block_device.npy", block_host.data(), {block_size, block_size}, "w");

  // // Test getting rowPtr and colInd
  // auto rowPtr = hmat_cuda->getFRBlockRowPtrHost();
  // cnpy::npy_save("row_ptr.npy", rowPtr.data(), {static_cast<size_t>(rowPtr.size())}, "w");
  // auto colInd = hmat_cuda->getFRBlockColIndHost();
  // cnpy::npy_save("col_ind.npy", colInd.data(), {static_cast<size_t>(colInd.size())}, "w");



//   #ifdef USE_ITT
//   __itt_resume();
// #endif

  // int N_matvec = 1;
  // auto start = std::chrono::high_resolution_clock::now(); 
  // for (int i=0; i<N_matvec; i++){
  //   t = hmat_io.MatVec(dd_view);
  // }
  // auto end = std::chrono::high_resolution_clock::now(); 
  // std::chrono::duration<double> duration = end - start; // Compute duration
  // std::cout << "Matvec time = " << duration.count()/N_matvec << " seconds\n";
 
// #ifdef USE_ITT
//   __itt_pause();
// #endif



  // SquareMatrixGenerator<double> M(my_mesh, ker, hr);

  // il::Array<double> dd{, 0.0};
  // il::Array<double> dd_perm{M.size(1), 0.0};
  // il::Array<double> trac{M.size(0), 0.0};
  // il::Array<double> trac_perm{M.size(0), 0.0};

  // double pre_fac = (8 * (1 - nu * nu)) / (pi * E);

  // double rsq;
  // for (il::int_t i = 0; i < M.sizeAsBlocks(1); i++) {
  //   rsq = xcol(i, 0) * xcol(i, 0) + xcol(i, 1) * xcol(i, 1);
  //   dd[dim * i + 2] = pre_fac * std::sqrt(radius * radius - rsq);
  // }

  // for (il::int_t i = 0; i < M.sizeAsBlocks(1); i++) {
  //   il::int_t j = hr->permutation_1_[i];
  //   for (uint d = 0; d < dim; d++) {
  //     dd_perm[dim * i + d] = dd[dim * j + d];
  //   }
  // }

  // double y0 = 0.;
  // auto start = omp_get_wtime();
  // // for (il::int_t i = 0; i < 1000; ++i) {
  // trac_perm = hmat.matvec(dd_perm.view());
  // y0 += trac_perm[0];
  // //}
  // auto end = omp_get_wtime();
  // std::cout << "Hmat matvec: " << (end - start) / 1000 << "s - y0: " << y0
  //           << std::endl;

  // for (il::int_t i = 0; i < M.sizeAsBlocks(0); i++) {
  //   il::int_t j = hr->permutation_0_[i];
  //   for (uint d = 0; d < dim; d++) {
  //     trac[dim * j + d] = trac_perm[dim * i + d];
  //   }
  // }

  // il::Array<double> rel_err{M.sizeAsBlocks(0), 0.};
  // for (il::int_t i = 0; i < M.sizeAsBlocks(0); i++) {
  //   rel_err[i] += trac[dim * i + 0] * trac[dim * i + 0];
  //   rel_err[i] += trac[dim * i + 1] * trac[dim * i + 1];
  //   rel_err[i] +=
  //       (trac[dim * i + 2] - pressure) * (trac[dim * i + 2] - pressure);
  //   rel_err[i] = std::sqrt(rel_err[i]);
  // }

  // std::cout << "Mean rel error " << il::mean(rel_err) << std::endl;
  // std::cout << "L2 rel error "
  //           << il::norm(rel_err, il::Norm::L2) / M.sizeAsBlocks(0) << std::endl;

  // return 0;
}
