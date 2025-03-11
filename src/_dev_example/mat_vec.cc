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
#include <chrono>
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

  std::cout << "Reading mesh coordinates from " << f_coor << std::endl; 
  std::cout << "Reading mesh connectivity from " << f_conn << std::endl; 

  auto coor_npy = cnpy::npy_load(f_coor);
  auto conn_npy = cnpy::npy_load(f_conn);
  
  int num_points = coor_npy.shape[0];
  int dim = coor_npy.shape[1];
  int num_elemts = conn_npy.shape[0];
  int type_elemts = conn_npy.shape[1];
  int num_dof = num_elemts * dim;

  std::cout << "Dimension = " << dim << std::endl; 
  std::cout << "Number of nodes = " << num_points << std::endl; 
  std::cout << "Number of elements = " << num_elemts << std::endl; 
  std::cout << "Number of dof = " << num_dof << std::endl; 

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
  il::int_t max_leaf_size = 64;
  double eta = 3.0;
  double eps_aca = 1.e-3;

  // BigWhamIO(const std::vector<double> &coor, const std::vector<int> &conn,
  // const std::string &kernel, const std::vector<double> &properties, const int n_openMP_threads=8, const bool verbose=true)
  BigWhamIO hmat_io(coor_vec, conn_vec, kernel, properties); //, max_leaf_size, eta, eps_aca);

  //   void BuildPattern(const int max_leaf_size, const double eta);
  hmat_io.BuildPattern(max_leaf_size, eta);

  // void BuildHierarchicalMatrix(const int max_leaf_size, const double eta, const double eps_aca); // construct Hierarchical matrix
  hmat_io.BuildHierarchicalMatrix(max_leaf_size, eta, eps_aca);

  std::cout << "Succesfully built the hmat" << std::endl; 

  // MATVEC
  // std::vector<double> dd(num_dof, 1.0);

  il::Array<double> dd{num_dof, il::align_t(), 64};
  auto dd_edit = dd.Edit();
  for (int i(0); i<num_dof; i++){
    dd_edit[i] = 1.0;
  }

  auto dd_view = dd.view();

  #ifdef USE_ITT
  __itt_resume();
#endif

  int N_matvec = 100;
  auto start = std::chrono::high_resolution_clock::now(); 
  for (int i=0; i<N_matvec; i++){
    hmat_io.MatVec(dd_view);
  }
  auto end = std::chrono::high_resolution_clock::now(); 
  std::chrono::duration<double> duration = end - start; // Compute duration
  std::cout << "Matvec time = " << duration.count()/N_matvec << " seconds\n";
 
#ifdef USE_ITT
  __itt_pause();
#endif

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
