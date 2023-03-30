#ifndef BIGWHAM_IO_GEN_H
#define BIGWHAM_IO_GEN_H

#include "core/be_mesh.h"
#include "core/bie_kernel.h"
#include "hmat/hmatrix/hmat.h"

namespace bie {

class BigWhamIOGen {
private:
  std::string kernel_name_;
  int spatial_dimension_;
  int dof_dimension_;

  Hmat<double> hmat_; // the  Hmat object
  std::shared_ptr<Mesh> mesh_src_;
  std::shared_ptr<Mesh> mesh_rec_;
  std::shared_ptr<BieKernel<double>> ker_obj_;

public:
  BigWhamIOGen(const std::vector<double> &coor, const std::vector<int> &conn,
               const std::string &kernel, const std::vector<double> &properties,
               const int max_leaf_size, const double eta, const double eps_aca);
  BigWhamIOGen(const std::vector<double> &coor_src,
               const std::vector<int> &conn_src,
               const std::vector<double> &coor_rec,
               const std::vector<int> &conn_rec, const std::string &kernel,
               const std::vector<double> &properties, const int max_leaf_size,
               const double eta, const double eps_aca);
  ~BigWhamIOGen();

  std::vector<double> GetCollocationPoints() const;
  std::vector<long> GetPermutation0();
  std::vector<long> GetPermutation1();
  std::vector<long> GetHPattern();
  void GetFullBlocks(std::vector<double> &val_list, std::vector<int> &pos_list);
  void GetDiagonal(std::vector<double> &val_list) ;
  std::vector<double> MatVec(const std::vector<double> &x) ;
  std::vector<double> MatVecPerm(const std::vector<double> &x) ;
  std::vector<double> ConvertToGlobal(const std::vector<double> &x_local) ;
  std::vector<double> ConvertToLocal(const std::vector<double> &x_global) ;
  double GetCompressionRatio();
  std::string kernel_name() { return kernel_name_; }
};

} // namespace bie

#endif // BIGWHAM_IO_GEN_H
