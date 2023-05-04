#include "hmat.h"

#if defined(IL_OPENMP)
#include <omp.h>
#endif

#if defined(BIGWHAM_HDF5)
#include "hdf5.h"
#endif

namespace bie {
/* -------------------------------------------------------------------------- */

template <typename T>
void Hmat<T>::fullBlocksOriginal(il::io_t, il::Array<T> &val_list,
                                 il::Array<int> &pos_list) {
  // return the full blocks in the permutted Original dof state
  // in the val_list and pos_list 1D arrays.
  IL_EXPECT_FAST(isBuilt_FR_);
  IL_EXPECT_FAST(hr_->permutation_1_.size() * dof_dimension_ == size_[1]);
  //  compute the number of  entries in the whole full rank blocks
  int nbfentry = 0;
  for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
    il::Array2DView<double> aux = ((*full_rank_blocks_[i]).view());
    nbfentry = nbfentry + static_cast<int>(aux.size(0) * aux.size(1));
  }

  // prepare outputs
  pos_list.Resize(nbfentry * 2);
  val_list.Resize(nbfentry);

  il::Array<int> permutDOF{dof_dimension_ * hr_->permutation_1_.size(), 0};
  IL_EXPECT_FAST(permutDOF.size() == size_[0]);
  for (il::int_t i = 0; i < hr_->permutation_1_.size(); i++) {
    for (il::int_t j = 0; j < dof_dimension_; j++) {
      permutDOF[i * dof_dimension_ + j] =
          hr_->permutation_1_[i] * dof_dimension_ + j;
    }
  }
  // loop on full rank and get i,j and val
  int nr = 0;
  int npos = 0;
  for (il::int_t k = 0; k < hr_->pattern_.n_FRB; k++) {
    il::Array2D<double> aux = *full_rank_blocks_[k];
    il::int_t i0 = hr_->pattern_.FRB_pattern(1, k);
    il::int_t j0 = hr_->pattern_.FRB_pattern(2, k);
    il::int_t index = 0;
    for (il::int_t j = 0; j < aux.size(1); j++) {
      for (il::int_t i = 0; i < aux.size(0); i++) {
        pos_list[npos + 2 * index] = permutDOF[(i + dof_dimension_ * i0)];
        pos_list[npos + 2 * index + 1] = permutDOF[(j + dof_dimension_ * j0)];
        val_list[nr + index] = aux(i, j);
        index++;
      }
    }
    nr = nr + static_cast<int>(aux.size(0) * aux.size(1));
    npos = npos + static_cast<int>(2 * aux.size(0) * aux.size(1));
  }
  std::cout << "Done Full Block: nval " << val_list.size() << " / "
            << pos_list.size() / 2 << " n^2 "
            << (this->size_[0]) * (this->size_[1]) << "\n";
}
/* -------------------------------------------------------------------------- */

template <typename T> std::vector<T> Hmat<T>::diagonalOriginal() {
  // return diagonal in original state....
  il::int_t diag_size = il::max(size_[0], size_[1]);
  il::int_t ncolpoints = diag_size / dof_dimension_;
  std::vector<T> diag_raw = this->diagonal();
  std::vector<T> diag(diag_raw.size(), 0.);
  // permut back
  for (il::int_t i = 0; i < ncolpoints; i++) {
    for (int j = 0; j < dof_dimension_; j++) {
      diag[dof_dimension_ * hr_->permutation_1_[i] + j] =
          diag_raw[dof_dimension_ * i + j];
    }
  }
  return diag;
}
/* -------------------------------------------------------------------------- */

template <typename T> std::vector<T> Hmat<T>::diagonal() {
  // return diagonal in permutted state....
  IL_EXPECT_FAST(isBuilt_FR_);
  il::int_t diag_size = il::max(size_[0], size_[1]);
  std::vector<T> diag(static_cast<long>(diag_size), 0.);

  for (il::int_t k = 0; k < hr_->pattern_.n_FRB; k++) {
    il::Array2D<double> aux = *full_rank_blocks_[k];
    il::int_t i0 = hr_->pattern_.FRB_pattern(1, k);
    il::int_t j0 = hr_->pattern_.FRB_pattern(2, k);
    il::int_t iend = hr_->pattern_.FRB_pattern(3, k);
    il::int_t jend = hr_->pattern_.FRB_pattern(4, k);
    // check if it intersects the diagonal
    //
    bool in_lower = (i0 > j0) && (i0 > jend) && (iend > j0) && (iend > jend);
    bool in_upper = (i0 < j0) && (i0 < jend) && (iend < j0) && (iend < jend);
    if ((!in_lower) && (!in_upper)) // this fb intersect the diagonal....
    {
      for (il::int_t j = 0; j < aux.size(1); j++) {
        for (il::int_t i = 0; i < aux.size(0); i++) {
          if ((i + dof_dimension_ * i0) ==
              (j + dof_dimension_ * j0)) { // on diagonal !
            diag[(i + dof_dimension_ * i0)] = aux(i, j);
          }
        }
      }
    }
  }
  return diag;
}
/* -------------------------------------------------------------------------- */

//   // simple constructor from pattern  -> todo remove
//   Hmat(const bie::HPattern& pattern){
//     pattern_=pattern;
//   };

// direct constructor
template <typename T>
Hmat<T>::Hmat(const bie::MatrixGenerator<T> &matrix_gen,
              const double epsilon_aca) {
  // construction directly
  this->toHmat(matrix_gen, epsilon_aca);
}
/* -------------------------------------------------------------------------- */

template <typename T>
Hmat<T>::Hmat(const std::string & filename){
  // construction directly
  il::Timer tt;
  tt.Start();
  this->readFromFile(filename);
  tt.Stop();
  std::cout << "Reading of hmat done in "  << tt.time() <<"\n";std::cout << "Compression ratio - " << this->compressionRatio() <<"\n";
  std::cout << "Hmat object - built " << "\n";
}
/* -------------------------------------------------------------------------- */

// Main constructor
template <typename T>
void Hmat<T>::toHmat(const MatrixGenerator<T> &matrix_gen,
                     const double epsilon_aca) {
  this->hr_ = matrix_gen.hr();
  // construction directly
  il::Timer tt;
  tt.Start();
  this->build(matrix_gen, epsilon_aca);
  tt.Stop();
  std::cout << "Creation of hmat done in " << tt.time() << "\n";
  std::cout << "Compression ratio - " << this->compressionRatio() << "\n";
  std::cout << "Hmat object - built "
            << "\n";
}
/* -------------------------------------------------------------------------- */

// delete the memory pointed by low_rank_blocks_ and  full_rank_blocks_
template <typename T> void Hmat<T>::hmatMemFree() {
  for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
    this->full_rank_blocks_[i].reset();
  }
  for (il::int_t i = 0; i < hr_->pattern_.n_LRB; i++) {
    this->low_rank_blocks_[i].reset();
  }
}
/* -------------------------------------------------------------------------- */

// getting the nb of entries of the hmatrix
template <typename T> il::int_t Hmat<T>::nbOfEntries() {
  IL_EXPECT_FAST(isBuilt_);
  il::int_t n = 0;
  for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
    il::Array2DView<double> a = (*full_rank_blocks_[i]).view();
    n += a.size(0) * a.size(1);
  }
  for (il::int_t i = 0; i < hr_->pattern_.n_LRB; i++) {
    il::Array2DView<double> a = (*low_rank_blocks_[i]).A.view();
    il::Array2DView<double> b = (*low_rank_blocks_[i]).B.view();
    n += a.size(0) * a.size(1) + b.size(0) * b.size(1);
  }
  return n;
}
/* -------------------------------------------------------------------------- */

template <typename T>
void Hmat<T>::buildFR(const bie::MatrixGenerator<T> &matrix_gen) {
  // construction of the full rank blocks
  std::cout << " Loop on full blocks construction  \n";
  std::cout << " N full blocks " << hr_->pattern_.n_FRB << " \n";

  full_rank_blocks_.resize(hr_->pattern_.n_FRB);
#pragma omp parallel for schedule(static, frb_chunk_size_)
  for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
    il::int_t i0 = hr_->pattern_.FRB_pattern(1, i);
    il::int_t j0 = hr_->pattern_.FRB_pattern(2, i);
    il::int_t iend = hr_->pattern_.FRB_pattern(3, i);
    il::int_t jend = hr_->pattern_.FRB_pattern(4, i);

    const il::int_t ni = matrix_gen.blockSize() * (iend - i0);
    const il::int_t nj = matrix_gen.blockSize() * (jend - j0);

    std::unique_ptr<il::Array2D<T>> a =
        std::make_unique<il::Array2D<T>>(ni, nj);
    matrix_gen.set(i0, j0, il::io, a->Edit());
    full_rank_blocks_[i] = std::move(a);
  }
  isBuilt_FR_ = true;
}
/* -------------------------------------------------------------------------- */

/// \param matrix_gen
/// \param epsilon
template <typename T>
template <il::int_t dim>
void Hmat<T>::buildLR(const bie::MatrixGenerator<T> &matrix_gen,
                      const double epsilon) {
  // constructing the low rank blocks
  dof_dimension_ = matrix_gen.blockSize();
  std::cout << "Loop on low rank blocks construction\n";
  std::cout << "N low rank blocks " << hr_->pattern_.n_LRB << "\n";
  std::cout << "dof_dimension: " << dof_dimension_ << "\n";

#if defined(_OPENMP)
  auto nthreads = omp_get_max_threads();
#endif
  low_rank_blocks_.resize(hr_->pattern_.n_LRB);
#pragma omp parallel for schedule(static, lrb_chunk_size_)
  for (il::int_t i = 0; i < hr_->pattern_.n_LRB; i++) {
    il::int_t i0 = hr_->pattern_.LRB_pattern(1, i);
    il::int_t j0 = hr_->pattern_.LRB_pattern(2, i);
    il::int_t iend = hr_->pattern_.LRB_pattern(3, i);
    il::int_t jend = hr_->pattern_.LRB_pattern(4, i);
    il::Range range0{i0, iend};
    il::Range range1{j0, jend};

    // we need 7a LRA generator virtual template similar to the Matrix
    // generator... here we have an if condition for the LRA call dependent on
    // dof_dimension_
    auto lra = bie::adaptiveCrossApproximation<dim>(matrix_gen, range0, range1,
                                                    epsilon);

    // store the rank in the low_rank pattern
    hr_->pattern_.LRB_pattern(5, i) = lra->A.size(1);
    low_rank_blocks_[i] =
        std::move(lra); // lra_p does not exist after such call
  }
  isBuilt_LR_ = true;
}
/* -------------------------------------------------------------------------- */

// filling up the h-matrix sub-blocks
template <typename T>
void Hmat<T>::build(const bie::MatrixGenerator<T> &matrix_gen,
                    const double epsilon) {
  dof_dimension_ = matrix_gen.blockSize();
  size_[0] = matrix_gen.size(0);
  size_[1] = matrix_gen.size(1);
  buildFR(matrix_gen);

  switch (matrix_gen.blockSize()) {
  case 1:
    buildLR<1>(matrix_gen, epsilon);
    break;
  case 2:
    buildLR<2>(matrix_gen, epsilon);
    break;
  case 3:
    buildLR<3>(matrix_gen, epsilon);
    break;
  }

  isBuilt_ = isBuilt_FR_ && isBuilt_LR_;
}
/* -------------------------------------------------------------------------- */

// getting the compression ratio of the hmatrix (double which is <=1)
template <typename T> double Hmat<T>::compressionRatio() {
  auto nb_elts = static_cast<double>(nbOfEntries());
  return nb_elts / static_cast<double>(size_[0] * size_[1]);
}
/* -------------------------------------------------------------------------- */

// H-Matrix vector multiplication without permutation
// in - il:Array<T>
// out - il:Array<T>
template <typename T> il::Array<T> Hmat<T>::matvec(il::ArrayView<T> x) {
  IL_EXPECT_FAST(x.size() == size_[1]);
  il::Array<T> y{size_[0], 0.};

#if defined(_OPENMP)
  auto nthreads = omp_get_max_threads();
  il::Array2D<T> yprivate_storage{size_[0], nthreads};
#endif

#pragma omp parallel shared(y)
  {
    il::Array<T> yprivate{size_[0], 0.};

#pragma omp for nowait schedule(static, frb_chunk_size_)
    for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
      auto i0 = hr_->pattern_.FRB_pattern(1, i);
      auto j0 = hr_->pattern_.FRB_pattern(2, i);
      auto iend = hr_->pattern_.FRB_pattern(3, i);
      auto jend = hr_->pattern_.FRB_pattern(4, i);

      auto a = (*full_rank_blocks_[i]).view();
      auto xs = x.view(il::Range{j0 * dof_dimension_, jend * dof_dimension_});
      auto ys =
          yprivate.Edit(il::Range{i0 * dof_dimension_, iend * dof_dimension_});

      il::blas(1.0, a, xs, 1.0, il::io, ys);
    }

#pragma omp for nowait schedule(static, lrb_chunk_size_)
    for (il::int_t ii = 0; ii < hr_->pattern_.n_LRB; ii++) {
      auto i0 = hr_->pattern_.LRB_pattern(1, ii);
      auto j0 = hr_->pattern_.LRB_pattern(2, ii);
      auto iend = hr_->pattern_.LRB_pattern(3, ii);
      auto jend = hr_->pattern_.LRB_pattern(4, ii);

      auto a = low_rank_blocks_[ii]->A.view();
      auto b = low_rank_blocks_[ii]->B.view();

      auto xs = x.view(il::Range{j0 * dof_dimension_, jend * dof_dimension_});
      auto ys =
          yprivate.Edit(il::Range{i0 * dof_dimension_, iend * dof_dimension_});
      auto r = a.size(1);
      il::Array<double> tmp{r, 0.0};

      il::blas(1.0, b, il::Dot::None, xs, 0.0, il::io,
               tmp.Edit()); // Note here we have stored b (not b^T)
      il::blas(1.0, a, tmp.view(), 1.0, il::io, ys);
    }

#pragma omp critical
    il::blas(1., yprivate.view(), il::io_t{}, y.Edit());
  }

  return y;
}
/* -------------------------------------------------------------------------- */

// H-Matrix vector multiplication with permutation for rectangular matrix
// cases (only 1 permutation) in & out as std::vector todo : write another one
// for the case of 2 permutations (rect. mat cases (for source != receivers)
template <typename T> il::Array<T> Hmat<T>::matvecOriginal(il::ArrayView<T> x) {

  il::Array<T> z{static_cast<il::int_t>(x.size())};
  // permutation of the dofs according to the re-ordering sue to clustering
  il::int_t ncolpoints = this->size(1) / dof_dimension_;
  il::int_t nrowpoints = this->size(0) / dof_dimension_;
  for (il::int_t i = 0; i < ncolpoints; i++) {
    for (int j = 0; j < dof_dimension_; j++) {
      z[dof_dimension_ * i + j] =
          x[dof_dimension_ * hr_->permutation_1_[i] + j];
    }
  }
  il::Array<T> y = this->matvec(z.view());
  il::Array<T> yout;
  yout.Resize(y.size(), 0.);
  // permut back
  for (il::int_t i = 0; i < nrowpoints; i++) {
    for (int j = 0; j < dof_dimension_; j++) {
      yout[dof_dimension_ * hr_->permutation_0_[i] + j] =
          y[dof_dimension_ * i + j];
    }
  }
  return yout;
}
/* -------------------------------------------------------------------------- */

// H-Matrix vector multiplication with permutation for rectangular matrix
// cases (only 1 permutation) in & out as std::vector todo : write another one
// for the case of 2 permutations (rect. mat cases (for source != receivers)
template <typename T>
std::vector<T> Hmat<T>::matvecOriginal(const std::vector<T> &x) {

  il::Array<T> z{static_cast<il::int_t>(x.size())};
  // permutation of the dofs according to the re-ordering sue to clustering
  il::int_t ncolpoints = this->size(1) / dof_dimension_;
  il::int_t nrowpoints = this->size(0) / dof_dimension_;
  for (il::int_t i = 0; i < ncolpoints; i++) {
    for (int j = 0; j < dof_dimension_; j++) {
      z[dof_dimension_ * i + j] =
          x[dof_dimension_ * hr_->permutation_1_[i] + j];
    }
  }
  il::Array<T> y = this->matvec(z.view());
  std::vector<T> yout;
  yout.assign(y.size(), 0.);
  // permut back
  for (il::int_t i = 0; i < nrowpoints; i++) {
    for (int j = 0; j < dof_dimension_; j++) {
      yout[dof_dimension_ * hr_->permutation_0_[i] + j] =
          y[dof_dimension_ * i + j];
    }
  }
  return yout;
}
/* -------------------------------------------------------------------------- */

// matvect in and outs as std::vector
template <typename T> std::vector<T> Hmat<T>::matvec(const std::vector<T> &x) {
  il::Array<T> xil{static_cast<il::int_t>(x.size())};
  // todo find a better way to convert il::Array to std::vect and vice versa !
  for (long i = 0; i < xil.size(); i++) {
    xil[i] = x[i];
  }
  il::Array<T> yil = this->matvec(xil.view());
  std::vector<T> y;
  y.reserve(static_cast<long>(yil.size()));
  for (long i = 0; i < yil.size(); i++) {
    y.push_back(yil[i]);
  }
  return y;
}
/* -------------------------------------------------------------------------- */

#if defined(BIGWHAM_HDF5)
namespace {
template<typename T> struct HDF5TypeHelper {
  static hid_t type() { throw std::runtime_error("Not implemented"); }
};

template<> struct HDF5TypeHelper<double> {
  static hid_t type() { return  H5T_NATIVE_DOUBLE; }
};

template<> struct HDF5TypeHelper<int64_t> {
  static hid_t type() { return  H5T_NATIVE_LLONG; }
};


template<typename T>
hid_t getHDF5Type(T /*unused*/) {
  return HDF5TypeHelper<T>::type();
}


template <class T> struct is_array_2d : public std::false_type {};
template <class T> struct is_array_2d<il::Array2D<T>> : public std::true_type {};
template <class T> struct is_array_2d<il::Array2DView<T>> : public std::true_type {};
template <class T> struct is_array_2d<il::Array2DEdit<T>> : public std::true_type {};

template <class T>
inline constexpr bool is_array_2d_v = is_array_2d<T>::value;

}
#endif

template <typename T> void Hmat<T>::writeToFile(const std::string &filename) {
#if defined(BIGWHAM_HDF5)
  auto file_id =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  auto hmat_gid =
      H5Gcreate(file_id, "/hmat", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  auto writeAttribute = [](auto &&name, auto &&gid, auto &&value) {
    auto aid = H5Screate(H5S_SCALAR);
    auto attr =
        H5Acreate(gid, name, getHDF5Type(value), aid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, getHDF5Type(value), &value);
    H5Aclose(attr);
    H5Sclose(aid);
  };

  auto writeArray = [](auto &&name, auto &&gid, auto &&array) {
    hsize_t dims[2];

    if constexpr (is_array_2d_v<std::decay_t<decltype(array)>>) {
      dims[0] = hsize_t(array.size(0));
      dims[1] = hsize_t(array.size(1)); // explicit conversion to avoid the compiler warning
    } else {
      dims[0] = hsize_t(array.size());
      dims[1] = 1;
    }
    
    auto datatype_id = getHDF5Type(std::decay_t<decltype(array.data())>{});
    auto dataspace_id = H5Screate_simple(2, dims, NULL);
    auto dataset_id =
        H5Dcreate(gid, std::string(name).c_str(), datatype_id,
                  dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, datatype_id, H5S_ALL, dataspace_id, H5P_DEFAULT,
             array.data());
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
  };

  writeAttribute("dof_dimension", hmat_gid, dof_dimension_);
  writeAttribute("m", hmat_gid, size_[0]);
  writeAttribute("n", hmat_gid, size_[1]);

  auto pattern_gid =
      H5Gcreate(hmat_gid, "pattern", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  writeAttribute("n_FRB", pattern_gid, hr_->pattern_.n_FRB);
  writeAttribute("n_LRB", pattern_gid, hr_->pattern_.n_LRB);
  writeArray("FRB", pattern_gid, hr_->pattern_.FRB_pattern);
  writeArray("LRB", pattern_gid, hr_->pattern_.LRB_pattern);

  H5Gclose(pattern_gid);

  auto permutation_gid =
      H5Gcreate(hmat_gid, "permutations", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  writeArray("rows", permutation_gid, hr_->permutation_0_);
  writeArray("cols", permutation_gid, hr_->permutation_1_);
  H5Gclose(permutation_gid);
    
  auto frb_gid =
      H5Gcreate(hmat_gid, "FRB", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  il::int_t i_frb = 0;
  for (auto &&frb : full_rank_blocks_) {
    writeArray("frb_" + std::to_string(i_frb), frb_gid, *frb);
    ++i_frb;
  }
  H5Gclose(frb_gid);

  auto lrbs_gid =
      H5Gcreate(hmat_gid, "LRB", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  il::int_t i_lrb = 0;
  for (auto &&lrb : low_rank_blocks_) {
    std::string group = "lrb_" + std::to_string(i_lrb);
    auto lrb_gid = H5Gcreate(lrbs_gid, group.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT);
    writeArray("A", lrb_gid, lrb->A);
    writeArray("B", lrb_gid, lrb->B);
    H5Gclose(lrb_gid);
    ++i_lrb;
  }
  H5Gclose(lrbs_gid);
  H5Gclose(hmat_gid);
  H5Fclose(file_id);
#else
  throw std::runtime_error("Recompile with HDF5 support");
#endif
}

template <typename T> void Hmat<T>::readFromFile(const std::string &filename) {
#if defined(BIGWHAM_HDF5)
  auto file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  auto hmat_gid = H5Gopen(file_id, "/hmat", H5P_DEFAULT);

  auto readAttribute = [](auto &&name, auto &&gid, auto &&value) {
    auto attr = H5Aopen(gid, name, H5P_DEFAULT);
    H5Aread(attr, getHDF5Type(value), &value);
    H5Aclose(attr);
  };

  auto readArray = [](auto &&name, auto &&gid, auto &&array) {
    hsize_t dims[2];
    auto datatype_id = getHDF5Type(std::decay_t<decltype(array.data())>{});
    auto dataset_id = H5Dopen(gid, std::string(name).c_str(), H5P_DEFAULT);

    auto dataspace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

    if constexpr (is_array_2d_v<std::decay_t<decltype(array)>>) {
      array.Resize(dims[0], dims[1]);
    } else {
      array.Resize(dims[0]);
    }
    H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array.Data());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
  };

  readAttribute("dof_dimension", hmat_gid, dof_dimension_);
  readAttribute("m", hmat_gid, size_[0]);
  readAttribute("n", hmat_gid, size_[1]);

  auto pattern_gid = H5Gopen(hmat_gid, "pattern", H5P_DEFAULT);
  readAttribute("n_FRB", pattern_gid, hr_->pattern_.n_FRB);
  readAttribute("n_LRB", pattern_gid, hr_->pattern_.n_LRB);

  readArray("FRB", pattern_gid, hr_->pattern_.FRB_pattern);
  readArray("LRB", pattern_gid, hr_->pattern_.LRB_pattern);
  H5Gclose(pattern_gid);

  auto permutation_gid =
      H5Gopen(hmat_gid, "permutations", H5P_DEFAULT);
  readArray("rows", permutation_gid, hr_->permutation_0_);
  readArray("cols", permutation_gid, hr_->permutation_1_);
  H5Gclose(permutation_gid);
  
  std::cout << "Reading HMat from \"" << filename << "\" - " << size_[0] << "x"
            << size_[1] << " - " << dof_dimension_ << std::endl;
  std::cout << " Number of blocks = " << hr_->pattern_.n_FRB + hr_->pattern_.n_LRB
            << std::endl;
  std::cout << " Number of full blocks = " << hr_->pattern_.n_FRB << std::endl;

  auto frbs_gid = H5Gopen(hmat_gid, "FRB", H5P_DEFAULT);
  full_rank_blocks_.resize(hr_->pattern_.n_FRB);
  il::int_t i_frb = 0;
  il::int_t i_frb_10_per = hr_->pattern_.n_FRB * 0.1;
  for (auto &&frb : full_rank_blocks_) {
    if (i_frb % i_frb_10_per == 0)
      std::cout << "." << std::flush;
    frb = std::make_unique<il::Array2D<T>>();
    readArray("frb_" + std::to_string(i_frb), frbs_gid, *frb);
    ++i_frb;
  }
  std::cout << std::endl;
  H5Gclose(frbs_gid);

  std::cout << " Number of low rank blocks = " << hr_->pattern_.n_LRB << std::endl;
  auto lrbs_gid = H5Gopen(hmat_gid, "LRB", H5P_DEFAULT);
  low_rank_blocks_.resize(hr_->pattern_.n_LRB);
  il::int_t i_lrb = 0;
  il::int_t i_lrb_10_per = hr_->pattern_.n_LRB * 0.1;
  for (auto &&lrb : low_rank_blocks_) {
    if (i_lrb % i_lrb_10_per == 0)
      std::cout << "." << std::flush;
    lrb = std::make_unique<LowRank<T>>();
    std::string group = "lrb_" + std::to_string(i_lrb);
    auto lrb_gid = H5Gopen(lrbs_gid, group.c_str(), H5P_DEFAULT);
    readArray("A", lrb_gid, lrb->A);
    readArray("B", lrb_gid, lrb->B);
    H5Gclose(lrb_gid);
    ++i_lrb;
  }
  std::cout << std::endl;
  H5Gclose(lrbs_gid);

  H5Gclose(hmat_gid);
  H5Fclose(file_id);

  isBuilt_FR_ = isBuilt_LR_ = isBuilt_ = true;
#else
  throw std::runtime_error("Recompile with HDF5 support");
#endif
}
/* -------------------------------------------------------------------------- */

template class Hmat<double>;

} // namespace bie
