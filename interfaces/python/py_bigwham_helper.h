#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

#include "il/Array.h"

#include "io/bigwham_io.h"
#include "io/bigwham_io_helper.h"

using namespace bigwham;
namespace py = pybind11;
template <typename T>
using pbarray = py::array_t<T, py::array::c_style | py::array::forcecast>;

/* -------------------------------------------------------------------------- */

template <typename T> inline py::array_t<T> as_pyarray(il::Array<T> &&seq) {
  auto size = seq.size();
  auto data = seq.data();
  std::unique_ptr<il::Array<T>> seq_ptr =
      std::make_unique<il::Array<T>>(std::move(seq));
  auto capsule = py::capsule(seq_ptr.get(), [](void *p) {
    std::unique_ptr<il::Array<T>>(reinterpret_cast<il::Array<T> *>(p));
  });
  seq_ptr.release();
  return pbarray<T>(size, data, capsule);
}
/* -------------------------------------------------------------------------- */
template <typename T>
inline il::ArrayEdit<T> as_array_edit(pbarray<T> &c) {
    py::buffer_info buf_info = c.request();
    double *ptr = static_cast<T *>(buf_info.ptr);

    int X = buf_info.shape[0];
    il::ArrayEdit<T>  d{ptr, X};//d{&c[0], static_cast<int>(c.size())};
    return d;
}
/* -------------------------------------------------------------------------- */

template <typename T>
inline il::ArrayView<T> as_array_view(const pbarray<T> &c) {
  il::ArrayView<T> d{c.data(), c.shape(0)};
  return d;
}
/* -------------------------------------------------------------------------- */

pbarray<double> PyGetCollocationPoints(const std::vector<double> &coor,
                                       const std::vector<int> &conn,
                                       const std::string &mesh_type) {

  std::shared_ptr<bigwham::Mesh> mesh;

  switch (hash_djb2a(mesh_type)) {
  case "3DT0"_sh: {
    int spatial_dimension = 3;
    int nvertices_per_elt = 3;
    using EltType = bigwham::Triangle<0>;
    mesh = CreateMeshFromVect<EltType>(spatial_dimension, nvertices_per_elt,
                                       coor, conn);
    break;
  }
  default: {
    std::cout << "Wrong Mesh Type, Not Implemented yet!!\n";
    il::abort();
  }
  }
  auto v = mesh->collocation_points();

  il::Array<double> pts{v.size(0) * v.size(1), 0.};

  int index = 0;
    for (il::int_t j = 0; j < v.size(1); j++) {
        for (il::int_t i = 0; i < v.size(0); i++) {
      pts[index] = v(i, j);
      index++;
    }
  }

  return as_pyarray<double>(std::move(pts));
}
/* -------------------------------------------------------------------------- */

pbarray<il::int_t> PyGetPermutation(const int dim, const pbarray<double> &Collocation_pts,
                                    const il::int_t leafsize) {

  il::Array2D<double> v{Collocation_pts.size() / dim, dim, 0.};
  int index = 0;
    for (il::int_t j = 0; j < v.size(1); j++) {
        for (il::int_t i = 0; i < v.size(0); i++) {
      v(i, j) = Collocation_pts.data()[index];
      index++;
    }
  }
  auto cluster = bigwham::cluster(leafsize, il::io, v);
  auto d = cluster.permutation;
  return as_pyarray<il::int_t>(std::move(d));
}