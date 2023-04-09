#include "hmat/hierarchical_representation.h"
#include "il/Array.h"
#include "io/bigwham_io_gen.h"
#include "io/bigwham_io_helper.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

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
inline il::ArrayView<T> as_array_view(const pbarray<T> &c) {
  il::ArrayView<T> d(c.data(), c.shape(0));
  return d;
}
/* -------------------------------------------------------------------------- */

pbarray<double> PyGetCollocationPoints(const std::vector<double> &coor,
                                       const std::vector<int> &conn,
                                       const std::string &kernel) {

  std::shared_ptr<bie::Mesh> mesh;

  switch (hash_djb2a(kernel)) {
  case "3DT0"_sh: {
    int spatial_dimension = 3;
    int nvertices_per_elt = 3;
    using EltType = bie::Triangle<0>;
    mesh = createMeshFromVect<EltType>(spatial_dimension, nvertices_per_elt,
                                       coor, conn);
    break;
  }
  default: {
    std::cout << "wrong inputs -abort \n";
    il::abort();
  }
  }
  auto v = mesh->collocation_points();

  il::Array<double> pts{v.size(0) * v.size(1), 0.};

  int index = 0;
  for (il::int_t i = 0; i < v.size(0); i++) {
    for (il::int_t j = 0; j < v.size(1); j++) {
      pts[index] = v(i, j);
      index++;
    }
  }

  return as_pyarray<double>(std::move(pts));
}
/* -------------------------------------------------------------------------- */

pbarray<il::int_t> PyGetPermutation(const int dim, const pbarray<double> &pts,
                                    const il::int_t leafsize) {

  il::Array2D<double> v{pts.size() / dim, dim, 0.};
  int index = 0;
  for (il::int_t i = 0; i < v.size(0); i++) {
    for (il::int_t j = 0; j < v.size(1); j++) {
      v(i, j) = pts.data()[index];
      index++;
    }
  }
  auto cluster = bie::cluster(leafsize, il::io, v);
  auto d = cluster.permutation;
  return as_pyarray<il::int_t>(std::move(d));
}
/* -------------------------------------------------------------------------- */

template <typename T>
void declare_array(py::module &m, const std::string &typestr) {
  using Class = il::Array2D<T>;
  std::string pyclass_name = typestr;
  py::class_<Class>(m, pyclass_name.c_str())
      .def(py::init<il::int_t, il::int_t, const T &>())
      .def("shape", &Class::size)
      .def("__setitem__", [](Class &instance, std::array<il::int_t, 2> i,
                             const T &value) { instance(i[0], i[1]) = value; })
      .def("__getitem__",
           [](const Class &instance, std::array<il::int_t, 2> i) {
             return instance(i[0], i[1]);
           })
      .def("__repr__", [](const Class &instance) {
        std::string t1 = "Size : " + std::to_string(instance.size(0)) + " X " +
                         std::to_string(instance.size(1)) + "  \n\n";
        std::string t2 = "";
        for (il::int_t i = 0; i < instance.size(0); ++i) {
          for (il::int_t j = 0; j < instance.size(1); ++j) {
            t2 += std::to_string(instance(i, j)) + " ";
          }
          t2 += "\n";
        }
        return t1 + t2;
      });
}
