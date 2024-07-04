#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

#include "il/Array.h"

#include "io/bigwham_io_gen.h"
#include "io/bigwham_io_helper.h"

using namespace bie;
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
  il::ArrayView<T> d{c.data(), c.shape(0)};
  return d;
}
/* -------------------------------------------------------------------------- */

