#include "il/container/1d/Array.h"
#include "io/bigwham_io_gen.h"
#include "io/bigwham_io_helper.h"
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array PyGetCollocationPoints(const std::vector<double> &coor,
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

  return py::array(pts.size(), pts.data());
}
