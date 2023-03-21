#include "core/be_mesh.h"
#include "elements/triangle.h"



namespace bie {
/* -------------------------------------------------------------------------- */

template <class ElemType>
il::Array2D<double> BEMesh<ElemType>::get_vertices(il::int_t element_id) const {
  il::Array2D<double> vertices{this->num_vertices_, this->spatial_dimension_};
  // loop over the vertices
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    for (il::int_t i = 0; i < num_vertices_; i++) {
      vertices(i, j) = coordinates_(connectivity_(element_id, i), j);
    }
  }
  return vertices;
}
/* -------------------------------------------------------------------------- */

template <class ElemType> void BEMesh<ElemType>::construct_mesh() {
  il::int_t j = 0;
  for (il::int_t i = 0; i < element_list_.size(); i++) {
    // define element
    std::shared_ptr<BoundaryElement> elem = std::make_shared<ElemType>();
    auto vertices = this->get_vertices(i);
    elem->set_element(vertices);
    this->element_list_[i] = std::move(elem);

    // define collocation points
    auto elem_colloc_pts = elem->get_collocation_points();
    for (il::int_t k = 0; k < spatial_dimension_; k++) {
      for (il::int_t j1 = 0; j1 < num_colloc_pts_per_element_; j1++) {
        collocation_points_(j + j1, k) = elem_colloc_pts(j1, k);
      }
    }
    j = j + num_colloc_pts_per_element_;
  }
}
/* -------------------------------------------------------------------------- */

template class BEMesh<Triangle<0>>;
template class BEMesh<Triangle<1>>;
template class BEMesh<Triangle<2>>;

} // namespace bie
