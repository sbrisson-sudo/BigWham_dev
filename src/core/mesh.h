#ifndef BIGWHAM_BEMESH_H
#define BIGWHAM_BEMESH_H

#include <memory>
#include <vector>

#include "core/elements/boundary_element.h"

namespace bie {

class Mesh {
protected:
  il::int_t spatial_dimension_;
  il::int_t num_nodes_;
  il::int_t num_elements_;
  std::vector<std::unique_ptr<BoundaryElement>> element_list_;
  // Coordinates of the nodes - size: number of nodes x problem dimension
  il::Array2D<double> coordinates_;
  // Connectivity matrix - size: number of elements x   n_vertex per element
  il::Array2D<il::int_t> connectivity_;

public:
  Mesh(){};
  Mesh(const il::Array2D<double> &coordinates,
       const il::Array2D<il::int_t> &connectivity) {
    num_nodes_ = coordinates.size(0);
    spatial_dimension_ = coordinates.size(1);
    num_elements_ = connectivity.size(0);

    coordinates_ = coordinates;
    connectivity_ = connectivity;
  };
  ~Mesh(){};
  virtual std::unique_ptr<BoundaryElement>
  getElement(il::int_t node_id) const = 0;
  il::int_t get_num_elements() const { return num_elements_; };
  // nodal coordinates related.
  il::Array2D<double> get_coordinates() const { return coordinates_; };
  // Read a particular element of the coordinates coordinates
  double get_coordinates(il::int_t k, il::int_t i) {
    return coordinates_(k, i);
  };
  // Read a particular element of the coordinates coordinates
  il::StaticArray<double, 2> coordinates(il::int_t k) const {
    il::StaticArray<double, 2> temp;
    temp[0] = coordinates_(k, 0);
    temp[1] = coordinates_(k, 1);
    return temp;
  };
  // connectivity related
  il::Array2D<il::int_t> get_connectivity() const { return connectivity_; };

  // method to return all collocation points of the mesh.
  il::Array2D<double> get_collocation_points() {
    il::Array2D<double> xv{
        number_vertex_,
        spatial_dimension_,
        0,
    };
    il::Array2D<double> Xcol{nodes_per_element_, spatial_dimension_, 0.},
        colPoints{n_elts_ * nodes_per_element_, spatial_dimension_, 0};
    il::int_t j = 0;
    std::cout << "Number of elements : " << n_elts_
              << ", Nodes_per_element : " << nodes_per_element_
              << ", Spatial dim : " << spatial_dimension_ << "\n";
    std::cout << "Number of vertex per element : " << number_vertex_ << "\n";
    for (il::int_t e = 0; e < n_elts_; e++) {
      this->setCurrentElement(e);
      Xcol = this->element_def_.getCollocationPoints();
      auto centro = this->element_def_.getCentroid();
      //  std::cout << "centro  " << e << "  centroid " << centro[0]  <<" -" <<
      //  centro[1]   <<"\n";
      // std::cout << "xcol  " << e << "size " << Xcol(0,0) <<"-" << Xcol(0,1)
      // <<"\n";
      for (il::int_t k = 0; k < spatial_dimension_; k++) {
        for (il::int_t j1 = 0; j1 < nodes_per_element_; j1++) {
          colPoints(j + j1, k) = Xcol(j1, k);
        }
      }
      j = j + nodes_per_element_;
    };
    return colPoints;
  };

  il::int_t inElement(il::int_t point_n) const {
    return il::floor(point_n / (nodes_per_element_)); // element
  }

  il::int_t localCollocationPointId(il::int_t point_n) const {
    return (point_n % (nodes_per_element_));
  }
};

} // namespace bie

#endif // BIGWHAM_BEMESH_H
