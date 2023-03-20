#include "core/be_mesh.h"

namespace bie {

  template<class ElemType> void BEMesh<ElemType>::constuct_elements() {
  }

  il::Array2D<double> getVertices(il::int_t ne) const {
    il::Array2D<double> vertElt{number_vertex_, spatial_dimension_};
    // loop over the vertices
    for (il::int_t j = 0; j < spatial_dimension_; j++) {
      for (il::int_t i = 0; i < number_vertex_; i++) {
        vertElt(i, j) = coordinates_(connectivity_(ne, i), j);
      }
    }
    return vertElt;
  }

  void setCurrentElement(il::int_t ne) {
    il::Array2D<double> xv{
        number_vertex_,
        spatial_dimension_,
        0,
    };
    xv = this->getVertices(ne);
    this->element_def_.setElement(xv);
  }

  // method to return all collocation points of the mesh.
  il::Array2D<double> getCollocationPoints() {
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


/* -------------------------------------------------------------------------- */
  il::Array2D<double> getVertices(il::int_t ne) const {
    il::Array2D<double> vertElt{number_vertex_, spatial_dimension_};
    // loop over the vertices
    for (il::int_t j = 0; j < spatial_dimension_; j++) {
      for (il::int_t i = 0; i < number_vertex_; i++) {
        vertElt(i, j) = coordinates_(connectivity_(ne, i), j);
      }
    }
    return vertElt;
  }

  void setCurrentElement(il::int_t ne) {
    il::Array2D<double> xv{
        number_vertex_,
        spatial_dimension_,
        0,
    };
    xv = this->getVertices(ne);
    this->element_def_.setElement(xv);
  }

  // method to return all collocation points of the mesh.
  il::Array2D<double> getCollocationPoints() {
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
  void construct_elements();

} // namespace bie
