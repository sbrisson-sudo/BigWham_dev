//
// This file is part of HFP.
//
// Created by Brice Lecampion on 12.12.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIE_DOMAINMESH_H
#define BIE_DOMAINMESH_H

#include <limits>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/container/1d/StaticArray.h>
#include <il/linearAlgebra.h>

// todo : most of this is to be removed.
// this is possibly redundant and should be removed alltogether.
// if we interface with just a list of points to get

namespace bie {

class DomainMesh {
 private:
  il::Array2D<double> nodes_;

  il::Array2D<il::int_t> connectivity_;

  il::Array<il::int_t>
      matid_;  // usually in this kind of domain wellMesh: material
  // is different for each element.. never known

  il::int_t nodes_per_elt_;

 public:
  ///////////////////////////////////////////////////
  // constructor ....
  ///////////////////////////////////////////////////
  DomainMesh(){};

  DomainMesh(const il::Array2D<double> &nodes,
             const il::Array2D<il::int_t> &connectivity,
             const il::Array<il::int_t> &matid) {
    IL_EXPECT_FAST(connectivity.size(0) == matid.size());

    nodes_ = nodes;
    connectivity_ = connectivity;
    matid_ = matid;

    nodes_per_elt_ = connectivity.size(1);
  };

  ///////////////////////////////////////////////////
  // get functions
  ///////////////////////////////////////////////////

  // returning the matid corresponding to an element.
  il::int_t getmatid(il::int_t k) const { return matid_[k]; };

  il::int_t numberOfElts() const { return connectivity_.size(0); };

  ///////////////////////////////////////////////////
  // Methods
  ///////////////////////////////////////////////////

  il::StaticArray<double, 2> elementCentroid(il::int_t k) {
    // compute the xy coor of  element k centroid.
    il::StaticArray<double,2> centroid{0.};

    for (il::int_t i = 0; i < nodes_per_elt_; i++) {
      centroid[0] += nodes_(connectivity_(k, i), 0) / nodes_per_elt_;
      centroid[1] += nodes_(connectivity_(k, i), 1) / nodes_per_elt_;
    }
    return centroid;
  }

  il::Array2D<double> allElementCentroid() {
    // compute the xy coor of  element k centroid.
    il::Array2D<double> allcentroids{connectivity_.size(0), 2, 0.};
    for (il::int_t k = 0; k < connectivity_.size(0); k++) {
      allcentroids(k, 0) = 0.;
      allcentroids(k, 1) = 0.;
      for (il::int_t i = 0; i < nodes_per_elt_; i++) {
        allcentroids(k, 0) += nodes_(connectivity_(k, i), 0) / nodes_per_elt_;
        allcentroids(k, 1) += nodes_(connectivity_(k, i), 1) / nodes_per_elt_;
      }
    }
    return allcentroids;
  }

  // method locating the element in which the xy point is belonging too
  il::int_t locate(il::StaticArray<double, 2> &xy) {

//    il::Array2D<double> centroids = allElementCentroid();

    // compute distance from all centroids.
//    double dist;
//    double min = std::numeric_limits<double>::max();
    il::int_t min_pos = -1;

    for (il::int_t i = 0; i < connectivity_.size(0); i++) {
      double angle_a, sum_angles = 0.,
              x_a, y_a, x_b, y_b,
              size_s, d, sign_d;
      il::StaticArray<double, 2> s, n, p;
      il::Array<double> p_angles {connectivity_.size(1)};
      x_a = nodes_(connectivity_(i, connectivity_.size(1) - 1), 0);
      y_a = nodes_(connectivity_(i, connectivity_.size(1) - 1), 1);
      for (il::int_t j = 0; j < connectivity_.size(1); j++) {
          x_b = nodes_(connectivity_(i, j), 0);
          y_b = nodes_(connectivity_(i, j), 1);
          s[0] = x_b - x_a; // tangent vector
          s[1] = y_b - y_a;
          size_s = std::sqrt(s[0] * s[0] + s[1] * s[1]);
          s[0] /= size_s; // normalize it
          s[1] /= size_s;
          n[0] = -1. * s[1]; // normal, left to s
          n[1] = s[0];
          p[0] = (0.5 * x_b + 0.5 * x_a) - xy[0]; // edge middle - xy
          p[1] = (0.5 * y_b + 0.5 * y_a) - xy[1];
          d = il::dot(p, n);
          if (d == 0.) { // if xy is on the border
              min_pos = i; // take the cell closer to the top of the list
              return min_pos;
              // break;
          } else {
              sign_d = (d > 0.) ? 1. : -1.;
              p_angles[j] = std::atan2(n[1] * sign_d, n[0] * sign_d);
              x_a = x_b;
              y_a = y_b;
          }
      }
      angle_a = p_angles[connectivity_.size(1) - 1];
      for (il::int_t j = 0; j < connectivity_.size(1); j++) {
          sum_angles += p_angles[j] - angle_a;
          if (p_angles[j] - angle_a > il::pi) {
              sum_angles -= 2.0 * il::pi;
          } else if (p_angles[j] - angle_a <= -il::pi) {
              sum_angles += 2.0 * il::pi;
          }
          angle_a = p_angles[j];
      }
      if (std::fabs(sum_angles) >= il::pi) {
//      if (std::fabs(std::fabs(sum_angles) - 2.0 * il::pi) < 0.01) {
        min_pos = i;
        break;
      }
//      const double dx = centroids(i, 0) - xy[0];
//      const double dy = centroids(i, 1) - xy[1];
//      dist = dx * dx + dy * dy;
//      if (dist < min) {
//        min = dist;
//        min_pos = i;
//      }
    }

    // this will work for a regular Quad Mesh but NOT necessarily for
    // unstructured triangular Mesh
    return min_pos;  // let's be bold don t do any more checks  ! aie aie.
    // todo bulletproof this locate routine -> find all elements around via node
    // sharing, then check it is not in one of the element around.
    // this will return the first element found if xy is exactly a node of the
    // background wellMesh.
  }
};

// if Cartesian quad only : we could be more optimized.
//  il::Array<il::int_t>  LocateInQuad(il::Array2D<double> &xy){
//
//  IL_EXPECT_FAST(nodes_per_elt_==4); // we want quad here
//
//
//
//  };
};

#endif  // BIE_DOMAINMESH_H
