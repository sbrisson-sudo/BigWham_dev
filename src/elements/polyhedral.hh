//
// This file is part of BigWham.
//
// Created by Sylvain Brisson on 07/11/2025
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_POLYHEDRAL_H
#define BIGWHAM_POLYHEDRAL_H

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/linearAlgebra/dense/norm.h>
#include "elements/boundary_element.h"


#include <iomanip>

namespace bigwham {

/**
 * @brief General class for polyhedral elements (needed for eigenstrain kernels)
 * 
 * @tparam p 
 */
template <int p> class Polyhedral : public BoundaryElement {

public:
    Polyhedral() : BoundaryElement(3, p) {}
    ~Polyhedral() {}

    /**
     * @brief Set the vertices.
     * 
     * @tparam p 
     * @param xv vertices
     */
    virtual void SetElement(const il::Array2D<double> &coods_vertices) override;

    /**
     * @brief This must returns the face indices, has to be defined in the 
     * inherited methods.
     * 
     */
    virtual void setFaceIndices() = 0;

    /**
     * @brief Has to be overriden but not set for now
     * 
     */
    virtual void SetRotationMatrices() override; 

    // Geometrical relation to point primitives
    bool isPointOnBoundary(const std::array<double, 3>& xyz_obs) const;
    bool isPointInPolyhedron(const std::array<double, 3>& xyz_obs) const;

    // Getters
    double getTol() const { return tol_; };
    int getNumFaces() const { return numFaces_; };
    int getNumVerticesPerFaces() const { return numVerticesPerFace_; };
    il::Array2D<size_t> getFaceIndices() const {return faceIndices_; };
    il::Array2D<double> getFaceNormals() const {return faceNormals_; };
    il::Array2D<double> getFaceCentroids() const {return faceCentroids_; };

protected:
    double tol_;
    int numFaces_;
    int numVerticesPerFace_;
    il::Array2D<size_t> faceIndices_;
    il::Array2D<double> faceNormals_;
    il::Array2D<double> faceCentroids_;
};


template <int p> 
void Polyhedral<p>::SetRotationMatrices() {
    // Runtime exception
    throw std::runtime_error("SetRotationMatrices not defined for Polyhedral");
}


template <int p>
inline void Polyhedral<p>::SetElement(const il::Array2D<double> &xv) {

    IL_EXPECT_FAST(xv.size(1) == spatial_dimension_);
    IL_EXPECT_FAST(xv.size(0) == this->num_vertices_);
    this->vertices_.Resize(num_vertices_, spatial_dimension_);

    // Set the vertices
    for (il::int_t j = 0; j < spatial_dimension_; j++) {
        for (il::int_t i = 0; i < num_vertices_; i++) {
            this->vertices_(i, j) = xv(i, j);
        }
    }

    // Compute the centroid
    for (il::int_t j = 0; j < spatial_dimension_; j++) {
        // always reset centroid when setting the coordinates
        this->centroid_[j] =0;
        for (il::int_t i = 0; i < num_vertices_; i++) {
            this->centroid_[j] = this->centroid_[j] + vertices_(i, j) / num_vertices_;
        }
    }

    // Set collocation point / node to centroid 
    this->collocation_points_(0,0) = this->centroid_[0];
    this->collocation_points_(0,1) = this->centroid_[1];
    this->collocation_points_(0,2) = this->centroid_[2];

    this->nodes_(0,0) = this->centroid_[0];
    this->nodes_(0,1) = this->centroid_[1];
    this->nodes_(0,2) = this->centroid_[2];

    for (il::int_t j = 0; j < spatial_dimension_; j++) {
        this->tangent1_[j] = vertices_(1, j) - vertices_(0, j);
        this->tangent2_[j] = vertices_(num_vertices_ - 1, j) - vertices_(0, j);
    }

    // Set the face indices 
    setFaceIndices();

    IL_EXPECT_FAST(numFaces_ == faceIndices_.size(0));
    IL_EXPECT_FAST(numVerticesPerFace_ == faceIndices_.size(1));

    faceNormals_.Resize(numFaces_, 3);
    faceCentroids_.Resize(numFaces_, 3);

    // Loop on polygonal face
    for (int i_face(0); i_face<numFaces_; i_face++){

        // Get three vertices 
        size_t a_1_i = faceIndices_(i_face,0);
        size_t a_2_i = faceIndices_(i_face,1);
        size_t a_3_i = faceIndices_(i_face,2);

        il::StaticArray<double, 3> a_1 = {il::value, 
            {vertices_(a_1_i, 0), vertices_(a_1_i, 1), vertices_(a_1_i, 2)}
        };
        il::StaticArray<double, 3> a_2 = {il::value, 
            {vertices_(a_2_i, 0), vertices_(a_2_i, 1), vertices_(a_2_i, 2)}
        };
        il::StaticArray<double, 3> a_3 = {il::value, 
            {vertices_(a_3_i, 0), vertices_(a_3_i, 1), vertices_(a_3_i, 2)}
        };

        // Get the edges vectors
        il::StaticArray<double, 3> u_1 = {il::value, 
            {a_2[0] - a_1[0], a_2[1] - a_1[1], a_2[2] - a_1[2]}
        };
        il::StaticArray<double, 3> u_2 = {il::value, 
            {a_3[0] - a_1[0], a_3[1] - a_1[1], a_3[2] - a_1[2]}
        };

        // Compute normal (cross product u_1 x u_2)
        il::StaticArray<double, 3> n = {il::value, {
            u_1[1]*u_2[2] - u_1[2]*u_2[1], 
            u_1[2]*u_2[0] - u_1[0]*u_2[2], 
            u_1[0]*u_2[1] - u_1[1]*u_2[0], 
        }};

        // Normalize it 
        double size_n = il::norm(n, il::Norm::L2);
        for (int i(0); i<3; i++) n[i] /= size_n;

        // Compute the face centroid first (needed for orientation check)
        il::StaticArray<double, 3> face_centroid = {il::value, {0., 0., 0.}};
        for (int i_vert(0); i_vert<numVerticesPerFace_; i_vert++){
            for (int i(0); i<3; i++) 
                face_centroid[i] += vertices_(faceIndices_(i_face, i_vert), i) / numVerticesPerFace_;
        }

        // Vector from polyhedron centroid to face centroid
        il::StaticArray<double, 3> centroid_to_face = {il::value, {
            face_centroid[0] - this->centroid_[0],
            face_centroid[1] - this->centroid_[1],
            face_centroid[2] - this->centroid_[2]
        }};

        // Check if normal points outward by checking dot product
        // If dot product is negative, normal points inward, so flip it
        double dot_product = n[0] * centroid_to_face[0] + 
                            n[1] * centroid_to_face[1] + 
                            n[2] * centroid_to_face[2];

        if (dot_product < 0) {
            // Normal points inward, flip it
            for (int i(0); i<3; i++) n[i] = -n[i];
        }

        // Set the outward-pointing normal
        for (int i(0); i<3; i++) faceNormals_(i_face, i) = n[i];

        // Set the face centroid
        for (int i(0); i<3; i++) faceCentroids_(i_face, i) = face_centroid[i];
    }

    // We set the tolerance used
    // We are not looping on the edges but we dont care 
    double smallest_edge = 1e100;
    for (il::int_t i = 0; i < num_vertices_; i++) {

        double a_1_x = vertices_(i, 0);
        double a_1_y = vertices_(i, 1);
        double a_1_z = vertices_(i, 2);
        double a_2_x = vertices_((i+1)%num_vertices_, 0);
        double a_2_y = vertices_((i+1)%num_vertices_, 1);
        double a_2_z = vertices_((i+1)%num_vertices_, 2);

        double edge_length = std::sqrt( (a_2_x - a_1_x)*(a_2_x - a_1_x) + (a_2_y - a_1_y)*(a_2_y - a_1_y) + (a_2_z - a_1_z)*(a_2_z - a_1_z));
        if (edge_length < smallest_edge) smallest_edge = edge_length;
    } 

    tol_ = smallest_edge * 1e-5;
}

/**
 * @brief Check if a point is on the boundary of the polyhedron
 * 
 * A point is on the boundary if it lies on any face of the polyhedron
 * within the specified tolerance.
 * 
 * @param xyz_obs 3D coordinates of the observation point
 * @return true if point is on boundary (within tolerance)
 */
template <int p>
bool Polyhedral<p>::isPointOnBoundary(const std::array<double, 3>& xyz_obs) const {
    
    // Use a more generous tolerance for boundary detection
    // double boundary_tol = std::max(tol_, 1e-10);
    double boundary_tol = tol_;
    
    // Loop over all faces
    for (int i_face = 0; i_face < numFaces_; i_face++) {
        
        // Get face normal
        il::StaticArray<double, 3> n = {il::value, {
            faceNormals_(i_face, 0),
            faceNormals_(i_face, 1),
            faceNormals_(i_face, 2)
        }};
        
        // Use first vertex of the face as reference point
        size_t v_idx = faceIndices_(i_face, 0);
        il::StaticArray<double, 3> face_point = {il::value, {
            vertices_(v_idx, 0),
            vertices_(v_idx, 1),
            vertices_(v_idx, 2)
        }};
        
        // Vector from face point to observation point
        il::StaticArray<double, 3> v = {il::value, {
            xyz_obs[0] - face_point[0],
            xyz_obs[1] - face_point[1],
            xyz_obs[2] - face_point[2]
        }};
        
        // Distance from point to face plane
        double dist_to_plane = std::abs(v[0] * n[0] + v[1] * n[1] + v[2] * n[2]);
        
        // If point is NOT close to the plane, skip this face
        if (dist_to_plane > boundary_tol) {
            continue;
        }
        
        // Point is close to plane, now check if it's within the face boundary
        // Use a simpler approach: check if point is "inside" all edges when projected
        
        // Project point onto face plane
        double proj_dist = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
        il::StaticArray<double, 3> projected_point = {il::value, {
            xyz_obs[0] - proj_dist * n[0],
            xyz_obs[1] - proj_dist * n[1],
            xyz_obs[2] - proj_dist * n[2]
        }};
        
        // Check if projected point is inside the convex face polygon
        // For a convex polygon, point is inside if it's on the "left" side of all edges
        // (when traversing edges counter-clockwise as seen from the normal direction)
        
        bool inside_face = true;
        for (int i_vert = 0; i_vert < numVerticesPerFace_; i_vert++) {
            int next_vert = (i_vert + 1) % numVerticesPerFace_;
            
            size_t v1_idx = faceIndices_(i_face, i_vert);
            size_t v2_idx = faceIndices_(i_face, next_vert);
            
            il::StaticArray<double, 3> v1 = {il::value, {
                vertices_(v1_idx, 0),
                vertices_(v1_idx, 1),
                vertices_(v1_idx, 2)
            }};
            
            il::StaticArray<double, 3> v2 = {il::value, {
                vertices_(v2_idx, 0),
                vertices_(v2_idx, 1),
                vertices_(v2_idx, 2)
            }};
            
            // Edge vector
            il::StaticArray<double, 3> edge = {il::value, {
                v2[0] - v1[0],
                v2[1] - v1[1],
                v2[2] - v1[2]
            }};
            
            // Vector from edge start to point
            il::StaticArray<double, 3> to_point = {il::value, {
                projected_point[0] - v1[0],
                projected_point[1] - v1[1],
                projected_point[2] - v1[2]
            }};
            
            // Cross product: edge × to_point
            il::StaticArray<double, 3> cross = {il::value, {
                edge[1] * to_point[2] - edge[2] * to_point[1],
                edge[2] * to_point[0] - edge[0] * to_point[2],
                edge[0] * to_point[1] - edge[1] * to_point[0]
            }};
            
            // Dot product with face normal
            double dot = cross[0] * n[0] + cross[1] * n[1] + cross[2] * n[2];
            
            // If dot product is negative, point is on the "outside" of this edge
            // Use a small tolerance for numerical stability
            if (dot < -boundary_tol) {
                inside_face = false;
                break;
            }
        }
        
        if (inside_face) {
            return true;  // Point is on this face
        }
    }
    
    return false;
}


/**
 * @brief Check if a point is inside the convex polyhedron
 * 
 * For a convex polyhedron, a point is inside if it is on the "inside" side
 * of all faces (i.e., negative side of all face normal vectors pointing outward).
 * 
 * @param xyz_obs 3D coordinates of the observation point
 * @return true if point is inside the polyhedron
 */
template <int p>
bool Polyhedral<p>::isPointInPolyhedron(const std::array<double, 3>& xyz_obs) const {
    
    // For a convex polyhedron, check if point is on the negative side of all faces
    // Assumes face normals point outward
    
    for (int i_face = 0; i_face < numFaces_; i_face++) {
        
        // Get face normal
        il::StaticArray<double, 3> n = {il::value, {
            faceNormals_(i_face, 0),
            faceNormals_(i_face, 1),
            faceNormals_(i_face, 2)
        }};
        
        // Use first vertex of the face as reference point
        size_t v_idx = faceIndices_(i_face, 0);
        il::StaticArray<double, 3> face_point = {il::value, {
            vertices_(v_idx, 0),
            vertices_(v_idx, 1),
            vertices_(v_idx, 2)
        }};
        
        // Vector from face point to observation point
        il::StaticArray<double, 3> v = {il::value, {
            xyz_obs[0] - face_point[0],
            xyz_obs[1] - face_point[1],
            xyz_obs[2] - face_point[2]
        }};
        
        // Compute signed distance (dot product with normal)
        double signed_dist = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
        
        // If point is on the positive side (outside) of any face, it's not inside
        // Use tolerance to handle points very close to boundary


        // if (signed_dist > tol_) {
        if (signed_dist > 0) {
            return false;
        }
    }
    
    return true;
}


} // namespace bigwham


#endif // BIGWHAM_POLYHEDRAL_H

