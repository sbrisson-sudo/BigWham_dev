#include <cmath>
#include <numbers>
#include <iostream>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "elastic_2dT0S0_V_element.hh"

namespace bigwham {

using tensor2d = std::array<std::array<double, 2>, 2>;

double point_to_segment_distance(
    const std::array<double, 2>& P,
    const std::array<double, 2>& A,
    const std::array<double, 2>& B)
{
    double dx = B[0] - A[0];
    double dy = B[1] - A[1];
    
    if (dx == 0.0 && dy == 0.0) {
        // A and B are the same point
        dx = P[0] - A[0];
        dy = P[1] - A[1];
        return std::sqrt(dx*dx + dy*dy);
    }

    // Project point P onto the segment AB, computing parameter t
    double t = ((P[0] - A[0]) * dx + (P[1] - A[1]) * dy) / (dx*dx + dy*dy);
    t = std::max(0.0, std::min(1.0, t)); // clamp t to [0,1]

    double proj_x = A[0] + t * dx;
    double proj_y = A[1] + t * dy;

    double dist_x = P[0] - proj_x;
    double dist_y = P[1] - proj_y;
    
    return std::sqrt(dist_x*dist_x + dist_y*dist_y);
}

bool is_point_on_triangle_boundary(
    const std::array<std::array<double, 2>, 3>& tri_vertices,
    const std::array<double, 2>& xy_obs,
    double tol)
{
    for (int i = 0; i < 3; ++i) {
        const auto& A = tri_vertices[i];
        const auto& B = tri_vertices[(i + 1) % 3];
        if (point_to_segment_distance(xy_obs, A, B) <= tol) {
            return true;
        }
    }
    return false;
}

bool is_point_in_triangle(
    const std::array<std::array<double, 2>, 3>& tri_vertices,
    const std::array<double, 2>& P,
    double tol = 1e-12)
{
    const auto& A = tri_vertices[0];
    const auto& B = tri_vertices[1];
    const auto& C = tri_vertices[2];

    // Compute vectors
    double v0x = C[0] - A[0];
    double v0y = C[1] - A[1];
    double v1x = B[0] - A[0];
    double v1y = B[1] - A[1];
    double v2x = P[0] - A[0];
    double v2y = P[1] - A[1];

    // Compute dot products
    double dot00 = v0x*v0x + v0y*v0y;
    double dot01 = v0x*v1x + v0y*v1y;
    double dot02 = v0x*v2x + v0y*v2y;
    double dot11 = v1x*v1x + v1y*v1y;
    double dot12 = v1x*v2x + v1y*v2y;

    // Compute barycentric coordinates
    double denom = dot00 * dot11 - dot01 * dot01;
    if (std::abs(denom) < tol) return false; // Degenerate triangle
    double invDenom = 1.0 / denom;
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle (allowing small tolerance)
    return (u >= -tol) && (v >= -tol) && (u + v <= 1.0 + tol);
}


il::StaticArray<double, 2> V_twoD_triangle_0(
    std::array<std::array<double, 2>, 3>  tri_vertices,
    il::StaticArray<double, 2> xy_obs,
    il::StaticArray<double, 2> n_obs,
    double G, double nu,
    bool recursive_call
){

    // Physical parameters    
    double lambda = 2*G*nu / (1 - 2*nu);
    double factor = 1. / (2*std::numbers::pi) / (1-nu);

    // Return value 
    il::StaticArray<double, 2> t_i;

    // Determine if xy_obs is on the boundary of the triangle
    double tol = 1e-4;
    bool is_xy_on_boundary = is_point_on_triangle_boundary(tri_vertices, {xy_obs[0],xy_obs[1]}, tol);

    // If on the boundary : 
    if (is_xy_on_boundary){
        // Evaluate the traction at two points on each side and take the average

        // std::cerr << "x_obs = [" << xy_obs[0] << ", " << xy_obs[1] << "]\n";
        // std::cerr << "n_obs = [" << n_obs[0] << ", " << n_obs[1] << "]\n";

        if (recursive_call) {
            throw std::runtime_error("Error: the point is still on the face when moving along the provided normal.");
        }

        double eps = 2*tol; // This shouldn't be hard coded

        il::StaticArray<double, 2> xy_obs_p{il::value, {
            xy_obs[0] + eps*n_obs[0], 
            xy_obs[1] + eps*n_obs[1]}
        };
        il::StaticArray<double, 2> xy_obs_m{il::value, {
            xy_obs[0] - eps*n_obs[0], 
            xy_obs[1] - eps*n_obs[1]}
        };

        auto t_p = V_twoD_triangle_0(tri_vertices, xy_obs_p, n_obs, G, nu, true);
        auto t_m = V_twoD_triangle_0(tri_vertices, xy_obs_m, n_obs, G, nu, true);
        
        t_i[0] = (t_p[0] + t_m[0])/2;
        t_i[1] = (t_p[1] + t_m[1])/2;
    }
    // If not on the boundary :
    else {
        // Compute the phi_ij 
        tensor2d phi_ij = phi_ij_triangle(
            tri_vertices, xy_obs
        );

        // Compute strain
        tensor2d eps_ij;
        for (int i(0); i<2; i++){
            for (int j(0); j<2; j++){
                eps_ij[i][j] = phi_ij[i][j] * factor;
            }
        }

        bool is_xy_in_inclusion = is_point_in_triangle(tri_vertices, {xy_obs[0],xy_obs[1]});
        if (is_xy_in_inclusion){
            eps_ij[0][0] -= 1.0;
            eps_ij[1][1] -= 1.0;   
        }

        double eps_kk = eps_ij[0][0] + eps_ij[1][1];

        // Compute stress
        tensor2d sigma_ij;
        sigma_ij[0][0] = 2*G*eps_ij[0][0] + lambda*eps_kk;
        sigma_ij[1][0] = 2*G*eps_ij[1][0];
        sigma_ij[1][1] = 2*G*eps_ij[1][1] + lambda*eps_kk;
        sigma_ij[0][1] = 2*G*eps_ij[0][1]; // Same as 1,0

        // Project stress
        t_i[0] = sigma_ij[0][0]*n_obs[0] + sigma_ij[1][0]*n_obs[1];
        t_i[1] = sigma_ij[1][0]*n_obs[0] + sigma_ij[1][1]*n_obs[1];
    }

    return t_i;
}

std::array<std::array<double, 2>, 2> phi_ij_triangle(
    std::array<std::array<double, 2>, 3>  tri_vertices,
    il::StaticArray<double, 2> xy_obs
){

    std::array<std::array<double, 2>, 2> phi_ij = {};
    
    // Loop on edges 
    for (int i_edge(0); i_edge<3; i_edge++){

        // Compute normal vector 
        il::StaticArray<double, 2> A_I{il::value, {
            tri_vertices[i_edge][0], 
            tri_vertices[i_edge][1]}
        };
        il::StaticArray<double, 2> A_Ip1{il::value, {
            tri_vertices[(i_edge+1)%3][0], 
            tri_vertices[(i_edge+1)%3][1]}
        };

        il::StaticArray<double, 2> normal{il::value, {
            -(A_Ip1[1] - A_I[1]), 
             A_Ip1[0] - A_I[0]}
        };

        double normal_norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
        // normal[0] /= normal_norm;
        // normal[1] /= normal_norm;
        normal[0] = -normal[0] / normal_norm;
        normal[1] = -normal[1] / normal_norm;

        // Compute line integral 
        il::StaticArray<double, 2> Igrad_j{0.};
        for (int j(0); j<2; j++){
            auto [val, status] = IGrad_j_lineintegral(A_I, A_Ip1, xy_obs, j);
            Igrad_j[j] = val;

            if (status != GSL_SUCCESS){
                std::cerr << "Error with QAGS :\n";
                std::cerr << "tri verts = [[" << tri_vertices[0][0] << ", " << tri_vertices[0][1] << "],[" << tri_vertices[1][0] << ", " << tri_vertices[1][1] << "],[" << tri_vertices[2][0] << ", " << tri_vertices[2][1] << "]]\n";
                std::cerr << "x_obs = [" << xy_obs[0] << ", " << xy_obs[1] << "]\n";

                throw std::runtime_error("Error: QAGS failed, the integrand is likely singular.");
            } 
        }

        // Add outer product 
        phi_ij[0][0] += normal[0] * Igrad_j[0];
        phi_ij[1][0] += normal[1] * Igrad_j[0];
        phi_ij[0][1] += normal[0] * Igrad_j[1];
        phi_ij[1][1] += normal[1] * Igrad_j[1];
    }

    return phi_ij;
}

// parameters to pass to the integrand
struct IntegrandI1jParams {
    il::StaticArray<double, 2> b_i;
    il::StaticArray<double, 2> l_i;
    unsigned int j;
};

double integrandI1j(double xi, void* p) {
    // retrieve parameters
    auto* params = static_cast<IntegrandI1jParams*>(p);
    const auto& b_i = params->b_i;
    const auto& l_i = params->l_i;
    unsigned int j = params->j;

    double num = b_i[j] + 0.5 * l_i[j] * xi;

    double denom = 0.0;
    for (size_t k = 0; k < 2; ++k) {
        double tmp = b_i[k] + 0.5 * l_i[k] * xi;
        denom += tmp * tmp;
    }

    return num / denom;
}


std::pair<double,int> IGrad_j_lineintegral(
    const il::StaticArray<double, 2> A_I,
    const il::StaticArray<double, 2> A_Ip1,
    const il::StaticArray<double, 2> xy_obs,
    unsigned int j
){

    // Helper vectors 
    il::StaticArray<double, 2> a_1 {il::value, {
        A_I[0] - xy_obs[0], 
        A_I[1] - xy_obs[1]
    }};
    il::StaticArray<double, 2> a_2 {il::value, {
        A_Ip1[0] - xy_obs[0], 
        A_Ip1[1] - xy_obs[1]
    }};
    il::StaticArray<double, 2> l_i {il::value, {
        a_2[0] - a_1[0], 
        a_2[1] - a_1[1]
    }};
    il::StaticArray<double, 2> b_i {il::value, {
        (a_2[0] + a_1[0])/2, 
        (a_2[1] + a_1[1])/2
    }};

    // pack it into params struct
    IntegrandI1jParams params{b_i, l_i, j};

    // Calls quadmath
    gsl_function F;
    F.function = &integrandI1j;
    F.params = &params;

    gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
    double result, error;

    int status = gsl_integration_qags(&F, -1, 1, 1.49e-8, 1.49e-8, 1000, w, &result, &error);

    gsl_integration_workspace_free(w);

    // mult by jacobian 
    double l_i_norm = std::sqrt(l_i[0]*l_i[0] + l_i[1]*l_i[1]);
    result *= l_i_norm / 2.0;
    return std::make_pair(result, status);
}

} // namespace bigwham