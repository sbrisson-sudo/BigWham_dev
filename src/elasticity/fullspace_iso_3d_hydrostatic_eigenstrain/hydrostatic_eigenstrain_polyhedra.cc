#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numbers>

#include "hydrostatic_eigenstrain_polyhedra.hh"

// #define DEBUG_KUVSHINOV

namespace bigwham {

using tensor3d = std::array<std::array<double, 3>, 3>;

template <typename T, long int p>
void prettyPrintArray(const il::StaticArray<T, p> &v, std::ostream& stream=std::cout) {

    // Determine the maximum width needed for any element in the matrix
    int maxWidth = 0;
    for (il::int_t i = 0; i < p; ++i) {
        std::ostringstream oss;
        oss << v[i];
        int width = oss.str().length();
        if (width > maxWidth) {
            maxWidth = width;
        }
    }

    // Print the vector with proper formatting
    stream << "[";
    for (il::int_t i = 0; i < p; ++i) 
        stream << std::setw(maxWidth + 2) << v[i];
    stream << "]";
}

il::StaticArray<double, 3> V_threeD_polyhedra_0(
    const Polyhedral<0> &polyhedra,
    const il::StaticArray<double, 3> x_obs,
    const il::StaticArray<double, 3> n_obs,
    double G, double nu,
    bool recursive_call
){

    // Returned traction 
    il::StaticArray<double, 3> t_i;

    double lambd = 2*G*nu / (1 - 2*nu);
    double factor = -1./(4*std::numbers::pi) * (1+nu)/(1-nu);

    bool is_x_on_boundary = polyhedra.isPointOnBoundary({x_obs[0],x_obs[1],x_obs[2]});

    #ifdef DEBUG_KUVSHINOV
    if (is_x_on_boundary){
        std::cerr << "On the boundary : ";
        prettyPrintArray(x_obs, std::cerr);
        std::cerr << "\n";
    }
    #endif

    // If on the boundary : 
    if (is_x_on_boundary){

        if (recursive_call) {
            throw std::runtime_error("Error: the point is still on the face when moving along the provided normal.");
        }

        double eps = 1e3 * polyhedra.getTol();

        il::StaticArray<double, 3> x_obs_p{il::value, {
            x_obs[0] + eps*n_obs[0], 
            x_obs[1] + eps*n_obs[1],
            x_obs[2] + eps*n_obs[2]
        }};
        il::StaticArray<double, 3> x_obs_m{il::value, {
            x_obs[0] - eps*n_obs[0], 
            x_obs[1] - eps*n_obs[1],
            x_obs[2] - eps*n_obs[2]
        }};

        // Ensuring that not both are considered as within the inclusion (tolerance issue)
        if (polyhedra.isPointOnBoundary({x_obs_p[0],x_obs_p[1],x_obs_p[2]}) && 
            polyhedra.isPointOnBoundary({x_obs_m[0],x_obs_m[1],x_obs_m[2]}))
            throw std::runtime_error("In V_threeD_polyhedra_0 : both points are considered as within the inclusion.");

        auto t_p = V_threeD_polyhedra_0(polyhedra, x_obs_p, n_obs, G, nu, true);
        auto t_m = V_threeD_polyhedra_0(polyhedra, x_obs_m, n_obs, G, nu, true);

        for (int i(0); i<3; i++) t_i[i] = (t_p[i] + t_m[i])/2;
    }
    // Not on the polyhedra boundary 
    else {

        // Compute phi_ij
        tensor3d phi_ij = phi_ij_polyhedra(polyhedra, x_obs);

        // Compute strain
        tensor3d eps_ij{0.};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                eps_ij[i][j] = phi_ij[i][j] * factor;

        // Correction if within the inclusion
        bool is_x_in_inclusion = polyhedra.isPointInPolyhedron({x_obs[0],x_obs[1],x_obs[2]});

        #ifdef DEBUG_KUVSHINOV
        if (is_x_in_inclusion){
            std::cerr << "Within the inclusion : ";
            prettyPrintArray(x_obs, std::cerr);
            std::cerr << "\n";
        }
        #endif

        if (is_x_in_inclusion)
            for (int i = 0; i < 3; ++i)
                eps_ij[i][i] -= 1.0;

        // Compute stresses
        tensor3d sigma_ij{0.};
        double eps_kk = eps_ij[0][0] + eps_ij[1][1] + eps_ij[2][2];
        for (int i = 0; i < 3; ++i){
            sigma_ij[i][i] += lambd * eps_kk;
            for (int j = 0; j < 3; ++j)
                sigma_ij[i][j] += 2*G*eps_ij[i][j];
        }

        // Get traction 
        t_i[0] = sigma_ij[0][0]*n_obs[0] + sigma_ij[0][1]*n_obs[1] + sigma_ij[0][2]*n_obs[2];
        t_i[1] = sigma_ij[1][0]*n_obs[0] + sigma_ij[1][1]*n_obs[1] + sigma_ij[1][2]*n_obs[2];
        t_i[2] = sigma_ij[2][0]*n_obs[0] + sigma_ij[2][1]*n_obs[1] + sigma_ij[2][2]*n_obs[2];
    }

    return t_i;
}


std::array<std::array<double, 3>, 3>  phi_ij_polyhedra(
    const Polyhedral<0> &polyhedra,
    const il::StaticArray<double, 3> xy_obs
){

    std::array<std::array<double, 3>, 3> phi_ij = {};

    int numFaces = polyhedra.getNumFaces();
    int numVertPerFace = polyhedra.getNumVerticesPerFaces();
    auto faceNormals = polyhedra.getFaceNormals();
    auto vertices = polyhedra.vertices();
    auto faceIndices = polyhedra.getFaceIndices();

    // Tolerance: relative to element size (cbrt of volume)
    double eps = polyhedra.getTol();

    // Loop on polygonal faces
    for (int i_face(0); i_face<numFaces; i_face++){

        // Retrieve normal 
        il::StaticArray<double, 3> face_n = {il::value, 
            {faceNormals(i_face, 0), faceNormals(i_face, 1), faceNormals(i_face, 2)}
        };

        // Loop on edges 
        for (int i_edge(0); i_edge<numVertPerFace; i_edge++){

            // vertices indices
            size_t a_I_i = faceIndices(i_face, i_edge);
            size_t a_Ip1_i = faceIndices(i_face, (i_edge+1)%numVertPerFace);

            // Get the vertices
            il::StaticArray<double, 3> A_I{il::value, {
                vertices(a_I_i, 0), 
                vertices(a_I_i, 1),
                vertices(a_I_i, 2)
            }};
            il::StaticArray<double, 3> A_Ip1{il::value, {
                vertices(a_Ip1_i, 0), 
                vertices(a_Ip1_i, 1),
                vertices(a_Ip1_i, 2)
            }};

            // Compute v = A_I A_I+1
            il::StaticArray<double, 3> v{il::value, {
                A_Ip1[0] - A_I[0], 
                A_Ip1[1] - A_I[1], 
                A_Ip1[2] - A_I[2]
            }};
            double v_norm = il::norm(v, il::Norm::L2);
            for (int i(0); i<3; i++) v[i] /= v_norm;

            // Compute b = n x v 
            il::StaticArray<double, 3> b = {il::value, {
                face_n[1]*v[2] - face_n[2]*v[1], 
                face_n[2]*v[0] - face_n[0]*v[2], 
                face_n[0]*v[1] - face_n[1]*v[0],
            }};
            double b_norm = il::norm(b, il::Norm::L2);
            for (int i(0); i<3; i++) b[i] /= b_norm;

            // Compute the two distance vectors 
            il::StaticArray<double, 3> r_1{il::value, {
                xy_obs[0] - A_I[0], 
                xy_obs[1] - A_I[1], 
                xy_obs[2] - A_I[2]
            }};
            il::StaticArray<double, 3> r_2{il::value, {
                xy_obs[0] - A_Ip1[0], 
                xy_obs[1] - A_Ip1[1], 
                xy_obs[2] - A_Ip1[2]
            }};
            double r_1_norm = il::norm(r_1, il::Norm::L2);
            double r_2_norm = il::norm(r_2, il::Norm::L2);

            // Compute projections
            double r_n = r_1[0]*face_n[0] + r_1[1]*face_n[1] + r_1[2]*face_n[2];
            double r_b = r_1[0]*b[0] + r_1[1]*b[1] + r_1[2]*b[2];
            double r_1_v = r_1[0]*v[0] + r_1[1]*v[1] + r_1[2]*v[2];
            double r_2_v = r_2[0]*v[0] + r_2[1]*v[1] + r_2[2]*v[2];

            // No contribution if r_b = 0
            if (std::abs(r_b) < eps) continue;

            // Compute I_-1 
            double I1 = std::log(r_2_norm - r_2_v) - std::log(r_1_norm - r_1_v);

            // If r_n == 0 it simplifies to :
            if (std::abs(r_n) < eps) {

                // Outer product
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        phi_ij[i][j] += face_n[i] * b[j] * I1; 

                continue;
            }

            // apex angle 
            double theta = std::atan(r_2_v / r_b) - std::atan(r_1_v / r_b);

            // compute J_-1
            double J1 = std::atan((r_n * r_2_v) / (r_b * r_2_norm))
                        - std::atan((r_n * r_1_v) / (r_b * r_1_norm));

            // second derivative update: phi_ij += outer(normal, ...)
            il::StaticArray<double, 3> coeff_vec{0.0};
            double sgn_rn = (r_n >= 0.0) ? 1.0 : -1.0;
            for (int i = 0; i < 3; ++i)
                coeff_vec[i] = b[i] * I1 + face_n[i] * (J1 - sgn_rn * theta);

            // phi_ij += np.einsum("i,j->ij", normal, coeff_vec)
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    phi_ij[i][j] += face_n[i] * coeff_vec[j];

        } // Loop on edges 
    } // Loop on polygonal faces

    return phi_ij;
}

} // namespace bigwham