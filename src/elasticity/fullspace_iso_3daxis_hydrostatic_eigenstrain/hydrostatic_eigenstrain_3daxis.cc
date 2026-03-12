#include "hydrostatic_eigenstrain_3daxis.hh"

#include <iostream>
#include <cmath>
#include <numbers>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include "../fullspace_iso_axisymmetry_flat_unidirectional/elliptic_integral.h"

// #define DEBUG

namespace bigwham {

// Parameters to pass to the integrand
struct IntegrandNMParams {
    double R;
    double Z;
    double a_start;
    double z_start;
    double da_dt;
    double dz_dt;
};

// Gauss-Legendre quadrature nodes and weights (10-point rule on [-1, 1])
struct GaussLegendreQuadrature {
    static constexpr int n_points = 10;
    static constexpr double nodes[n_points] = {
        -0.9739065285171717,
        -0.8650633666889845,
        -0.6794095682990244,
        -0.4333953941292472,
        -0.1488743389816312,
         0.1488743389816312,
         0.4333953941292472,
         0.6794095682990244,
         0.8650633666889845,
         0.9739065285171717
    };
    static constexpr double weights[n_points] = {
        0.0666713443086881,
        0.1494513491505806,
        0.2190863625159820,
        0.2692667193099963,
        0.2955242247147529,
        0.2955242247147529,
        0.2692667193099963,
        0.2190863625159820,
        0.1494513491505806,
        0.0666713443086881
    };
};

// Helper function to select integration method based on distance criteria
IntegrationMethod selectIntegrationMethod(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs
){
    double R = rz_obs[0];
    double Z = rz_obs[1];

    // Get the polygon vertices
    const int num_vertices = polygon.num_vertices();
    auto vertices_il_array = polygon.vertices();

    // Get centroid
    auto centroid = polygon.centroid();
    double a_centroid = centroid[0];
    double z_centroid = centroid[1];

    // Compute distance to centroid 
    double d = std::sqrt(std::pow(R - a_centroid, 2) + std::pow(Z - z_centroid, 2));

    // Compute elements characteristic size = sqrt of area
    double a_elmt = std::sqrt(polygon.size());

    // Compute the selection ratio
    double ratio = d / a_elmt;  // Add small epsilon to avoid division by zero

    // Select method based on ratio
    if (ratio < 1.0){
        return IntegrationMethod::QUADPACK;
    } else if (ratio < 20.0){
        return IntegrationMethod::GAUSS_LEGENDRE;
    } else {
        return IntegrationMethod::FAR_FIELD;
    }
}

// Helper function to compute far-field approximation
double phi_3daxis_far_field(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs
){
    double R = rz_obs[0];
    double Z = rz_obs[1];

    // Get area (already computed via shoelace formula)
    double area = polygon.size();

    // Get centroid
    auto centroid = polygon.centroid();
    double a_centroid = centroid[0];
    double z_centroid = centroid[1];

    // Compute k at centroid
    double k = 2 * std::sqrt(a_centroid * R) / std::sqrt(std::pow(R + a_centroid, 2) + std::pow(Z - z_centroid, 2));

    // Complete elliptic integral of the first kind
    double K = elliptic_fk(k);

    // Far-field approximation: φ ≈ -area * 2 * √(a/R) * k * K(k)
    return -area * 2 * std::sqrt(a_centroid / R) * k * K;
}

// Integrand to integrate over the polygon edges
double integrandNM(double t, void* p) {

    // retrieve parameters
    auto* params = static_cast<IntegrandNMParams*>(p);
    double& R = params->R;
    double& Z = params->Z;
    double& da_dt = params->da_dt;
    double& dz_dt = params->dz_dt;
    double& a_start = params->a_start;
    double& z_start = params->z_start;

    // Source coordinates
    double a = a_start + t * da_dt;
    double z = z_start + t * dz_dt;

    // Zero if on the polar axis
    if ((R <= 0) || (a <= 0))
        return 0.0;

    // Compute k
    double k = 2 * std::sqrt(a * R) / std::sqrt(std::pow(R + a, 2) + std::pow(Z - z, 2));

    // Complete elliptic integral of the first kind: K(k)
    // Note: elliptic_fk and elliptic_ek take the modulus k, not k²
    double K = elliptic_fk(k);
    double E = elliptic_ek(k);

    // Compute kernel M
    double M = (Z - z) * std::sqrt(a / R) * k * K;

    // Compute kernel N
    double term1 = (a + R) * K;
    double term2 = (2 * R / (k * k)) * (K - E);
    double N = std::sqrt(a / R) * k * (term1 - term2);

    // Return result
    return M * da_dt + N * dz_dt;
}


double phi_3daxis(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs,
    IntegrationMethod method
){

    double R = rz_obs[0];
    double Z = rz_obs[1];

    if (R < 0){
        std::cerr << "phi_3daxis called with (R,Z)=(" << R << ", " << Z << ")\n";
        throw std::runtime_error("3D axisymmetrical eigenstrain potential: R must be positive");
    }

    // Use far-field approximation if requested
    if (method == IntegrationMethod::FAR_FIELD){
        return phi_3daxis_far_field(polygon, rz_obs);
    }

    // Result
    double phi = 0.0;

    // Get the polygon vertices
    const int num_vertices = polygon.num_vertices();
    auto vertices_il_array = polygon.vertices();

    // Loop on edges
    for (int i_edge(0); i_edge<num_vertices; i_edge++){

        // Get coor start and end
        double a_start = vertices_il_array(i_edge, 0);
        double z_start = vertices_il_array(i_edge, 1);
        double a_end = vertices_il_array((i_edge+1)%num_vertices, 0);
        double z_end = vertices_il_array((i_edge+1)%num_vertices, 1);

        // Parametrize the edge: (a(t), z(t)) = start + t * (end - start), t ∈ [0, 1]
        double da_dt = a_end - a_start;
        double dz_dt = z_end - z_start;

        // Put together the integration parameters
        IntegrandNMParams params{R, Z, a_start, z_start, da_dt, dz_dt};

        double result;

        if (method == IntegrationMethod::QUADPACK){
            // Adaptive QAGS integration (scipy quad equivalent)
            gsl_function F;
            F.function = &integrandNM;
            F.params = &params;

            gsl_integration_workspace* w  = gsl_integration_workspace_alloc(1000);
            double error;

            int status = gsl_integration_qags(&F, 0, 1, 1.49e-8, 1.49e-8, 1000, w, &result, &error);

            gsl_integration_workspace_free(w);

            if (status != GSL_SUCCESS){
                throw std::runtime_error("Error: QAGS failed, the integrand is likely singular.");
            }
        }
        else if (method == IntegrationMethod::GAUSS_LEGENDRE){
            // Gauss-Legendre quadrature
            // Transform from [-1, 1] to [0, 1]: t = (s + 1)/2, dt = ds/2
            result = 0.0;
            for (int i = 0; i < GaussLegendreQuadrature::n_points; i++){
                double s = GaussLegendreQuadrature::nodes[i];
                double w = GaussLegendreQuadrature::weights[i];
                double t = (s + 1.0) / 2.0;  // Transform to [0, 1]
                result += w * integrandNM(t, &params);
            }
            result *= 0.5;  // Jacobian for transformation
        }

        // Add edge contribution
        phi += result;
    }

    return -phi;
}

std::array<double, 4> strain_3daxis(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs,
    double nu,
    IntegrationMethod method
){

    // FD step
    double eps = 1e-4;

    // Compute function values for finite differences
    double phi_0 = phi_3daxis(polygon, rz_obs, method);

    double R = rz_obs[0];
    double Z = rz_obs[1];

    // Base stencil
    double phi_R_plus  = phi_3daxis(polygon, {il::value, {R+eps, Z}}, method);
    double phi_R_minus = phi_3daxis(polygon, {il::value, {R-eps, Z}}, method);
    double phi_Z_plus  = phi_3daxis(polygon, {il::value, {R, Z+eps}}, method);
    double phi_Z_minus = phi_3daxis(polygon, {il::value, {R, Z-eps}}, method);

    // First order derivative
    double dR = (phi_R_plus - phi_R_minus) / (2 * eps);
    // double dZ = (phi_Z_plus - phi_Z_minus) / (2 * eps);

    // Corner points for rz component
    double phi_RpZp = phi_3daxis(polygon, {il::value, {R+eps, Z+eps}}, method);
    double phi_RpZm = phi_3daxis(polygon, {il::value, {R+eps, Z-eps}}, method);
    double phi_RmZp = phi_3daxis(polygon, {il::value, {R-eps, Z+eps}}, method);
    double phi_RmZm = phi_3daxis(polygon, {il::value, {R-eps, Z-eps}}, method);

    // Second derivatives
    double dRR = (phi_R_plus - 2 * phi_0 + phi_R_minus) / (eps*eps);
    double dZZ = (phi_Z_plus - 2 * phi_0 + phi_Z_minus) / (eps*eps);
    double dRZ = (phi_RpZp - phi_RpZm - phi_RmZp + phi_RmZm) / (4 * eps * eps);

    // Prefactor
    double factor = -(1 / (4 * std::numbers::pi)) * (1 + nu) / (1 - nu);

    // Strain components
    double err = factor * dRR;
    double ezz = factor * dZZ;
    double ett = factor * dR / R;
    double erz = factor * dRZ;

    return {err, ezz, ett, erz};
}

std::array<double, 4> stress_3daxis(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs,
    double nu, double G,
    IntegrationMethod method
){
    // Lame parameter
    double lam = 2 * G * nu / (1 - 2*nu);

    // Compute strain
    auto e_ij = strain_3daxis(polygon, rz_obs, nu, method);

    // If field point in inclusion: taking eigenstrain into account
    if (polygon.isPointInPolygon({rz_obs[0],rz_obs[1]})){
        e_ij[0] += 1.0; // Should be consistent in the convention taken
        e_ij[1] += 1.0;
        e_ij[2] += 1.0;
    }

    // Volumetric strain
    double eps_kk = e_ij[0] + e_ij[1] + e_ij[2];

    double srr = lam * eps_kk + 2 * G * e_ij[0];
    double szz = lam * eps_kk + 2 * G * e_ij[1];
    double stt = lam * eps_kk + 2 * G * e_ij[2];
    double srz = 2 * G * e_ij[3];

    return {srr, szz, stt, srz};
}

il::StaticArray<double, 2> V_threeDAxis_polygon_0(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs,
    const il::StaticArray<double, 2> n_obs,
    double G, double nu,
    bool recursive_call
){

    #ifdef DEBUG
        std::cerr << "Calling V_threeDAxis_polygon_0 at (R,Z) = (" << rz_obs[0] << ", " <<  rz_obs[1] << ")\n";
    #endif

    // Return value
    il::StaticArray<double, 2> t_i;

    // Determine if xy_obs is on the boundary of the triangle
    bool is_rz_on_boundary = polygon.isPointOnBoundary({rz_obs[0], rz_obs[1]});

    // Select integration method automatically (only once, not in recursive calls)
    IntegrationMethod method = IntegrationMethod::QUADPACK;  // Default
    if (!recursive_call){
        method = selectIntegrationMethod(polygon, rz_obs);
    }

    // If on the boundary :
    if (is_rz_on_boundary){
        // Evaluate the traction at two points on each side and take the average

        // std::cerr << "Point on the boundary\n";

        if (recursive_call) {
            throw std::runtime_error("Error: the point is still on the face when moving along the provided normal.");
        }

        double eps = 1e3 * polygon.getTol();

        il::StaticArray<double, 2> rz_obs_p{il::value, {
            rz_obs[0] + eps*n_obs[0],
            rz_obs[1] + eps*n_obs[1]}
        };
        il::StaticArray<double, 2> rz_obs_m{il::value, {
            rz_obs[0] - eps*n_obs[0],
            rz_obs[1] - eps*n_obs[1]}
        };

        // Ensuring that not both are considered as within the inclusion (tolerance issue)
        if (polygon.isPointOnBoundary({rz_obs_p[0],rz_obs_p[1]}) && polygon.isPointOnBoundary({rz_obs_m[0],rz_obs_m[1]}))
            throw std::runtime_error("In V_twoD_polygon_0 : both points are considered as within the inclusion.");

        auto t_p = V_threeDAxis_polygon_0(polygon, rz_obs_p, n_obs, G, nu, true);
        auto t_m = V_threeDAxis_polygon_0(polygon, rz_obs_m, n_obs, G, nu, true);

        t_i[0] = (t_p[0] + t_m[0])/2;
        t_i[1] = (t_p[1] + t_m[1])/2;
    }
    else {

        auto s_ij = stress_3daxis(polygon, rz_obs, nu, G, method);

        t_i[0] = s_ij[0] * n_obs[0] + s_ij[3] * n_obs[1];
        t_i[1] = s_ij[3] * n_obs[0] + s_ij[1] * n_obs[1];
    }

    return t_i;
}





} // namespace bigwham 