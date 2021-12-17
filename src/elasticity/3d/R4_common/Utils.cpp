//
// Created by Federico Ciardo on 11.08.21.
//

#include "Utils.h"

namespace bie {
//
//// Rotation matrix in 3D
//inline il::Array2D<double> RotationMatrix3D(il::Array<double> &normal_vector,
//                                            double &theta) {
//  il::Array2D<double> R{3, 3, 0.};
//
//  double n1, n2, n3, n1square, n2square, n3square;
//  n1 = normal_vector[0];
//  n2 = normal_vector[1];
//  n3 = normal_vector[2];
//  n1square = n1 * n1;
//  n2square = n2 * n2;
//  n3square = n3 * n3;
//
//  R(0, 0) = cos(theta) + (n1square * (1 - cos(theta)));
//  R(0, 1) = (n1 * n2 * (1 - cos(theta))) - (n3 * sin(theta));
//  R(0, 2) = (n1 * n3 * (1 - cos(theta))) + (n2 * sin(theta));
//
//  R(1, 0) = (n1 * n2 * (1 - cos(theta))) + n3 * sin(theta);
//  R(1, 1) = cos(theta) + (n2square * (1 - cos(theta)));
//  R(1, 2) = (n2 * n3 * (1 - cos(theta))) - (n1 * sin(theta));
//
//  R(2, 0) = (n1 * n3 * (1 - cos(theta))) - (n2 * sin(theta));
//  R(2, 1) = (n2 * n3 * (1 - cos(theta))) + (n1 * sin(theta));
//  R(2, 2) = cos(theta) + (n3square * (1 - cos(theta)));
//
//  return R;
//}

/////////
// Auxiliary function for assembly process
// It returns a given row (vector - specified by idx) of a 2D array
    template <typename T>
    inline il::Array<T> row_selection(const il::Array2D<T> &arr, il::int_t idx) {
        il::int_t arr_size1 = arr.size(1);
        il::Array<T> vect{arr_size1, 0};
        for (il::int_t i = 0; i < arr_size1; ++i) {
            vect[i] = arr(idx, i);
        }
        return vect;
    }

/////////
// Auxiliary function that calculates the Arc Cotangent hyperbolic
inline double ArcCoth(double x) {
  double res = atanh(1. / x);
  return res;
}

/////////
// Auxiliary function that calculates the natural log giving in output a complex
// number.
inline std::complex<double> LogComplex(double x) {
  std::complex<double> log_Arg = x;
  return log(log_Arg);
}

}  // namespace bie