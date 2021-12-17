//
// Created by Federico Ciardo on 11.08.21.
//

// Inclusion from IL library
#include <cluster/cluster.h>
#include <il/Array.h>
#include <il/Array2D.h>

#ifndef INC_3DEQSIM_SRC_UTILS_H
#define INC_3DEQSIM_SRC_UTILS_H

namespace bie{

//inline il::Array2D<double> RotationMatrix3D(il::Array<double> &normal_vector,
//                                            double &theta);

inline double ArcCoth(double x);
inline std::complex<double> LogComplex(double x);

template <class T>
inline il::Array<T> row_selection(const il::Array2D<T> &arr, il::int_t idx);

}  // namespace bie

#endif  // INC_3DEQSIM_SRC_UTILS_H
