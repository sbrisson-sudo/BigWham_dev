//
// Created by Federico Ciardo on 11.08.21.
//

// Inclusion from IL library
#include <cluster/cluster.h>
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from the project
#include "ElementData.h"
#include "Mesh.h"

#ifndef INC_3DEQSIM_SRC_UTILS_H
#define INC_3DEQSIM_SRC_UTILS_H

namespace EQSim {

inline il::Array2D<double> RotationMatrix3D(il::Array<double> &normal_vector,
                                            double &theta);

inline double ArcCoth(double x);
inline std::complex<double> LogComplex(double x);

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_UTILS_H
