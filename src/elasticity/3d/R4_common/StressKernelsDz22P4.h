//
// Created by Federico Ciardo on 17.08.21.
//

// Inclusion from IL library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linearAlgebra.h>

#ifndef INC_3DEQSIM_SRC_STRESSKERNELSDZ22P4_H
#define INC_3DEQSIM_SRC_STRESSKERNELSDZ22P4_H

namespace bie{

il::Array<double> StressComponentsDueToDDz22P4(double &a, double &b,
                                               double &x1, double &x2,
                                               double &x3, double &Nu,
                                               double &G);
}

#endif  // INC_3DEQSIM_SRC_STRESSKERNELSDZ22P4_H
