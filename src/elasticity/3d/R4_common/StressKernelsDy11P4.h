//
// Created by Federico Ciardo on 16.08.21.
//

// Inclusion from IL library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linearAlgebra.h>

#ifndef INC_3DEQSIM_SRC_STRESSKERNELSDY11P4_H
#define INC_3DEQSIM_SRC_STRESSKERNELSDY11P4_H

namespace EQSim{

il::Array<double> StressComponentsDueToDDy11P4(double &a, double &b,
                                               double &x1, double &x2,
                                               double &x3, double &Nu,
                                               double &G);

}



#endif  // INC_3DEQSIM_SRC_STRESSKERNELSDY11P4_H
