//
// Created by Federico Ciardo on 16.08.21.
//

// Inclusion from IL library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linearAlgebra.h>

#ifndef INC_3DEQSIM_SRC_ELASTICKERNELSP4_H
#define INC_3DEQSIM_SRC_ELASTICKERNELSP4_H

namespace bie{

il::Array2D<double> StressTensorDueToDDx11P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDx12P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDx13P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDx21P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDx22P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDx23P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDx31P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDx32P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDx33P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy11P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy12P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy13P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy21P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy22P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy23P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy31P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy32P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDy33P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz11P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz12P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz13P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz21P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz22P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz23P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz31P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz32P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

il::Array2D<double> StressTensorDueToDDz33P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G);

}

#endif  // INC_3DEQSIM_SRC_ELASTICKERNELSP4_H
