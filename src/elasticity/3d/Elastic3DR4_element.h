//
// Created by Federico Ciardo on 11.08.21.
//

// Inclusion from IL library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linearAlgebra.h>

// Inclusion from the project
#include <src/core/ElasticProperties.h>

#ifndef INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H
#define INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H

namespace bie{

il::Array2D<double> TractionsDueToDDsOnSingleEltP4(
    il::Array<double> &a_mapped, il::Array<double> &b_mapped, il::Array2D<double> &Xe_mapped,
    il::Array2D<double> &shear1_vector_mapped, il::Array2D<double> &shear2_vector_mapped,
    il::Array2D<double> &normal_vector_mapped,
    ElasticProperties Matrix_Prop, il::io_t);

}  // namespace bie

#endif  // INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H
