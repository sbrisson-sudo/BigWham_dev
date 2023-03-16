//
// Created by Alexis Sáez Uribe on 30/05/2022.
// Edited by Ankit on 16 March 2023
//

#ifndef BIGWHAM_ELASTICAXI3DP0_ELEMENT_H
#define BIGWHAM_ELASTICAXI3DP0_ELEMENT_H


#include <il/StaticArray2D.h>
#include <il/StaticArray.h>

#include <core/ElasticProperties.h>
#include <core/Mesh2D.h>
#include <core/SegmentData.h>

namespace bie {

    il::StaticArray2D<double, 2, 2> traction_influence_Axi3DP0(
            SegmentData &source_elt, SegmentData &receiver_elt, const ElasticProperties &Elas);

    double stress_disk_dislocation( double rObs, double rSrc );

}

#endif //BIGWHAM_ELASTICAXI3DP0_ELEMENT_H
