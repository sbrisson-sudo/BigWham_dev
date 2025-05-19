//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 04.07.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <il/math.h>

#include "elastic_axi3dP0_element.h"
#include "elliptic_integral.h"

namespace bigwham{

// radial coordinate of the observation point
// radius of the dislocation disk
    double stress_disk_dislocation(double &rObs,double &rSrc)  {
        double IF;
        if (rSrc > 0.0) {
            double rho = rObs / rSrc;
            double k = 2.0 * pow(rho, 0.5) / (1.0 + rho);
            double k1 = (1.0 - rho) / (1.0 + rho);
            double P101 =
                    (k / (2.0 * il::pi * pow(rho, 0.5))) *
                    ((pow(k, 2.0) * (1.0 - pow(rho, 2.0)) / (4.0 * rho * pow(k1, 2.0))) *
                     elliptic_em(pow(k, 2.0)) +
                     elliptic_fm(pow(k, 2.0)));
            IF = P101 / (4.0 * rSrc);
        }  else {
            IF = 0.0;
        }
        return IF;
    }


}