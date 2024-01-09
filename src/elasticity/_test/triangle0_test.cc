//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 01.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <cstdlib>
#include <memory>

#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "core/elastic_properties.h"
#include "elasticity/bie_elastostatic.h"
#include "elements/triangle.h"
#include "elements/point.h"


TEST(Triangle0, test_H_1) {
  // test self-effect triangle 0
  il::Array2D<double> xy{3, 3, 0.};
  xy(1, 0) = 1;
  xy(2, 1) = 1.0;
  bie::Triangle<0> source;
  source.SetElement(xy);
  bie::ElasticProperties elas(1, 0.3);
  bie::BieElastostatic<bie::Triangle<0>, bie::Triangle<0>,
                       bie::ElasticKernelType::H>
      test(elas, xy.size(1));
  std::vector<double> test_self = test.influence(source, 0, source, 0);
  // std::cout << "test self effect "
  //           << "\n";
  // for (int i = 0; i < 9; i++) {
  //   std::cout << test_self[i] << "\n";
  // }
  ASSERT_NEAR(test_self[0], 0.656304, 1.e-3);
}


TEST(Triangle0, test_stress_1) {
    // test stress observation single element
    il::Array2D<double> xy{3, 3, 0.};
    xy(1, 0) = 1;
    xy(2, 1) = 1.0;

    bie::Triangle<0> source;
    source.SetElement(xy);
    il::Array2D<double> xobs{1,3,1.};
    bie::Point<3> obs;
    obs.SetElement(xobs);
    bie::ElasticProperties elas(1, 0.3);
    bie::BieElastostatic<bie::Triangle<0>, bie::Point<3>,
          bie::ElasticKernelType::W> singleT(elas, xy.size(1));
    std::vector<double> test_stress = singleT.influence(source, 0, obs, 0);
    std::cout << "test stress "      << "\n";
    bool test=true;

    for (int i = 0; i < 18; i++) {
       std::cout << "i: " << i << " - " << test_stress[i] << "\n";
     }
    for (int i = 12; i < 18; i++) {
        test = test && (test_stress[i]==0);
    }

    ASSERT_TRUE(test );
    std::cout <<"end test Triangle0.test_stress_1" <<"\n";
}



