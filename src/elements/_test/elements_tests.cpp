//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 11.01.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2024.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <gtest/gtest.h>
#include <il/Array.h>


#include "elements/boundary_element.h"
#include "elements/point.h"

TEST(single_ELT_pt,point_2_1){

    il::Array2D<double> my_pts{1,2,0.};

    bie::Point<2> test;
    test.SetElement(my_pts);

    ASSERT_TRUE(test.num_vertices()==1);

}


TEST(single_ELT_pt,point_2_2){

    il::Array2D<double> my_pts{1,2,0.};

    bie::Point<2> test;
    test.SetElement(my_pts);

    ASSERT_TRUE(test.spatial_dimension()==2);
}

TEST(single_ELT_pt,point_2_3){

    il::Array2D<double> my_pts{1,2,0.};

    bie::Point<2> test;
    test.SetElement(my_pts);
    auto pts = test.nodes();
    bool tt =true;
    for (int i=0;i<test.num_vertices();i++){
        for (int j=0;j<test.spatial_dimension();j++){
           tt= tt && ( pts(i,j)==my_pts(i,j));
        }
    }
    ASSERT_TRUE(tt);
}


TEST(single_ELT_pt,point_2_4){

    il::Array2D<double> my_pts{1,2,0.};

    bie::Point<2> test;
    test.SetElement(my_pts);
    auto pts = test.collocation_points();
    bool tt =true;
    for (int i=0;i<test.num_vertices();i++){
        for (int j=0;j<test.spatial_dimension();j++){
            tt= tt && ( pts(i,j)==my_pts(i,j));
        }
    }
    ASSERT_TRUE(tt);
}



TEST(single_ELT_pt,point_3_1){

    il::Array2D<double> my_pts{1,3,2.};

    bie::Point<3> test;
    test.SetElement(my_pts);

    ASSERT_TRUE(test.num_vertices()==1);

}


TEST(single_ELT_pt,point_3_2){

    il::Array2D<double> my_pts{1,3,3.};

    bie::Point<3> test;
    test.SetElement(my_pts);

    ASSERT_TRUE(test.spatial_dimension()==3);
}

TEST(single_ELT_pt,point_3_3){

    il::Array2D<double> my_pts{1,3,0.};

    bie::Point<3> test;
    test.SetElement(my_pts);
    auto pts = test.nodes();
    bool tt =true;
    for (int i=0;i<test.num_vertices();i++){
        for (int j=0;j<test.spatial_dimension();j++){
            tt= tt && ( pts(i,j)==my_pts(i,j));
        }
    }
    ASSERT_TRUE(tt);
}


TEST(single_ELT_pt,point_3_4){

    il::Array2D<double> my_pts{1,3,0.};

    bie::Point<3> test;
    test.SetElement(my_pts);
    auto pts = test.collocation_points();
    bool tt =true;
    for (int i=0;i<test.num_vertices();i++){
        for (int j=0;j<test.spatial_dimension();j++){
            tt= tt && ( pts(i,j)==my_pts(i,j));
        }
    }
    ASSERT_TRUE(tt);
}

