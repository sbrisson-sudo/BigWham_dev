//
// Created by Federico Ciardo on 11.08.21.
//

// Inclusion from the project
#include "Elastic3DR4_element.h"

#include "R4_common/Stress_tensors_3DP4.h"
#include "R4_common/Utils.cpp"
#include <src/core/ElasticProperties.h>

namespace bie{

/// P4 elements
il::Array2D<double> TractionsDueToDDsOnSingleEltP4(
    il::Array<double> &a_mapped, il::Array<double> &b_mapped,
    il::Array2D<double> &Xe_mapped, il::Array2D<double> &shear1_vector_mapped,
    il::Array2D<double> &shear2_vector_mapped,
    il::Array2D<double> &normal_vector_mapped,
    ElasticProperties &Matrix_Prop,
    il::io_t) {

  // Poisson ratio
  double nu = Matrix_Prop.getNu();;
  // Shear modulus
  double Shear_Mod = Matrix_Prop.getG();

  il::Array<double> x1{Xe_mapped.size(0), 0.};
  il::Array<double> x2{Xe_mapped.size(0), 0.};
  il::Array<double> x3{Xe_mapped.size(0), 0.};

  for (il::int_t I = 0; I < Xe_mapped.size(0); ++I) {
    x1[I] = Xe_mapped(I, 0);
    x2[I] = Xe_mapped(I, 1);
    x3[I] = Xe_mapped(I, 2);
  }

  il::Array2D<double> StressTensorDueToDx11{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDx12{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDx13{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDx21{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDx22{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDx23{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDx31{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDx32{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDx33{3, 3, 0.};

  il::Array2D<double> StressTensorDueToDy11{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy12{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy13{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy21{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy22{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy23{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy31{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy32{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDy33{3, 3, 0.};

  il::Array2D<double> StressTensorDueToDz11{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz12{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz13{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz21{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz22{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz23{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz31{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz32{3, 3, 0.};
  il::Array2D<double> StressTensorDueToDz33{3, 3, 0.};

  StressTensorDueToDx11 = bie::StressTensorDueToDDx11P4(
      a_mapped[0], b_mapped[0], x1[0], x2[0], x3[0], nu, Shear_Mod);
  StressTensorDueToDx12 = bie::StressTensorDueToDDx12P4(
      a_mapped[1], b_mapped[1], x1[1], x2[1], x3[1], nu, Shear_Mod);
  StressTensorDueToDx13 = bie::StressTensorDueToDDx13P4(
      a_mapped[2], b_mapped[2], x1[2], x2[2], x3[2], nu, Shear_Mod);
  StressTensorDueToDx21 = bie::StressTensorDueToDDx21P4(
      a_mapped[3], b_mapped[3], x1[3], x2[3], x3[3], nu, Shear_Mod);
  StressTensorDueToDx22 = bie::StressTensorDueToDDx22P4(
      a_mapped[4], b_mapped[4], x1[4], x2[4], x3[4], nu, Shear_Mod);
  StressTensorDueToDx23 = bie::StressTensorDueToDDx23P4(
      a_mapped[5], b_mapped[5], x1[5], x2[5], x3[5], nu, Shear_Mod);
  StressTensorDueToDx31 = bie::StressTensorDueToDDx31P4(
      a_mapped[6], b_mapped[6], x1[6], x2[6], x3[6], nu, Shear_Mod);
  StressTensorDueToDx32 = bie::StressTensorDueToDDx32P4(
      a_mapped[7], b_mapped[7], x1[7], x2[7], x3[7], nu, Shear_Mod);
  StressTensorDueToDx33 = bie::StressTensorDueToDDx33P4(
      a_mapped[8], b_mapped[8], x1[8], x2[8], x3[8], nu, Shear_Mod);

  StressTensorDueToDy11 = bie::StressTensorDueToDDy11P4(
      a_mapped[0], b_mapped[0], x1[0], x2[0], x3[0], nu, Shear_Mod);
  StressTensorDueToDy12 = bie::StressTensorDueToDDy12P4(
      a_mapped[1], b_mapped[1], x1[1], x2[1], x3[1], nu, Shear_Mod);
  StressTensorDueToDy13 = bie::StressTensorDueToDDy13P4(
      a_mapped[2], b_mapped[2], x1[2], x2[2], x3[2], nu, Shear_Mod);
  StressTensorDueToDy21 = bie::StressTensorDueToDDy21P4(
      a_mapped[3], b_mapped[3], x1[3], x2[3], x3[3], nu, Shear_Mod);
  StressTensorDueToDy22 = bie::StressTensorDueToDDy22P4(
      a_mapped[4], b_mapped[4], x1[4], x2[4], x3[4], nu, Shear_Mod);
  StressTensorDueToDy23 = bie::StressTensorDueToDDy23P4(
      a_mapped[5], b_mapped[5], x1[5], x2[5], x3[5], nu, Shear_Mod);
  StressTensorDueToDy31 = bie::StressTensorDueToDDy31P4(
      a_mapped[6], b_mapped[6], x1[6], x2[6], x3[6], nu, Shear_Mod);
  StressTensorDueToDy32 = bie::StressTensorDueToDDy32P4(
      a_mapped[7], b_mapped[7], x1[7], x2[7], x3[7], nu, Shear_Mod);
  StressTensorDueToDy33 = bie::StressTensorDueToDDy33P4(
      a_mapped[8], b_mapped[8], x1[8], x2[8], x3[8], nu, Shear_Mod);

  StressTensorDueToDz11 = bie::StressTensorDueToDDz11P4(
      a_mapped[0], b_mapped[0], x1[0], x2[0], x3[0], nu, Shear_Mod);
  StressTensorDueToDz12 = bie::StressTensorDueToDDz12P4(
      a_mapped[1], b_mapped[1], x1[1], x2[1], x3[1], nu, Shear_Mod);
  StressTensorDueToDz13 = bie::StressTensorDueToDDz13P4(
      a_mapped[2], b_mapped[2], x1[2], x2[2], x3[2], nu, Shear_Mod);
  StressTensorDueToDz21 = bie::StressTensorDueToDDz21P4(
      a_mapped[3], b_mapped[3], x1[3], x2[3], x3[3], nu, Shear_Mod);
  StressTensorDueToDz22 = bie::StressTensorDueToDDz22P4(
      a_mapped[4], b_mapped[4], x1[4], x2[4], x3[4], nu, Shear_Mod);
  StressTensorDueToDz23 = bie::StressTensorDueToDDz23P4(
      a_mapped[5], b_mapped[5], x1[5], x2[5], x3[5], nu, Shear_Mod);
  StressTensorDueToDz31 = bie::StressTensorDueToDDz31P4(
      a_mapped[6], b_mapped[6], x1[6], x2[6], x3[6], nu, Shear_Mod);
  StressTensorDueToDz32 = bie::StressTensorDueToDDz32P4(
      a_mapped[7], b_mapped[7], x1[7], x2[7], x3[7], nu, Shear_Mod);
  StressTensorDueToDz33 = bie::StressTensorDueToDDz33P4(
      a_mapped[8], b_mapped[8], x1[8], x2[8], x3[8], nu, Shear_Mod);


  il::Array<double> n_1 = bie::row_selection(normal_vector_mapped, 0);
  il::Array<double> n_2 = bie::row_selection(normal_vector_mapped, 1);
  il::Array<double> n_3 = bie::row_selection(normal_vector_mapped, 2);
  il::Array<double> n_4 = bie::row_selection(normal_vector_mapped, 3);
  il::Array<double> n_5 = bie::row_selection(normal_vector_mapped, 4);
  il::Array<double> n_6 = bie::row_selection(normal_vector_mapped, 5);
  il::Array<double> n_7 = bie::row_selection(normal_vector_mapped, 6);
  il::Array<double> n_8 = bie::row_selection(normal_vector_mapped, 7);
  il::Array<double> n_9 = bie::row_selection(normal_vector_mapped, 8);

  il::Array<double> s1_1 = bie::row_selection(shear1_vector_mapped, 0);
  il::Array<double> s1_2 = bie::row_selection(shear1_vector_mapped, 1);
  il::Array<double> s1_3 = bie::row_selection(shear1_vector_mapped, 2);
  il::Array<double> s1_4 = bie::row_selection(shear1_vector_mapped, 3);
  il::Array<double> s1_5 = bie::row_selection(shear1_vector_mapped, 4);
  il::Array<double> s1_6 = bie::row_selection(shear1_vector_mapped, 5);
  il::Array<double> s1_7 = bie::row_selection(shear1_vector_mapped, 6);
  il::Array<double> s1_8 = bie::row_selection(shear1_vector_mapped, 7);
  il::Array<double> s1_9 = bie::row_selection(shear1_vector_mapped, 8);

  il::Array<double> s2_1 = bie::row_selection(shear2_vector_mapped, 0);
  il::Array<double> s2_2 = bie::row_selection(shear2_vector_mapped, 1);
  il::Array<double> s2_3 = bie::row_selection(shear2_vector_mapped, 2);
  il::Array<double> s2_4 = bie::row_selection(shear2_vector_mapped, 3);
  il::Array<double> s2_5 = bie::row_selection(shear2_vector_mapped, 4);
  il::Array<double> s2_6 = bie::row_selection(shear2_vector_mapped, 5);
  il::Array<double> s2_7 = bie::row_selection(shear2_vector_mapped, 6);
  il::Array<double> s2_8 = bie::row_selection(shear2_vector_mapped, 7);
  il::Array<double> s2_9 = bie::row_selection(shear2_vector_mapped, 8);

  // Traction vector on (x1,x2,x3) due to Dx in 11 etc..
  il::Array<double> tractionVectorDx11 = il::dot(StressTensorDueToDx11, n_1);
  il::Array<double> tractionVectorDx12 = il::dot(StressTensorDueToDx12, n_2);
  il::Array<double> tractionVectorDx13 = il::dot(StressTensorDueToDx13, n_3);
  il::Array<double> tractionVectorDx21 = il::dot(StressTensorDueToDx21, n_4);
  il::Array<double> tractionVectorDx22 = il::dot(StressTensorDueToDx22, n_5);
  il::Array<double> tractionVectorDx23 = il::dot(StressTensorDueToDx23, n_6);
  il::Array<double> tractionVectorDx31 = il::dot(StressTensorDueToDx31, n_7);
  il::Array<double> tractionVectorDx32 = il::dot(StressTensorDueToDx32, n_8);
  il::Array<double> tractionVectorDx33 = il::dot(StressTensorDueToDx33, n_9);

  // Traction vector on (x1,x2,x3) due to Dy in 11 etc..
  il::Array<double> tractionVectorDy11 = il::dot(StressTensorDueToDy11, n_1);
  il::Array<double> tractionVectorDy12 = il::dot(StressTensorDueToDy12, n_2);
  il::Array<double> tractionVectorDy13 = il::dot(StressTensorDueToDy13, n_3);
  il::Array<double> tractionVectorDy21 = il::dot(StressTensorDueToDy21, n_4);
  il::Array<double> tractionVectorDy22 = il::dot(StressTensorDueToDy22, n_5);
  il::Array<double> tractionVectorDy23 = il::dot(StressTensorDueToDy23, n_6);
  il::Array<double> tractionVectorDy31 = il::dot(StressTensorDueToDy31, n_7);
  il::Array<double> tractionVectorDy32 = il::dot(StressTensorDueToDy32, n_8);
  il::Array<double> tractionVectorDy33 = il::dot(StressTensorDueToDy33, n_9);

  // Traction vector on (x1,x2,x3) due to Dz in 11 etc..
  il::Array<double> tractionVectorDz11 = il::dot(StressTensorDueToDz11, n_1);
  il::Array<double> tractionVectorDz12 = il::dot(StressTensorDueToDz12, n_2);
  il::Array<double> tractionVectorDz13 = il::dot(StressTensorDueToDz13, n_3);
  il::Array<double> tractionVectorDz21 = il::dot(StressTensorDueToDz21, n_4);
  il::Array<double> tractionVectorDz22 = il::dot(StressTensorDueToDz22, n_5);
  il::Array<double> tractionVectorDz23 = il::dot(StressTensorDueToDz23, n_6);
  il::Array<double> tractionVectorDz31 = il::dot(StressTensorDueToDz31, n_7);
  il::Array<double> tractionVectorDz32 = il::dot(StressTensorDueToDz32, n_8);
  il::Array<double> tractionVectorDz33 = il::dot(StressTensorDueToDz33, n_9);

  // Initialization
  double ts1Dx11 = 0., ts1Dx12 = 0., ts1Dx13 = 0., ts1Dx21 = 0., ts1Dx22 = 0.,
         ts1Dx23 = 0., ts1Dx31 = 0., ts1Dx32 = 0., ts1Dx33 = 0.;
  double ts1Dy11 = 0., ts1Dy12 = 0., ts1Dy13 = 0., ts1Dy21 = 0., ts1Dy22 = 0.,
         ts1Dy23 = 0., ts1Dy31 = 0., ts1Dy32 = 0., ts1Dy33 = 0.;
  double ts1Dz11 = 0., ts1Dz12 = 0., ts1Dz13 = 0., ts1Dz21 = 0., ts1Dz22 = 0.,
         ts1Dz23 = 0., ts1Dz31 = 0., ts1Dz32 = 0., ts1Dz33 = 0.;

  double ts2Dx11 = 0., ts2Dx12 = 0., ts2Dx13 = 0., ts2Dx21 = 0., ts2Dx22 = 0.,
         ts2Dx23 = 0., ts2Dx31 = 0., ts2Dx32 = 0., ts2Dx33 = 0.;
  double ts2Dy11 = 0., ts2Dy12 = 0., ts2Dy13 = 0., ts2Dy21 = 0., ts2Dy22 = 0.,
         ts2Dy23 = 0., ts2Dy31 = 0., ts2Dy32 = 0., ts2Dy33 = 0.;
  double ts2Dz11 = 0., ts2Dz12 = 0., ts2Dz13 = 0., ts2Dz21 = 0., ts2Dz22 = 0.,
         ts2Dz23 = 0., ts2Dz31 = 0., ts2Dz32 = 0., ts2Dz33 = 0.;

  double tnDx11 = 0., tnDx12 = 0., tnDx13 = 0., tnDx21 = 0., tnDx22 = 0.,
         tnDx23 = 0., tnDx31 = 0., tnDx32 = 0., tnDx33 = 0.;
  double tnDy11 = 0., tnDy12 = 0., tnDy13 = 0., tnDy21 = 0., tnDy22 = 0.,
         tnDy23 = 0., tnDy31 = 0., tnDy32 = 0., tnDy33 = 0.;
  double tnDz11 = 0., tnDz12 = 0., tnDz13 = 0., tnDz21 = 0., tnDz22 = 0.,
         tnDz23 = 0., tnDz31 = 0., tnDz32 = 0., tnDz33 = 0.;

  for (int I = 0; I < tractionVectorDx11.size(); ++I) {
    // Shear traction ts1 on (x1,x2,x3) due to Dy in 11 etc...
    ts1Dx11 += s1_1[I] * tractionVectorDx11[I];
    ts1Dx12 += s1_2[I] * tractionVectorDx12[I];
    ts1Dx13 += s1_3[I] * tractionVectorDx13[I];
    ts1Dx21 += s1_4[I] * tractionVectorDx21[I];
    ts1Dx22 += s1_5[I] * tractionVectorDx22[I];
    ts1Dx23 += s1_6[I] * tractionVectorDx23[I];
    ts1Dx31 += s1_7[I] * tractionVectorDx31[I];
    ts1Dx32 += s1_8[I] * tractionVectorDx32[I];
    ts1Dx33 += s1_9[I] * tractionVectorDx33[I];

    // Shear traction ts1 on (x1,x2,x3) due to Dy in 11 etc...
    ts1Dy11 += s1_1[I] * tractionVectorDy11[I];
    ts1Dy12 += s1_2[I] * tractionVectorDy12[I];
    ts1Dy13 += s1_3[I] * tractionVectorDy13[I];
    ts1Dy21 += s1_4[I] * tractionVectorDy21[I];
    ts1Dy22 += s1_5[I] * tractionVectorDy22[I];
    ts1Dy23 += s1_6[I] * tractionVectorDy23[I];
    ts1Dy31 += s1_7[I] * tractionVectorDy31[I];
    ts1Dy32 += s1_8[I] * tractionVectorDy32[I];
    ts1Dy33 += s1_9[I] * tractionVectorDy33[I];

    // Shear traction ts1 on (x1,x2,x3) due to Dz in 11 etc...
    ts1Dz11 += s1_1[I] * tractionVectorDz11[I];
    ts1Dz12 += s1_2[I] * tractionVectorDz12[I];
    ts1Dz13 += s1_3[I] * tractionVectorDz13[I];
    ts1Dz21 += s1_4[I] * tractionVectorDz21[I];
    ts1Dz22 += s1_5[I] * tractionVectorDz22[I];
    ts1Dz23 += s1_6[I] * tractionVectorDz23[I];
    ts1Dz31 += s1_7[I] * tractionVectorDz31[I];
    ts1Dz32 += s1_8[I] * tractionVectorDz32[I];
    ts1Dz33 += s1_9[I] * tractionVectorDz33[I];

    // Shear traction ts2 on (x1,x2,x3) due to Dx in 11 etc...
    ts2Dx11 += s2_1[I] * tractionVectorDx11[I];
    ts2Dx12 += s2_2[I] * tractionVectorDx12[I];
    ts2Dx13 += s2_3[I] * tractionVectorDx13[I];
    ts2Dx21 += s2_4[I] * tractionVectorDx21[I];
    ts2Dx22 += s2_5[I] * tractionVectorDx22[I];
    ts2Dx23 += s2_6[I] * tractionVectorDx23[I];
    ts2Dx31 += s2_7[I] * tractionVectorDx31[I];
    ts2Dx32 += s2_8[I] * tractionVectorDx32[I];
    ts2Dx33 += s2_9[I] * tractionVectorDx33[I];

    // Shear traction ts2 on (x1,x2,x3) due to Dy in 11 etc...
    ts2Dy11 += s2_1[I] * tractionVectorDy11[I];
    ts2Dy12 += s2_2[I] * tractionVectorDy12[I];
    ts2Dy13 += s2_3[I] * tractionVectorDy13[I];
    ts2Dy21 += s2_4[I] * tractionVectorDy21[I];
    ts2Dy22 += s2_5[I] * tractionVectorDy22[I];
    ts2Dy23 += s2_6[I] * tractionVectorDy23[I];
    ts2Dy31 += s2_7[I] * tractionVectorDy31[I];
    ts2Dy32 += s2_8[I] * tractionVectorDy32[I];
    ts2Dy33 += s2_9[I] * tractionVectorDy33[I];

    // Shear traction ts2 on (x1,x2,x3) due to Dz in 11 etc...
    ts2Dz11 += s2_1[I] * tractionVectorDz11[I];
    ts2Dz12 += s2_2[I] * tractionVectorDz12[I];
    ts2Dz13 += s2_3[I] * tractionVectorDz13[I];
    ts2Dz21 += s2_4[I] * tractionVectorDz21[I];
    ts2Dz22 += s2_5[I] * tractionVectorDz22[I];
    ts2Dz23 += s2_6[I] * tractionVectorDz23[I];
    ts2Dz31 += s2_7[I] * tractionVectorDz31[I];
    ts2Dz32 += s2_8[I] * tractionVectorDz32[I];
    ts2Dz33 += s2_9[I] * tractionVectorDz33[I];

    // Normal traction tn on (x1,x2,x3) due to Dx in 11 etc...
    tnDx11 += n_1[I] * tractionVectorDx11[I];
    tnDx12 += n_2[I] * tractionVectorDx12[I];
    tnDx13 += n_3[I] * tractionVectorDx13[I];
    tnDx21 += n_4[I] * tractionVectorDx21[I];
    tnDx22 += n_5[I] * tractionVectorDx22[I];
    tnDx23 += n_6[I] * tractionVectorDx23[I];
    tnDx31 += n_7[I] * tractionVectorDx31[I];
    tnDx32 += n_8[I] * tractionVectorDx32[I];
    tnDx33 += n_9[I] * tractionVectorDx33[I];

    // Normal traction tn on (x1,x2,x3) due to Dy in 11 etc...
    tnDy11 += n_1[I] * tractionVectorDy11[I];
    tnDy12 += n_2[I] * tractionVectorDy12[I];
    tnDy13 += n_3[I] * tractionVectorDy13[I];
    tnDy21 += n_4[I] * tractionVectorDy21[I];
    tnDy22 += n_5[I] * tractionVectorDy22[I];
    tnDy23 += n_6[I] * tractionVectorDy23[I];
    tnDy31 += n_7[I] * tractionVectorDy31[I];
    tnDy32 += n_8[I] * tractionVectorDy32[I];
    tnDy33 += n_9[I] * tractionVectorDy33[I];

    // Normal traction tn on (x1,x2,x3) due to Dz in 11 etc...
    tnDz11 += n_1[I] * tractionVectorDz11[I];
    tnDz12 += n_2[I] * tractionVectorDz12[I];
    tnDz13 += n_3[I] * tractionVectorDz13[I];
    tnDz21 += n_4[I] * tractionVectorDz21[I];
    tnDz22 += n_5[I] * tractionVectorDz22[I];
    tnDz23 += n_6[I] * tractionVectorDz23[I];
    tnDz31 += n_7[I] * tractionVectorDz31[I];
    tnDz32 += n_8[I] * tractionVectorDz32[I];
    tnDz33 += n_9[I] * tractionVectorDz33[I];
  }

  il::Array2D<double> TractionsDueToDDsOnSingleEltP4{3, 3, 0.};

  // (ts1_Dx, ts1_Dy, ts1_Dz) -> ts1 due to Dx, Dy, Dz, respectively in all the
  // elements (superimposition)
  TractionsDueToDDsOnSingleEltP4(0, 0) = ts1Dx11 + ts1Dx12 + ts1Dx13 + ts1Dx21 +
                                         ts1Dx22 + ts1Dx23 + ts1Dx31 + ts1Dx32 +
                                         ts1Dx33;
  TractionsDueToDDsOnSingleEltP4(0, 1) = ts1Dy11 + ts1Dy12 + ts1Dy13 + ts1Dy21 +
                                         ts1Dy22 + ts1Dy23 + ts1Dy31 + ts1Dy32 +
                                         ts1Dy33;
  TractionsDueToDDsOnSingleEltP4(0, 2) = ts1Dz11 + ts1Dz12 + ts1Dz13 + ts1Dz21 +
                                         ts1Dz22 + ts1Dz23 + ts1Dz31 + ts1Dz32 +
                                         ts1Dz33;

  // (ts2_Dx, ts2_Dy, ts2_Dz) -> ts2 due to Dx, Dy, Dz, respectively in all the
  // elements (superimposition)
  TractionsDueToDDsOnSingleEltP4(1, 0) = ts2Dx11 + ts2Dx12 + ts2Dx13 + ts2Dx21 +
                                         ts2Dx22 + ts2Dx23 + ts2Dx31 + ts2Dx32 +
                                         ts2Dx33;
  TractionsDueToDDsOnSingleEltP4(1, 1) = ts2Dy11 + ts2Dy12 + ts2Dy13 + ts2Dy21 +
                                         ts2Dy22 + ts2Dy23 + ts2Dy31 + ts2Dy32 +
                                         ts2Dy33;
  TractionsDueToDDsOnSingleEltP4(1, 2) = ts2Dz11 + ts2Dz12 + ts2Dz13 + ts2Dz21 +
                                         ts2Dz22 + ts2Dz23 + ts2Dz31 + ts2Dz32 +
                                         ts2Dz33;

  // (tn_Dx, tn_Dy, tn_Dz) -> tn due to Dx, Dy, Dz, respectively in all the
  // elements (superimposition)
  TractionsDueToDDsOnSingleEltP4(2, 0) = tnDx11 + tnDx12 + tnDx13 + tnDx21 +
                                         tnDx22 + tnDx23 + tnDx31 + tnDx32 +
                                         tnDx33;
  TractionsDueToDDsOnSingleEltP4(2, 1) = tnDy11 + tnDy12 + tnDy13 + tnDy21 +
                                         tnDy22 + tnDy23 + tnDy31 + tnDy32 +
                                         tnDy33;
  TractionsDueToDDsOnSingleEltP4(2, 2) = tnDz11 + tnDz12 + tnDz13 + tnDz21 +
                                         tnDz22 + tnDz23 + tnDz31 + tnDz32 +
                                         tnDz33;

  return TractionsDueToDDsOnSingleEltP4;
}

}  // namespace bie