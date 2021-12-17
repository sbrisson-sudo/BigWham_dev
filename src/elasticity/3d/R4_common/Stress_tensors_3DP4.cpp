//
// Created by Federico Ciardo on 16.08.21.
//

// Inclusion from the project
#include "StressKernelsDx11P4.h"
#include "StressKernelsDx12P4.h"
#include "StressKernelsDx13P4.h"
#include "StressKernelsDx21P4.h"
#include "StressKernelsDx22P4.h"
#include "StressKernelsDx23P4.h"
#include "StressKernelsDx31P4.h"
#include "StressKernelsDx32P4.h"
#include "StressKernelsDx33P4.h"
#include "StressKernelsDy11P4.h"
#include "StressKernelsDy12P4.h"
#include "StressKernelsDy13P4.h"
#include "StressKernelsDy21P4.h"
#include "StressKernelsDy22P4.h"
#include "StressKernelsDy23P4.h"
#include "StressKernelsDy31P4.h"
#include "StressKernelsDy32P4.h"
#include "StressKernelsDy33P4.h"
#include "StressKernelsDz11P4.h"
#include "StressKernelsDz12P4.h"
#include "StressKernelsDz13P4.h"
#include "StressKernelsDz21P4.h"
#include "StressKernelsDz22P4.h"
#include "StressKernelsDz23P4.h"
#include "StressKernelsDz31P4.h"
#include "StressKernelsDz32P4.h"
#include "StressKernelsDz33P4.h"

namespace bie{

il::Array2D<double> StressTensorDueToDDx11P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx11 = 0., sigmayyDx11 = 0., sigmazzDx11 = 0., sigmaxyDx11 = 0.,
         sigmaxzDx11 = 0., sigmayzDx11 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDx11 = StressComponentsDueToDDx11P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx11;
  res(0, 1) = sigmaxyDx11;
  res(0, 2) = sigmaxzDx11;
  res(1, 0) = sigmaxyDx11;
  res(1, 1) = sigmayyDx11;
  res(1, 2) = sigmayzDx11;
  res(2, 0) = sigmaxzDx11;
  res(2, 1) = sigmayzDx11;
  res(2, 2) = sigmazzDx11;

  return res;
}

il::Array2D<double> StressTensorDueToDDx12P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx12 = 0., sigmayyDx12 = 0., sigmazzDx12 = 0., sigmaxyDx12 = 0.,
         sigmaxzDx12 = 0., sigmayzDx12 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
   sigmaxxDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[6];
   sigmayyDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[7];
   sigmazzDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[8];
   sigmaxyDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[9];
   sigmaxzDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[10];
   sigmayzDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
   sigmaxxDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[0];
   sigmayyDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[1];
   sigmazzDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[2];
   sigmaxyDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[3];
   sigmaxzDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[4];
   sigmayzDx12 = StressComponentsDueToDDx12P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx12;
  res(0, 1) = sigmaxyDx12;
  res(0, 2) = sigmaxzDx12;
  res(1, 0) = sigmaxyDx12;
  res(1, 1) = sigmayyDx12;
  res(1, 2) = sigmayzDx12;
  res(2, 0) = sigmaxzDx12;
  res(2, 1) = sigmayzDx12;
  res(2, 2) = sigmazzDx12;

  return res;
}

il::Array2D<double> StressTensorDueToDDx13P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx13 = 0., sigmayyDx13 = 0., sigmazzDx13 = 0., sigmaxyDx13 = 0.,
  sigmaxzDx13 = 0., sigmayzDx13 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDx13 = StressComponentsDueToDDx13P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx13;
  res(0, 1) = sigmaxyDx13;
  res(0, 2) = sigmaxzDx13;
  res(1, 0) = sigmaxyDx13;
  res(1, 1) = sigmayyDx13;
  res(1, 2) = sigmayzDx13;
  res(2, 0) = sigmaxzDx13;
  res(2, 1) = sigmayzDx13;
  res(2, 2) = sigmazzDx13;

  return res;
}

il::Array2D<double> StressTensorDueToDDx21P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx21 = 0., sigmayyDx21 = 0., sigmazzDx21 = 0., sigmaxyDx21 = 0.,
  sigmaxzDx21 = 0., sigmayzDx21 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDx21 = StressComponentsDueToDDx21P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx21;
  res(0, 1) = sigmaxyDx21;
  res(0, 2) = sigmaxzDx21;
  res(1, 0) = sigmaxyDx21;
  res(1, 1) = sigmayyDx21;
  res(1, 2) = sigmayzDx21;
  res(2, 0) = sigmaxzDx21;
  res(2, 1) = sigmayzDx21;
  res(2, 2) = sigmazzDx21;

  return res;
}

il::Array2D<double> StressTensorDueToDDx22P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx22 = 0., sigmayyDx22 = 0., sigmazzDx22 = 0., sigmaxyDx22 = 0.,
  sigmaxzDx22 = 0., sigmayzDx22 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDx22 = StressComponentsDueToDDx22P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx22;
  res(0, 1) = sigmaxyDx22;
  res(0, 2) = sigmaxzDx22;
  res(1, 0) = sigmaxyDx22;
  res(1, 1) = sigmayyDx22;
  res(1, 2) = sigmayzDx22;
  res(2, 0) = sigmaxzDx22;
  res(2, 1) = sigmayzDx22;
  res(2, 2) = sigmazzDx22;

  return res;
}

il::Array2D<double> StressTensorDueToDDx23P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx23 = 0., sigmayyDx23 = 0., sigmazzDx23 = 0., sigmaxyDx23 = 0.,
  sigmaxzDx23 = 0., sigmayzDx23 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDx23 = StressComponentsDueToDDx23P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx23;
  res(0, 1) = sigmaxyDx23;
  res(0, 2) = sigmaxzDx23;
  res(1, 0) = sigmaxyDx23;
  res(1, 1) = sigmayyDx23;
  res(1, 2) = sigmayzDx23;
  res(2, 0) = sigmaxzDx23;
  res(2, 1) = sigmayzDx23;
  res(2, 2) = sigmazzDx23;

  return res;
}

il::Array2D<double> StressTensorDueToDDx31P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx31 = 0., sigmayyDx31 = 0., sigmazzDx31 = 0., sigmaxyDx31 = 0.,
  sigmaxzDx31 = 0., sigmayzDx31 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDx31 = StressComponentsDueToDDx31P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx31;
  res(0, 1) = sigmaxyDx31;
  res(0, 2) = sigmaxzDx31;
  res(1, 0) = sigmaxyDx31;
  res(1, 1) = sigmayyDx31;
  res(1, 2) = sigmayzDx31;
  res(2, 0) = sigmaxzDx31;
  res(2, 1) = sigmayzDx31;
  res(2, 2) = sigmazzDx31;

  return res;
}

il::Array2D<double> StressTensorDueToDDx32P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx32 = 0., sigmayyDx32 = 0., sigmazzDx32 = 0., sigmaxyDx32 = 0.,
  sigmaxzDx32 = 0., sigmayzDx32 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDx32 = StressComponentsDueToDDx32P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx32;
  res(0, 1) = sigmaxyDx32;
  res(0, 2) = sigmaxzDx32;
  res(1, 0) = sigmaxyDx32;
  res(1, 1) = sigmayyDx32;
  res(1, 2) = sigmayzDx32;
  res(2, 0) = sigmaxzDx32;
  res(2, 1) = sigmayzDx32;
  res(2, 2) = sigmazzDx32;

  return res;
}

il::Array2D<double> StressTensorDueToDDx33P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDx33 = 0., sigmayyDx33 = 0., sigmazzDx33 = 0., sigmaxyDx33 = 0.,
  sigmaxzDx33 = 0., sigmayzDx33 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDx33 = StressComponentsDueToDDx33P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDx33;
  res(0, 1) = sigmaxyDx33;
  res(0, 2) = sigmaxzDx33;
  res(1, 0) = sigmaxyDx33;
  res(1, 1) = sigmayyDx33;
  res(1, 2) = sigmayzDx33;
  res(2, 0) = sigmaxzDx33;
  res(2, 1) = sigmayzDx33;
  res(2, 2) = sigmazzDx33;

  return res;
}

il::Array2D<double> StressTensorDueToDDy11P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy11 = 0., sigmayyDy11 = 0., sigmazzDy11 = 0., sigmaxyDy11 = 0.,
  sigmaxzDy11 = 0., sigmayzDy11 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDy11 = StressComponentsDueToDDy11P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy11;
  res(0, 1) = sigmaxyDy11;
  res(0, 2) = sigmaxzDy11;
  res(1, 0) = sigmaxyDy11;
  res(1, 1) = sigmayyDy11;
  res(1, 2) = sigmayzDy11;
  res(2, 0) = sigmaxzDy11;
  res(2, 1) = sigmayzDy11;
  res(2, 2) = sigmazzDy11;

  return res;
}

il::Array2D<double> StressTensorDueToDDy12P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy12 = 0., sigmayyDy12 = 0., sigmazzDy12 = 0., sigmaxyDy12 = 0.,
  sigmaxzDy12 = 0., sigmayzDy12 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDy12 = StressComponentsDueToDDy12P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy12;
  res(0, 1) = sigmaxyDy12;
  res(0, 2) = sigmaxzDy12;
  res(1, 0) = sigmaxyDy12;
  res(1, 1) = sigmayyDy12;
  res(1, 2) = sigmayzDy12;
  res(2, 0) = sigmaxzDy12;
  res(2, 1) = sigmayzDy12;
  res(2, 2) = sigmazzDy12;

  return res;
}

il::Array2D<double> StressTensorDueToDDy13P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy13 = 0., sigmayyDy13 = 0., sigmazzDy13 = 0., sigmaxyDy13 = 0.,
  sigmaxzDy13 = 0., sigmayzDy13 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDy13 = StressComponentsDueToDDy13P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy13;
  res(0, 1) = sigmaxyDy13;
  res(0, 2) = sigmaxzDy13;
  res(1, 0) = sigmaxyDy13;
  res(1, 1) = sigmayyDy13;
  res(1, 2) = sigmayzDy13;
  res(2, 0) = sigmaxzDy13;
  res(2, 1) = sigmayzDy13;
  res(2, 2) = sigmazzDy13;

  return res;
}

il::Array2D<double> StressTensorDueToDDy21P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy21 = 0., sigmayyDy21 = 0., sigmazzDy21 = 0., sigmaxyDy21 = 0.,
  sigmaxzDy21 = 0., sigmayzDy21 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDy21 = StressComponentsDueToDDy21P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy21;
  res(0, 1) = sigmaxyDy21;
  res(0, 2) = sigmaxzDy21;
  res(1, 0) = sigmaxyDy21;
  res(1, 1) = sigmayyDy21;
  res(1, 2) = sigmayzDy21;
  res(2, 0) = sigmaxzDy21;
  res(2, 1) = sigmayzDy21;
  res(2, 2) = sigmazzDy21;

  return res;
}

il::Array2D<double> StressTensorDueToDDy22P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy22 = 0., sigmayyDy22 = 0., sigmazzDy22 = 0., sigmaxyDy22 = 0.,
  sigmaxzDy22 = 0., sigmayzDy22 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDy22 = StressComponentsDueToDDy22P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy22;
  res(0, 1) = sigmaxyDy22;
  res(0, 2) = sigmaxzDy22;
  res(1, 0) = sigmaxyDy22;
  res(1, 1) = sigmayyDy22;
  res(1, 2) = sigmayzDy22;
  res(2, 0) = sigmaxzDy22;
  res(2, 1) = sigmayzDy22;
  res(2, 2) = sigmazzDy22;

  return res;
}

il::Array2D<double> StressTensorDueToDDy23P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy23 = 0., sigmayyDy23 = 0., sigmazzDy23 = 0., sigmaxyDy23 = 0.,
  sigmaxzDy23 = 0., sigmayzDy23 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDy23 = StressComponentsDueToDDy23P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy23;
  res(0, 1) = sigmaxyDy23;
  res(0, 2) = sigmaxzDy23;
  res(1, 0) = sigmaxyDy23;
  res(1, 1) = sigmayyDy23;
  res(1, 2) = sigmayzDy23;
  res(2, 0) = sigmaxzDy23;
  res(2, 1) = sigmayzDy23;
  res(2, 2) = sigmazzDy23;

  return res;
}

il::Array2D<double> StressTensorDueToDDy31P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy31 = 0., sigmayyDy31 = 0., sigmazzDy31 = 0., sigmaxyDy31 = 0.,
  sigmaxzDy31 = 0., sigmayzDy31 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDy31 = StressComponentsDueToDDy31P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy31;
  res(0, 1) = sigmaxyDy31;
  res(0, 2) = sigmaxzDy31;
  res(1, 0) = sigmaxyDy31;
  res(1, 1) = sigmayyDy31;
  res(1, 2) = sigmayzDy31;
  res(2, 0) = sigmaxzDy31;
  res(2, 1) = sigmayzDy31;
  res(2, 2) = sigmazzDy31;

  return res;
}

il::Array2D<double> StressTensorDueToDDy32P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy32 = 0., sigmayyDy32 = 0., sigmazzDy32 = 0., sigmaxyDy32 = 0.,
  sigmaxzDy32 = 0., sigmayzDy32 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmaxzDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDy32 = StressComponentsDueToDDy32P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy32;
  res(0, 1) = sigmaxyDy32;
  res(0, 2) = sigmaxzDy32;
  res(1, 0) = sigmaxyDy32;
  res(1, 1) = sigmayyDy32;
  res(1, 2) = sigmayzDy32;
  res(2, 0) = sigmaxzDy32;
  res(2, 1) = sigmayzDy32;
  res(2, 2) = sigmazzDy32;

  return res;
}

il::Array2D<double> StressTensorDueToDDy33P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDy33 = 0., sigmayyDy33 = 0., sigmazzDy33 = 0., sigmaxyDy33 = 0.,
  sigmaxzDy33 = 0., sigmayzDy33 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxyDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmaxzDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[10];
    sigmayzDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[11];

  } else {
    sigmaxxDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxyDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[3];
    sigmayzDy33 = StressComponentsDueToDDy33P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDy33;
  res(0, 1) = sigmaxyDy33;
  res(0, 2) = sigmaxzDy33;
  res(1, 0) = sigmaxyDy33;
  res(1, 1) = sigmayyDy33;
  res(1, 2) = sigmayzDy33;
  res(2, 0) = sigmaxzDy33;
  res(2, 1) = sigmayzDy33;
  res(2, 2) = sigmazzDy33;

  return res;
}

il::Array2D<double> StressTensorDueToDDz11P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz11 = 0., sigmayyDz11 = 0., sigmazzDz11 = 0., sigmaxyDz11 = 0.,
  sigmaxzDz11 = 0., sigmayzDz11 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz11 = StressComponentsDueToDDz11P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz11;
  res(0, 1) = sigmaxyDz11;
  res(0, 2) = sigmaxzDz11;
  res(1, 0) = sigmaxyDz11;
  res(1, 1) = sigmayyDz11;
  res(1, 2) = sigmayzDz11;
  res(2, 0) = sigmaxzDz11;
  res(2, 1) = sigmayzDz11;
  res(2, 2) = sigmazzDz11;

  return res;
}

il::Array2D<double> StressTensorDueToDDz12P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz12 = 0., sigmayyDz12 = 0., sigmazzDz12 = 0., sigmaxyDz12 = 0.,
  sigmaxzDz12 = 0., sigmayzDz12 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz12 = StressComponentsDueToDDz12P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz12;
  res(0, 1) = sigmaxyDz12;
  res(0, 2) = sigmaxzDz12;
  res(1, 0) = sigmaxyDz12;
  res(1, 1) = sigmayyDz12;
  res(1, 2) = sigmayzDz12;
  res(2, 0) = sigmaxzDz12;
  res(2, 1) = sigmayzDz12;
  res(2, 2) = sigmazzDz12;

  return res;
}

il::Array2D<double> StressTensorDueToDDz13P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz13 = 0., sigmayyDz13 = 0., sigmazzDz13 = 0., sigmaxyDz13 = 0.,
  sigmaxzDz13 = 0., sigmayzDz13 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz13 = StressComponentsDueToDDz13P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz13;
  res(0, 1) = sigmaxyDz13;
  res(0, 2) = sigmaxzDz13;
  res(1, 0) = sigmaxyDz13;
  res(1, 1) = sigmayyDz13;
  res(1, 2) = sigmayzDz13;
  res(2, 0) = sigmaxzDz13;
  res(2, 1) = sigmayzDz13;
  res(2, 2) = sigmazzDz13;

  return res;
}

il::Array2D<double> StressTensorDueToDDz21P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz21 = 0., sigmayyDz21 = 0., sigmazzDz21 = 0., sigmaxyDz21 = 0.,
  sigmaxzDz21 = 0., sigmayzDz21 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz21 = StressComponentsDueToDDz21P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz21;
  res(0, 1) = sigmaxyDz21;
  res(0, 2) = sigmaxzDz21;
  res(1, 0) = sigmaxyDz21;
  res(1, 1) = sigmayyDz21;
  res(1, 2) = sigmayzDz21;
  res(2, 0) = sigmaxzDz21;
  res(2, 1) = sigmayzDz21;
  res(2, 2) = sigmazzDz21;

  return res;
}

il::Array2D<double> StressTensorDueToDDz22P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz22 = 0., sigmayyDz22 = 0., sigmazzDz22 = 0., sigmaxyDz22 = 0.,
  sigmaxzDz22 = 0., sigmayzDz22 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz22 = StressComponentsDueToDDz22P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz22;
  res(0, 1) = sigmaxyDz22;
  res(0, 2) = sigmaxzDz22;
  res(1, 0) = sigmaxyDz22;
  res(1, 1) = sigmayyDz22;
  res(1, 2) = sigmayzDz22;
  res(2, 0) = sigmaxzDz22;
  res(2, 1) = sigmayzDz22;
  res(2, 2) = sigmazzDz22;

  return res;
}

il::Array2D<double> StressTensorDueToDDz23P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz23 = 0., sigmayyDz23 = 0., sigmazzDz23 = 0., sigmaxyDz23 = 0.,
  sigmaxzDz23 = 0., sigmayzDz23 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz23 = StressComponentsDueToDDz23P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz23;
  res(0, 1) = sigmaxyDz23;
  res(0, 2) = sigmaxzDz23;
  res(1, 0) = sigmaxyDz23;
  res(1, 1) = sigmayyDz23;
  res(1, 2) = sigmayzDz23;
  res(2, 0) = sigmaxzDz23;
  res(2, 1) = sigmayzDz23;
  res(2, 2) = sigmazzDz23;

  return res;
}

il::Array2D<double> StressTensorDueToDDz31P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz31 = 0., sigmayyDz31 = 0., sigmazzDz31 = 0., sigmaxyDz31 = 0.,
  sigmaxzDz31 = 0., sigmayzDz31 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz31 = StressComponentsDueToDDz31P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz31;
  res(0, 1) = sigmaxyDz31;
  res(0, 2) = sigmaxzDz31;
  res(1, 0) = sigmaxyDz31;
  res(1, 1) = sigmayyDz31;
  res(1, 2) = sigmayzDz31;
  res(2, 0) = sigmaxzDz31;
  res(2, 1) = sigmayzDz31;
  res(2, 2) = sigmazzDz31;

  return res;
}

il::Array2D<double> StressTensorDueToDDz32P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz32 = 0., sigmayyDz32 = 0., sigmazzDz32 = 0., sigmaxyDz32 = 0.,
  sigmaxzDz32 = 0., sigmayzDz32 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz32 = StressComponentsDueToDDz32P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz32;
  res(0, 1) = sigmaxyDz32;
  res(0, 2) = sigmaxzDz32;
  res(1, 0) = sigmaxyDz32;
  res(1, 1) = sigmayyDz32;
  res(1, 2) = sigmayzDz32;
  res(2, 0) = sigmaxzDz32;
  res(2, 1) = sigmayzDz32;
  res(2, 2) = sigmazzDz32;

  return res;
}

il::Array2D<double> StressTensorDueToDDz33P4(double &a, double &b, double &x1,
                                             double &x2, double &x3, double &Nu,
                                             double &G) {
  double sigmaxxDz33 = 0., sigmayyDz33 = 0., sigmazzDz33 = 0., sigmaxyDz33 = 0.,
  sigmaxzDz33 = 0., sigmayzDz33 = 0.;

  // Planar fault on xy plane
  if (x3 == 0) {
    sigmaxxDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[6];
    sigmayyDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[7];
    sigmazzDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[8];
    sigmaxzDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[9];
    sigmayzDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[10];

  } else {
    sigmaxxDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[0];
    sigmayyDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[1];
    sigmazzDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[2];
    sigmaxzDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[4];
    sigmayzDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[5];
  }

  sigmaxyDz33 = StressComponentsDueToDDz33P4(a, b, x1, x2, x3, Nu, G)[3];

  il::Array2D<double> res{3, 3, 0.};
  res(0, 0) = sigmaxxDz33;
  res(0, 1) = sigmaxyDz33;
  res(0, 2) = sigmaxzDz33;
  res(1, 0) = sigmaxyDz33;
  res(1, 1) = sigmayyDz33;
  res(1, 2) = sigmayzDz33;
  res(2, 0) = sigmaxzDz33;
  res(2, 1) = sigmayzDz33;
  res(2, 2) = sigmazzDz33;

  return res;
}

}  // namespace EQSim