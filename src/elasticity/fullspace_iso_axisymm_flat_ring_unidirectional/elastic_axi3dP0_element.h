//
// Created by Alexis Sáez Uribe on 30/05/2022.
// Edited by Ankit on 16 March 2023
//

#ifndef BIGWHAM_ELASTICAXI3DP0_ELEMENT_H
#define BIGWHAM_ELASTICAXI3DP0_ELEMENT_H

#include <vector>

#include <core/elements/Segment.h>
#include <core/ElasticProperties.h>

#include <elasticity/BieElastostatic.h>
#include <elasticity/FsIsoAxiFlatRingUnidirectional/elliptic_integral.h>

using ElemType = bie::Segment<0>;

namespace bie {

class ElasticAxiSymmRingKernel
    : public BieElastostatic<ElemType, ElemType, H> {

public:
  ElasticAxiSymmRingKernel(ElasticProperties &elas, il::int_t dim)
      : BieElastostatic<ElemType, ElemType, H>(elas, dim){};
  std::vector<double> influence(ElemType, il::int_t, ElemType, il::int_t) const;

private:
  double stress_disk_dislocation(double &rObs, double &rSrc) const;
};

std::vector<double>
ElasticAxiSymmRingKernel::influence(ElemType src_el, il::int_t src_colid,
                                    ElemType rec_el,
                                    il::int_t rec_colid) const {

  // get constitutive parameters. NOTE:   Poisson's ratio is taken as zero
  double G = this->elas_.getG();

  auto src_vertices = src_el.getVertices();

  // std::cout << src_vertices.size(0) << "  " <<  src_vertices.size(1) << "\n";

  double rExt = std::sqrt(src_vertices(1, 0) * src_vertices(1, 0) +
                          src_vertices(1, 1) * src_vertices(1, 1));
  double rInt = std::sqrt(src_vertices(0, 0) * src_vertices(0, 0) +
                          src_vertices(0, 1) * src_vertices(0, 1));

  // std::cout << "Vertices good \n";

  // double rExt = source_elt.Xmid(0) + source_elt.size() / 2.0;
  // double rInt = source_elt.Xmid(0) - source_elt.size() / 2.0;

  auto rec_colpts = rec_el.collocation_points();
  double rObs = std::sqrt(rec_colpts(0, 0) * rec_colpts(0, 0) +
                          rec_colpts(0, 1) * rec_colpts(0, 1));

  // std::cout << "Collocation good \n";

  // double rObs = receiver_elt.Xmid(0);

  double IF = this->stress_disk_dislocation(rObs, rExt) -
              stress_disk_dislocation(rObs, rInt);

  il::StaticArray2D<double, 2, 2> traction_vector; // (shear, normal)

  // shear stress
  //  effect of shear dd
  traction_vector(0, 0) = 2.0 * G * IF;
  //  effect of normal dd
  traction_vector(0, 1) = 0.0;

  // normal stress
  //  effect of shear dd
  traction_vector(1, 0) = 0.0;
  //  effect of normal dd
  traction_vector(1, 1) = 2.0 * G * IF; // 2G=E when nu=0

  // std vector output in column major format
  std::vector<double> stnl(4, 0.);
  int k = 0;
  for (int j = 0; j < 2; j++) {
    for (int i = 0; i < 2; i++) {
      stnl[k] = traction_vector(i, j);
      k++;
    }
  }
  return stnl;
}

double ElasticAxiSymmRingKernel::stress_disk_dislocation(
    double &rObs, // radial coordinate of the observation point
    double &rSrc  // radius of the dislocation disk
) const {
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
  } else {
    IF = 0.0;
  }

  return IF;
}

} // namespace bie

#endif // BIGWHAM_ELASTICAXI3DP0_ELEMENT_H
