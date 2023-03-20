//
// This file is part of HFP.
//
// Created by Brice Lecampion on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIE_ELASTICPROPERTIES_H
#define BIE_ELASTICPROPERTIES_H

#include <il/core.h>

namespace bie {

class ElasticProperties {
 private:
  double young_;    // Young's modulus
  double poiss_;    // Poisson's ratio
  double G_;    // Lame's second parameter or shear modulus
  double bulkm_;    // Bulk modulus
  double youngPS_;  // Plane-strain Young's modulus
  bool iso_;

 public:

  ElasticProperties(){};

  ~ElasticProperties() {};

  // Creation of elastic properties class from Young's modulus and Poisson's
  // Ratio
  ElasticProperties(double YoungModulus, double PoissonRatio) {
    IL_ASSERT(YoungModulus>0.);
    IL_ASSERT(PoissonRatio>-1.);
    IL_ASSERT(PoissonRatio<=0.5);

    young_ = YoungModulus;
    poiss_ = PoissonRatio;
    bulkm_ = young_ / (3.0 * (1.0 - 2.0 * poiss_));
    G_ = young_ / (2.0 * (1 + poiss_));
    youngPS_ = young_ / (1.0 - poiss_ * poiss_);
    iso_=true;

  }

  // Explicit creation from bulk and shear modulus
  void setElasKG(double BulkModulus, double ShearModulus) {
    IL_ASSERT(ShearModulus>0.);

    bulkm_ = BulkModulus;
    G_ = ShearModulus;
    young_ = 9.0 * bulkm_ * G_ / (3.0 * bulkm_ + G_);
    poiss_ = (3.0 * bulkm_ - 2.0 * G_) / (2.0 * (3.0 * bulkm_ + G_));
  //  lame1_ = bulkm_ - 2.0 * G_ / 3.0;
    youngPS_ = young_ / (1.0 - poiss_ * poiss_);
  }

  // Methods for accessing the elastic properties

  double getE() const { return young_; }
  double getNu() const { return poiss_; }
  double getK() const { return bulkm_; }

  double getG() const { return G_; }
  double getEp() const { return youngPS_; }
  double isIsotropic() const {return iso_;}

};


}

#endif  // BIE_ELASTICPROPERTIES_H
