//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE file for more details.
//

#ifndef BIE_ELASTICPROPERTIES_H
#define BIE_ELASTICPROPERTIES_H

#include <il/core.h>

namespace bigwham {

class ElasticProperties {
private:
  double young_modulus_; // Young's modulus
  double poisson_ratio_; // Poisson's ratio

  double shear_modulus_;             // Lame's second parameter or shear modulus
  double bulk_modulus_;              // Bulk modulus
  double young_modulus_plane_strain_; // Plane-strain Young's modulus
  bool isotropic_;

public:
  ElasticProperties(){};

  ~ElasticProperties(){};

  // Creation of elastic properties class from Young's modulus and Poisson's
  // Ratio
  ElasticProperties(const double E, const double nu) {
    IL_ASSERT(E > 0.);
    IL_ASSERT(nu > -1.);
    IL_ASSERT(nu <= 0.5);

    young_modulus_ = E;
    poisson_ratio_ = nu;
    bulk_modulus_ = young_modulus_ / (3.0 * (1.0 - 2.0 * poisson_ratio_));
    shear_modulus_ = young_modulus_ / (2.0 * (1 + poisson_ratio_));
    young_modulus_plane_strain_ =
        young_modulus_ / (1.0 - poisson_ratio_ * poisson_ratio_);
    isotropic_ = true;
  }

  // Explicit creation from bulk and shear modulus
  void SetBulkShearModulus(double K, double mu) {
    IL_ASSERT(mu > 0.);

    bulk_modulus_ = K;
    shear_modulus_ = mu;
    young_modulus_ = 9.0 * bulk_modulus_ * shear_modulus_ /
                     (3.0 * bulk_modulus_ + shear_modulus_);
    poisson_ratio_ = (3.0 * bulk_modulus_ - 2.0 * shear_modulus_) /
                     (2.0 * (3.0 * bulk_modulus_ + shear_modulus_));
    //  lame1_ = bulk_modulus_ - 2.0 * shear_modulus_ / 3.0;
    young_modulus_plane_strain_ =
        young_modulus_ / (1.0 - poisson_ratio_ * poisson_ratio_);
    isotropic_ = true;
  }

  // Methods for accessing the elastic properties

  double young_modulus() const { return young_modulus_; }
  double poisson_ratio() const { return poisson_ratio_; }
  double bulk_modulus() const { return bulk_modulus_; }

  double shear_modulus() const { return shear_modulus_; }
  double young_modulus_plane_strain() const {
    return young_modulus_plane_strain_;
  }
  bool isotropic() const { return isotropic_; }
};

} // namespace bigwham

#endif // BIE_ELASTICPROPERTIES_H
