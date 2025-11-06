//
// This file is part of BigWham.
//
// Created by Sylvain Brisson on  05.November.2025
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef BIGWHAM_BIE_MATRIX_GENERATOR_BY_DOF_H
#define BIGWHAM_BIE_MATRIX_GENERATOR_BY_DOF_H

#include <il/core/core.h>

#include "core/bie_kernel.h"
#include "core/elastic_properties.h"
#include "core/mesh.h"
#include "hmat/arrayFunctor/matrix_generator.h"
#include "hmat/hierarchical_representation.h"
#include "bie_matrix_generator.h"

namespace bigwham {

template <typename T>
class BieMatrixGeneratorByDof : public BieMatrixGenerator<T> {
public:

    // Delete the default constructor
    BieMatrixGeneratorByDof() = delete;
    
    // Copy constructor to copy from an object of the base class
    BieMatrixGeneratorByDof(const BieMatrixGenerator<T>& other)
        : BieMatrixGenerator<T>(other) {}

    // Override the set method
    void set(il::int_t b0, il::int_t b1, il::io_t io, il::Array2DEdit<T> M) const override;
}; // class BieMatrixGeneratorByDof


template <typename T>
void BieMatrixGeneratorByDof<T>::set(il::int_t b0, il::int_t b1, il::io_t,
                                       il::Array2DEdit<T> M) const
{
    // IL_EXPECT_MEDIUM(M.size(0) % this->block_size_[0] == 0);
    // IL_EXPECT_MEDIUM(M.size(1) % this->block_size_[1] == 0);
    IL_EXPECT_MEDIUM(b0 + M.size(0) <= this->num_row_points_);
    IL_EXPECT_MEDIUM(b1 + M.size(1) <= this->num_col_points_);

// #pragma omp parallel if (M.size(1) / this->block_size_[1] >= 32)
//   {
// #pragma omp for

    for (il::int_t j1_dof = 0; j1_dof < M.size(1); ++j1_dof) { // Loop over columns = src dofs

        // get the associated element 
        il::int_t k1 = (b1 + j1_dof) / this->block_size_[1]; // rounding via integer division

        // from k1 - permute back to original mesh ordering using permutation
        // of
        // the clusters.
        il::int_t old_k1 = this->hr_->permutation_1_[k1];
        il::int_t e_k1 = this->mesh_src_->GetElementId(old_k1);
        il::int_t is_l = this->mesh_src_->GetElementCollocationId(old_k1);

        auto source_element = this->mesh_src_->GetElement(e_k1);

      // Loop over rows / receiver elementss
      for (il::int_t j0_dof = 0; j0_dof < M.size(0); ++j0_dof) {

            // get the associated element 
            il::int_t k0 = (b0 + j0_dof) / this->block_size_[0];  // rounding via integer division

            il::int_t old_k0 = this->hr_->permutation_0_[k0];
            il::int_t e_k0 = this->mesh_rec_->GetElementId(old_k0); //  receiver element
            il::int_t ir_l = this->mesh_rec_->GetElementCollocationId(old_k0);

            auto receiver_element = this->mesh_rec_->GetElement(e_k0);

            // std::cout << "Influence for elements " << k0 << " * " << k1 << " (dof = " << j0_dof << " * " << j1_dof << "\n";

            // Compute influence
            std::vector<double> st = this->bie_kernel_->influence(
                *source_element, is_l,
                *receiver_element,ir_l
            ); // column major

            // Copy terms to matrix
            IL_EXPECT_FAST(st.size() == this->block_size_[0] * this->block_size_[1]);

            // We only extract one term of what the kernel returns 
            il::int_t j0_dof_local = (b0 + j0_dof) % this->block_size_[0];
            il::int_t j1_dof_local = (b1 + j1_dof) % this->block_size_[1];

            // std::cout << "Extracting (" << j0_dof_local << " * " << j1_dof_local << ") from element pair (" << k0 << " * " << k1 << ")\n";

            M(j0_dof, j1_dof) = st[j1_dof_local * this->block_size_[1] + j0_dof_local];
      }
    }

//   }

}


} // namespace bigwham


#endif // BIGWHAM_BIE_MATRIX_GENERATOR_BY_DOF_H