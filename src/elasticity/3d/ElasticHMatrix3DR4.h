//
// Created by Federico Ciardo on 11.08.21.
//

#ifndef INC_3DEQSIM_SRC_ELASTICHMATRIX3D_H
#define INC_3DEQSIM_SRC_ELASTICHMATRIX3D_H

// Inclusion from the project
#include <hmat/arrayFunctor/MatrixGenerator.h>

#include "FullSpaceElasticity.h"
#include <core/Mesh3D.h>
#include "Utils.cpp"

namespace EQSim {

// We do not need a Template class since we are going to use always doubles for
// elastic coefficients
class ElasticHMatrix3D : public il::MatrixGenerator<double> {
 private:
  il::Array2D<double> centroids_;
  il::Array<il::int_t> permutation_;
  EQSim::Mesh mesh_;
  EQSim::SolidMatrixProperties matrixProperties_;
  il::Array2D<il::int_t> neigh_elts_;

 public:
  // Constructor
  ElasticHMatrix3D(il::Array2D<double> &points,
                   const il::Array<il::int_t> &permutation, EQSim::Mesh &mesh,
                   il::Array2D<il::int_t> &neigh_elts,
                   EQSim::SolidMatrixProperties &matrix_prop) {
    centroids_ = points;
    permutation_ = permutation;
    mesh_ = mesh;
    matrixProperties_ = matrix_prop;
    neigh_elts_ = neigh_elts;
    IL_EXPECT_FAST(centroids_.size(1) == 3);
  }

  // Methods to be overrided -> see MatrixGenerator.h
  il::int_t size(il::int_t d) const override {
    IL_EXPECT_MEDIUM(d == 0 || d == 1);

    il::int_t size = mesh_.getNumberOfDofs();

    return size;
  };

  il::int_t blockSize() const override { return 3; };

  il::int_t sizeAsBlocks(il::int_t d) const override {
    IL_EXPECT_MEDIUM(d == 0 || d == 1);

    il::int_t sizeAsBlocks = (mesh_.getNumberOfDofs() / 3);

    return sizeAsBlocks;
  };

  void set(il::int_t b0, il::int_t b1, il::io_t,
           il::Array2DEdit<double> M) const override {
    // Some check conditions
    IL_EXPECT_MEDIUM(M.size(0) % blockSize() == 0);
    IL_EXPECT_MEDIUM(M.size(1) % blockSize() == 0);
    IL_EXPECT_MEDIUM(b0 + M.size(0) / blockSize() <= centroids_.size(0));
    IL_EXPECT_MEDIUM(b1 + M.size(1) / blockSize() <= centroids_.size(0));

    // Get interpolation order
    il::int_t interp_order = mesh_.getInterpolationOrder();

#pragma omp parallel
      {
#pragma omp for collapse(2)
        /// Loop over the source elements
        for (il::int_t k = 0; k < M.size(1) / blockSize(); ++k) {
          /// Loop over the receiver elements
          for (il::int_t k2 = 0; k2 < M.size(0) / blockSize(); ++k2) {
            // Initialization of some variables needed for P4 elements
            il::Array<double> normal_source_elt_i{3, 0.};
            double theta_source_elt_i, Elast_Coeff;
            il::Array<double> relative_distance_receiv_source_i{3, 0.};
            il::Array2D<double> R{3, 3, 0.};   // Rotation matrix
            il::Array2D<double> Rt{3, 3, 0.};  // Transpose of rotation matrix
            il::Array2D<double> ne{9, 3, 0.};
            il::Array2D<double> se1{9, 3, 0.};
            il::Array2D<double> se2{9, 3, 0.};
            il::Array2D<double> Xe{9, 3, 0.};
            il::Array<double> auxXe{3, 0.};
            il::Array<double> auxNe{3, 0.};
            il::Array<double> auxSe1{3, 0.};
            il::Array<double> auxSe2{3, 0.};
            il::Array<double> a_elts_i{9, 0.};
            il::Array<double> b_elts_i{9, 0.};
            il::Array2D<double> neMapped{9, 3, 0.};
            il::Array2D<double> se1Mapped{9, 3, 0.};
            il::Array2D<double> se2Mapped{9, 3, 0.};
            il::Array2D<double> XeMapped{9, 3, 0.};
            il::Array<double> aMapped{9, 0.};
            il::Array<double> bMapped{9, 0.};
            il::Array2D<double> TractionsDueToDDsOnSingleElt{3, 3, 0.};
            EQSim::ElementData elt_s_i;  // Data of source element i
            EQSim::ElementData elt_r;    // Data of receiver element

            // Mapping based on the function I implemented for finding the
            // neighbour elements of a source element and the convention of the
            // nine-element patch of Shou et al.
            il::Array<il::int_t> MapElementPatch{9, 0};
            MapElementPatch[0] = 6;
            MapElementPatch[1] = 7;
            MapElementPatch[2] = 8;
            MapElementPatch[3] = 3;
            MapElementPatch[4] = 4;
            MapElementPatch[5] = 5;
            MapElementPatch[6] = 0;
            MapElementPatch[7] = 1;
            MapElementPatch[8] = 2;

            il::int_t k1 = b1 + k;
            il::int_t elt_k1 = permutation_[k1];  // original source element

            il::int_t k3 = b0 + k2;
            il::int_t elt_k3 = permutation_[k3];  // original receiver element

            // Make sure that we consider only nine-element patches for source &
            // receiver elements
            il::Array<il::int_t> neigh_elts_eltk1 =
                row_selection(neigh_elts_, elt_k1);
            il::Array<il::int_t> neigh_elts_eltk3 =
                row_selection(neigh_elts_, elt_k3);
            il::int_t counter_1 = 0;
            il::int_t counter_2 = 0;
            for (il::int_t I = 0; I < neigh_elts_.size(1); ++I) {
              if (neigh_elts_eltk1[I] >= 0) {
                counter_1 += 1;
              }
              if (neigh_elts_eltk3[I] >= 0) {
                counter_2 += 1;
              }
            }
            if (counter_1 != 9) {
              for (il::int_t j = 0; j < TractionsDueToDDsOnSingleElt.size(0);
                   ++j) {
                for (il::int_t i = 0; i < TractionsDueToDDsOnSingleElt.size(1);
                     ++i) {
                  M(k2 * 3 + i, k * 3 + j) = 0.;
                }
              }
              continue;
            }

            if (counter_2 != 9) {
              for (il::int_t j = 0; j < TractionsDueToDDsOnSingleElt.size(0);
                   ++j) {
                for (il::int_t i = 0; i < TractionsDueToDDsOnSingleElt.size(1);
                     ++i) {
                  M(k2 * 3 + i, k * 3 + j) = 0.;
                }
              }
              continue;
            }

            // Get receiver element data
            elt_r = mesh_.getElementData(elt_k3);

            /// Loop over the nine source elements composing one source element
            /// patch
            for (il::int_t i = 0; i < neigh_elts_.size(1); ++i) {
              // Get source i element data
              il::int_t neigh_inner_elt_i = neigh_elts_eltk1[i];
              elt_s_i = mesh_.getElementData(neigh_inner_elt_i);
              normal_source_elt_i = elt_s_i.getN();
              theta_source_elt_i = elt_s_i.getTheta();

              // Get rotation matrix for the source element i and its transpose
              R = EQSim::RotationMatrix3D(normal_source_elt_i,
                                          theta_source_elt_i);
              Rt = R;
              Rt(0, 1) = R(1, 0);
              Rt(0, 2) = R(2, 0);
              Rt(1, 0) = R(0, 1);
              Rt(1, 2) = R(2, 1);
              Rt(2, 0) = R(0, 2);
              Rt(2, 1) = R(1, 2);

              // Calculate the relative distance receiver - source
              for (il::int_t I = 0;
                   I < relative_distance_receiv_source_i.size(); ++I) {
                relative_distance_receiv_source_i[I] =
                    elt_r.getCentroidElt(I) - elt_s_i.getCentroidElt(I);
              }
              // Switch to frame of source element I
              auxXe = il::dot(Rt, relative_distance_receiv_source_i);
              auxSe1 = il::dot(Rt, elt_r.getS1());
              auxSe2 = il::dot(Rt, elt_r.getS2());
              auxNe = il::dot(Rt, elt_r.getN());

              // Fill 2D array Xe, Se1, Se2, Ne
              for (il::int_t I = 0; I < auxXe.size(); ++I) {
                Xe(i, I) = auxXe[I];
                se1(i, I) = auxSe1[I];
                se2(i, I) = auxSe2[I];
                ne(i, I) = auxNe[I];
              }

              a_elts_i[i] = elt_s_i.getA();
              b_elts_i[i] = elt_s_i.getB();
            }

            // Apply mapping
            for (il::int_t I = 0; I < MapElementPatch.size(); ++I) {
              aMapped[I] = a_elts_i[MapElementPatch[I]];
              bMapped[I] = b_elts_i[MapElementPatch[I]];

              XeMapped(I, 0) = Xe(MapElementPatch[I], 0);
              neMapped(I, 0) = ne(MapElementPatch[I], 0);
              se1Mapped(I, 0) = se1(MapElementPatch[I], 0);
              se2Mapped(I, 0) = se2(MapElementPatch[I], 0);

              XeMapped(I, 1) = Xe(MapElementPatch[I], 1);
              neMapped(I, 1) = ne(MapElementPatch[I], 1);
              se1Mapped(I, 1) = se1(MapElementPatch[I], 1);
              se2Mapped(I, 1) = se2(MapElementPatch[I], 1);

              XeMapped(I, 2) = Xe(MapElementPatch[I], 2);
              neMapped(I, 2) = ne(MapElementPatch[I], 2);
              se1Mapped(I, 2) = se1(MapElementPatch[I], 2);
              se2Mapped(I, 2) = se2(MapElementPatch[I], 2);
            }

            // Call to the elastic kernel
            TractionsDueToDDsOnSingleElt =
                EQSim::TractionsDueToDDsOnSingleEltP4(
                    aMapped, bMapped, XeMapped, se1Mapped, se2Mapped, neMapped,
                    matrixProperties_, il::io);

            Elast_Coeff =
                ((2 * matrixProperties_.getShearModulus()) /
                 (8. * il::pi * (1 - matrixProperties_.getPoissonRatio())));

            // Fill the submatrix
            for (il::int_t j = 0; j < TractionsDueToDDsOnSingleElt.size(0);
                 ++j) {
              for (il::int_t i = 0; i < TractionsDueToDDsOnSingleElt.size(1);
                   ++i) {
                M(k2 * 3 + i, k * 3 + j) =
                    Elast_Coeff * TractionsDueToDDsOnSingleElt(i, j);
              }
            }
          }
        }
      }



  };

  /// Non parallelized but optimized version of set method!
  /*void set(il::int_t b0, il::int_t b1, il::io_t,
           il::Array2DEdit<double> M) const override {
    // Some check conditions
    IL_EXPECT_MEDIUM(M.size(0) % blockSize() == 0);
    IL_EXPECT_MEDIUM(M.size(1) % blockSize() == 0);
    IL_EXPECT_MEDIUM(b0 + M.size(0) / blockSize() <= centroids_.size(0));
    IL_EXPECT_MEDIUM(b1 + M.size(1) / blockSize() <= centroids_.size(0));

    // Get interpolation order
    il::int_t interp_order = mesh_.getInterpolationOrder();

    // Initialization of Traction Matrix, i.e. the matrix of tractions that all
    // the DDs (Dx, Dy, Dz) induce on another element.
    // Column 0 is the effect of Dx
    // Column 1 is the effect of Dy
    // Column 2 is the effect of Dz
    il::Array2D<double> TractionsDueToDDsOnSingleElt{3, 3, 0.};

    // Initialization of some variables
    il::int_t k1, elt_k1;
    il::int_t k3, elt_k3;

    il::Array2D<double> R{3, 3, 0.};   // Rotation matrix
    il::Array2D<double> Rt{3, 3, 0.};  // Transpose of rotation matrix

    EQSim::ElementData elt_r;  // Data of receiver element

    ////// P4 elements //////
    if (interp_order == 4) {
      // Initialization of some variables needed for P4 elements
      EQSim::ElementData elt_s_i;
      il::Array<double> normal_source_elt_i{3, 0.};
      double theta_source_elt_i, Elast_Coeff;
      il::Array<double> relative_distance_receiv_source_i{3, 0.};

      il::Array2D<double> ne{9, 3, 0.};
      il::Array2D<double> se1{9, 3, 0.};
      il::Array2D<double> se2{9, 3, 0.};
      il::Array2D<double> Xe{9, 3, 0.};
      il::Array<double> auxXe{3, 0.};
      il::Array<double> auxNe{3, 0.};
      il::Array<double> auxSe1{3, 0.};
      il::Array<double> auxSe2{3, 0.};
      il::Array<double> a_elts_i{9, 0.};
      il::Array<double> b_elts_i{9, 0.};
      il::Array2D<double> neMapped{9, 3, 0.};
      il::Array2D<double> se1Mapped{9, 3, 0.};
      il::Array2D<double> se2Mapped{9, 3, 0.};
      il::Array2D<double> XeMapped{9, 3, 0.};
      il::Array<double> aMapped{9, 0.};
      il::Array<double> bMapped{9, 0.};

      // Mapping based on the function I implemented for finding the neighbour
      // elements of a source element and the convention of the nine-element
      // patch of Shou et al.
      il::Array<il::int_t> MapElementPatch{9, 0};
      MapElementPatch[0] = 6;
      MapElementPatch[1] = 7;
      MapElementPatch[2] = 8;
      MapElementPatch[3] = 3;
      MapElementPatch[4] = 4;
      MapElementPatch[5] = 5;
      MapElementPatch[6] = 0;
      MapElementPatch[7] = 1;
      MapElementPatch[8] = 2;

      /// Loop over the source elements
      for (il::int_t k = 0; k < M.size(1) / blockSize(); ++k) {
        k1 = b1 + k;
        elt_k1 = permutation_[k1];  // original source element

        /// Loop over the receiver elements
        for (il::int_t k2 = 0; k2 < M.size(0) / blockSize(); ++k2) {
          k3 = b0 + k2;
          elt_k3 = permutation_[k3];  // original receiver element

          // Get receiver element data
          auto elt_k3_mapped = mesh_.getIndexesInnerCentroids(elt_k3);
          elt_r = mesh_.getElementData(elt_k3_mapped);

          /// Loop over the nine source elements composing one source element
          /// patch
          for (il::int_t i = 0; i < neigh_elts_.size(1); ++i) {
            // Get source i element data
            auto elt_k1_mapped = mesh_.getIndexesInnerCentroids(elt_k1);
            il::int_t neigh_inner_elt_i = neigh_elts_(elt_k1_mapped, i);
            elt_s_i = mesh_.getElementData(neigh_inner_elt_i);
            normal_source_elt_i = elt_s_i.getN();
            theta_source_elt_i = elt_s_i.getTheta();

            // Get rotation matrix for the source element i and its transpose
            R = EQSim::RotationMatrix3D(normal_source_elt_i,
                                        theta_source_elt_i);
            Rt = R;
            Rt(0, 1) = R(1, 0);
            Rt(0, 2) = R(2, 0);
            Rt(1, 0) = R(0, 1);
            Rt(1, 2) = R(2, 1);
            Rt(2, 0) = R(0, 2);
            Rt(2, 1) = R(1, 2);

            // Calculate the relative distance receiver - source
            for (il::int_t I = 0; I < relative_distance_receiv_source_i.size();
                 ++I) {
              relative_distance_receiv_source_i[I] =
                  elt_r.getCentroidElt(I) - elt_s_i.getCentroidElt(I);
            }
            // Switch to frame of source element I
            auxXe = il::dot(Rt, relative_distance_receiv_source_i);
            auxSe1 = il::dot(Rt, elt_r.getS1());
            auxSe2 = il::dot(Rt, elt_r.getS2());
            auxNe = il::dot(Rt, elt_r.getN());

            // Fill 2D array Xe, Se1, Se2, Ne
            for (il::int_t I = 0; I < auxXe.size(); ++I) {
              Xe(i, I) = auxXe[I];
              se1(i, I) = auxSe1[I];
              se2(i, I) = auxSe2[I];
              ne(i, I) = auxNe[I];
            }

            a_elts_i[i] = elt_s_i.getA();
            b_elts_i[i] = elt_s_i.getB();
          }

          // Apply mapping
          for (il::int_t I = 0; I < MapElementPatch.size(); ++I) {
            aMapped[I] = a_elts_i[MapElementPatch[I]];
            bMapped[I] = b_elts_i[MapElementPatch[I]];

            XeMapped(I, 0) = Xe(MapElementPatch[I], 0);
            neMapped(I, 0) = ne(MapElementPatch[I], 0);
            se1Mapped(I, 0) = se1(MapElementPatch[I], 0);
            se2Mapped(I, 0) = se2(MapElementPatch[I], 0);

            XeMapped(I, 1) = Xe(MapElementPatch[I], 1);
            neMapped(I, 1) = ne(MapElementPatch[I], 1);
            se1Mapped(I, 1) = se1(MapElementPatch[I], 1);
            se2Mapped(I, 1) = se2(MapElementPatch[I], 1);

            XeMapped(I, 2) = Xe(MapElementPatch[I], 2);
            neMapped(I, 2) = ne(MapElementPatch[I], 2);
            se1Mapped(I, 2) = se1(MapElementPatch[I], 2);
            se2Mapped(I, 2) = se2(MapElementPatch[I], 2);
          }

          // Call to the elastic kernel
          TractionsDueToDDsOnSingleElt = EQSim::TractionsDueToDDsOnSingleEltP4(
              aMapped, bMapped, XeMapped, se1Mapped, se2Mapped, neMapped,
              matrixProperties_, il::io);

          Elast_Coeff =
              ((2 * matrixProperties_.getShearModulus()) /
               (8. * il::pi * (1 - matrixProperties_.getPoissonRatio())));

          // Fill the submatrix
          for (il::int_t j = 0; j < TractionsDueToDDsOnSingleElt.size(0); ++j) {
            for (il::int_t i = 0; i < TractionsDueToDDsOnSingleElt.size(1);
                 ++i) {
              M(k2 * 3 + i, k * 3 + j) =
                  Elast_Coeff * TractionsDueToDDsOnSingleElt(i, j);
            }
          }
        }
      }

      ////// P0 elements //////
    } else if (interp_order == 0) {
      // Initialization of some variables
      EQSim::ElementData elt_s;  // Data of source element
      il::Array<double> normal_source_elt{};
      il::Array<double> normal_receiver_elt{};
      il::Array<double> s1_receiver_elt{};
      il::Array<double> s2_receiver_elt{};
      il::Array<double> ne{};
      il::Array<double> se1{};
      il::Array<double> se2{};
      il::Array<double> relative_distance_receiv_source{3, 0.};
      il::Array<double> Xe{3, 0.};
      double theta_source_elt;

      /// Loop over the source elements
      for (il::int_t k = 0; k < M.size(1) / blockSize(); ++k) {
        k1 = b1 + k;
        elt_k1 = permutation_[k1];  // original source element

        // Get source element data
        elt_s = mesh_.getElementData(elt_k1);

        /// Loop over the receiver elements
        for (il::int_t k2 = 0; k2 < M.size(0) / blockSize(); ++k2) {
          k3 = b0 + k2;
          elt_k3 = permutation_[k3];  // original receiver element

          // Get receiver element data
          elt_r = mesh_.getElementData(elt_k3);

          // Get normal and theta of source element
          normal_source_elt = elt_s.getN();
          theta_source_elt = elt_s.getTheta();

          // Get rotation matrix for the source element and its transpose
          R = EQSim::RotationMatrix3D(normal_source_elt, theta_source_elt);
          Rt = R;
          Rt(0, 1) = R(1, 0);
          Rt(0, 2) = R(2, 0);
          Rt(1, 0) = R(0, 1);
          Rt(1, 2) = R(2, 1);
          Rt(2, 0) = R(0, 2);
          Rt(2, 1) = R(1, 2);

          // Calculate the relative distance receiver - source
          for (il::int_t I = 0; I < relative_distance_receiv_source.size();
               ++I) {
            relative_distance_receiv_source[I] =
                elt_r.getCentroidElt(I) - elt_s.getCentroidElt(I);
          }
          // Switch to frame of source element
          Xe = il::dot(Rt, relative_distance_receiv_source);

          normal_receiver_elt = elt_r.getN();
          s1_receiver_elt = elt_r.getS1();
          s2_receiver_elt = elt_r.getS2();

          se1 = il::dot(Rt, s1_receiver_elt);
          se2 = il::dot(Rt, s2_receiver_elt);
          ne = il::dot(Rt, normal_receiver_elt);

          auto a_elt_s = elt_s.getA();
          auto b_elt_s = elt_s.getB();

          // Call to the elastic kernel
          TractionsDueToDDsOnSingleElt = EQSim::TractionsDueToDDsOnSingleEltP0(
              a_elt_s, b_elt_s, Xe, se1, se2, ne, matrixProperties_, il::io);

          // Fill the submatrix
          for (il::int_t j = 0; j < TractionsDueToDDsOnSingleElt.size(0); ++j) {
            for (il::int_t i = 0; i < TractionsDueToDDsOnSingleElt.size(1);
                 ++i) {
              M(k2 * 3 + i, k * 3 + j) = TractionsDueToDDsOnSingleElt(i, j);
            }
          }
        }
      }

    }
  };*/
};

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_ELASTICHMATRIX3D_H
