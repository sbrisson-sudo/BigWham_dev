//
// Created by Alexis Sáez Uribe on 24.01.21.
//

#ifndef BIGWHAMELASTUNITTEST_ELASTICHMATRIX3DT0_H
#define BIGWHAMELASTUNITTEST_ELASTICHMATRIX3DT0_H

#pragma once

#include <Hmat-lib/arrayFunctor/MatrixGenerator.h>
#include <elasticity/3d/Elastic3DT0_element.h>
#include <src/core/ElasticProperties.h>
#include <src/core/FaceData.h>
#include <string>
#include <core/Mesh3D.h>
#include <iostream>

namespace bie {

    template <typename T>
    class ElasticHMatrix3DT0 : public il::MatrixGenerator<T> {
    private:
        il::Array2D<double> point_;
        il::Array<il::int_t> permutation_;
        const bie::Mesh3D mesh_;

    public:
        il::int_t I_want_global_DD;
        il::int_t I_want_global_codomain;
        il::int_t I_want_DD_to_traction_kernel;
        bie::ElasticProperties const elas_;
        ElasticHMatrix3DT0(il::Array2D<double> &point, const il::Array<il::int_t> &permutation,
                           bie::Mesh3D &i_meshtools, bie::ElasticProperties &elas,
                           il::int_t I_want_global_DD, // 0 is local and 1 is global
                           il::int_t I_want_global_codomain, // 0 is local and 1 is global
                           il::int_t I_want_DD_to_traction_kernel); // 0 is displacement kernel and 1 is traction kernel

        il::int_t size(il::int_t d) const override;
        il::int_t blockSize() const override;
        il::int_t sizeAsBlocks(il::int_t d) const override;
        void set(il::int_t b0, il::int_t b1, il::io_t, il::Array2DEdit<T> M) const override;
    };

    template <typename T>
    ElasticHMatrix3DT0<T>::ElasticHMatrix3DT0(il::Array2D<double> &point, const il::Array<il::int_t> &permutation,
                                              bie::Mesh3D &i_meshtools,bie::ElasticProperties &elas, il::int_t I_want_global_DD,
                                              il::int_t I_want_global_codomain,
                                              il::int_t I_want_DD_to_traction_kernel)
            : point_{point},
              permutation_{permutation},
              mesh_{i_meshtools},
              elas_{elas},
              I_want_global_DD{I_want_global_DD},
              I_want_global_codomain{I_want_global_codomain},
              I_want_DD_to_traction_kernel{I_want_DD_to_traction_kernel}
    {
        IL_EXPECT_FAST(point_.size(1) == 3);  // size(1) point==3 in 3D
    };

    template <typename T>
    il::int_t ElasticHMatrix3DT0<T>::size(il::int_t d) const {
        IL_EXPECT_MEDIUM(d == 0 || d == 1);

        return mesh_.numberCollPts() * 3;  // num of nodes * (# of degree of freedom per node)
    }

    template <typename T>
    il::int_t ElasticHMatrix3DT0<T>::blockSize() const {
        return 3;  // # of degree of freedom per node
    }

    template <typename T>
    il::int_t ElasticHMatrix3DT0<T>::sizeAsBlocks(il::int_t d) const {
        IL_EXPECT_MEDIUM(d == 0 || d == 1);

        return (mesh_.numberCollPts());
    }

    template <typename T>
    void ElasticHMatrix3DT0<T>::set(il::int_t b0,
                                    il::int_t b1,
                                    il::io_t,
                                    il::Array2DEdit<T> M) const
    {
        /*  Preliminaries:
         *  Divide the full Hmat in blocks of 3x3 units. 3 is the # of DD for each node of the mesh.
         *  Consider a submatrix M that has certain sizes =< sizes of Hmat. The submatrix M considers a number of couples
         *  node - collocation points thus its dimensions are multiple of 3.
         *  If we think in terms of nodes rather than degree of freedom the natural choice is
         *  to loop over the columns or rows jumping of 3 rows or columns per time. Thus we can define a (block)-column index
         *  and a (block)-row index as to span in such the columns and the rows of Hmat or the matrix M.
         *
         *  The purpose of this routine is to fill M.
         *  Dictionary:
         *  b1 := is the (block)-column number of the top left 3x3 submatrix of M in the numeration of the HMatrix
         *  b0 := is the (block)-row number of the top left 3x3 submatrix of M in the numeration of the HMatrix
         *  j1 := iterator going from 0 to the number of 3x3 submatricies of M in a column - it represents the source node number - it is local
         *  j0 := iterator going from 0 to the number of 3x3 submatricies of M in a row - it represents the receiver node number - it is local
         *  k1 := b1 + j1, is the (block)-column number of the "current" source node expressed in the numeration of the HMatrix
         *  k0 := b1 + j1, is the (block)-column number of the "current" receiver node expressed in the numeration of the HMatrix
         *  e_k1 := taking the floor of k1/(nodes per element) you obtain the element ID
         *  e_k0 := taking the floor of k0/(collocation point per element) you obtain the element ID
         *  -----------
         *  ( the above expression "current" in k0 and k1, refers to the double loop below)
         */

        IL_EXPECT_MEDIUM(M.size(0) % blockSize() == 0);
        IL_EXPECT_MEDIUM(M.size(1) % blockSize() == 0);
        IL_EXPECT_MEDIUM(b0 + M.size(0) / blockSize() <= point_.size(0));
        IL_EXPECT_MEDIUM(b1 + M.size(1) / blockSize() <= point_.size(0));

        //    il::int_t old_k1;
        //    il::int_t old_k0;
        //    il::int_t e_k1, e_k0;
        //    il::Array2D<double> stnl{3,3,0.0};

#ifndef NUMBEROFTHREADS
#define NUMBEROFTHREADS 4
#endif
#pragma omp parallel for num_threads(NUMBEROFTHREADS)
        for (il::int_t j1 = 0; j1 < M.size(1) / blockSize(); ++j1)  { // Loop over a subset of source nodes

            il::int_t old_k1;
            il::int_t old_k0;
            il::int_t e_k1, e_k0;
            il::Array2D<double> stnl{3, 3, 0.0};


            il::int_t k1 = b1 + j1;
            old_k1 = permutation_[k1];
            e_k1 = old_k1; // il::floor(old_k1 / number of nodes per element);

            il::Array2D<double> xv = mesh_.getVerticesElt(e_k1); // get vertices' coordinates of source element
            bie::FaceData elem_data_s(xv, 0); // 0 = interpolation order

            // Loop over a subset of collocation points

            for (il::int_t j0 = 0; j0 < M.size(0) / blockSize(); ++j0)
            {
                il::int_t k0 = b0 + j0;
                old_k0 = permutation_[k0];
                e_k0 = old_k0; // il::floor(old_k1 / number of nodes per element);

                xv = mesh_.getVerticesElt(e_k0); // get vertices' coordinates of receiver element
                bie::FaceData elem_data_r(xv, 0); // 0 = interpolation order

                switch (I_want_DD_to_traction_kernel) {
                    case 0: { //false
                        stnl = bie::NodeDDtriplet_to_CPdisplacement_influence_3DT0(elem_data_s, // source element
                                                                              elem_data_r, // receiver element
                                                                              elas_, // elastic properties
                                                                              I_want_global_DD ,
                                                                              I_want_global_codomain); //https://en.wikipedia.org/wiki/Codomain
                        break;
                    }
                    case 1: { //true
                        stnl = bie::NodeDDtriplet_to_CPtraction_influence_3DT0(elem_data_s,
                                                                          elem_data_r,
                                                                          elas_,
                                                                          I_want_global_DD,
                                                                          I_want_global_codomain); //https://en.wikipedia.org/wiki/Codomain
                        break;
                    }
                    default: { std::cout << "ERROR: bad options given for switch in routine: ElasticHMatrix3DT0 = " << I_want_DD_to_traction_kernel << "\n";}
                }

                for (il::int_t j = 0; j < 3; j++)
                {
                    for (il::int_t i = 0; i < 3; i++)
                    {
                        M(j0 * 3 + i, j1 * 3 + j) = stnl(i, j);
                        // I'm writing on
                        // M( direction , number of DD )
                    }
                }
            }
        }
    }

}

#endif //BIGWHAMELASTUNITTEST_ELASTICHMATRIX3DT0_H
