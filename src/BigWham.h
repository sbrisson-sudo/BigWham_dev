//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 15.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Nov. 12 2020

#include <iostream>

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/Dynamic.h>
#include <il/Map.h>
#include <il/SparseMatrixCSR.h>

#include <Hmat-lib/cluster/cluster.h>
#include <Hmat-lib/compression/toHMatrix.h>
#include <Hmat-lib/hmatrix/HMatrix.h>
#include <Hmat-lib/hmatrix/HMatrixUtils.h>
#include <Hmat-lib/linearAlgebra/blas/hdot.h>

#include <src/core/ElasticProperties.h>
#include <src/core/Mesh2D.h>
#include <src/core/Mesh3D.h>
#include <src/elasticity/PostProcessDDM_2d.h>
#include <src/elasticity/PostProcessDDM_3d.h>

// kernels.
#include <elasticity/2d/ElasticHMatrix2DP0.h>
#include <elasticity/2d/ElasticHMatrix2DP1.h>
#include <elasticity/3d/ElasticHMatrix3DT6.h>
#include <elasticity/3d/ElasticHMatrix3DR0.h>
#include <elasticity/3d/ElasticHMatrix3DT0.h>

//#pragma once

//#include <il/Gmres.h>
//#include <Hmat-lib/linearAlgebra/factorization/luDecomposition.h>
//#include <src/elasticity/jacobi_prec_assembly.h>  // for diagonal
//#include <src/solvers/HIterativeSolverUtilities.h>


class Bigwhamio
        {
          private:
          il::HMatrix<double> h_; // dd to traction
          //   arrays storing the pattern (to speed up the hdot ...).
          il::Array2D<il::int_t> lr_pattern_;  // low rank block pattern
          il::Array2D<il::int_t> fr_pattern_;  // full rank block pattern

          il::Array<il::int_t> permutation_;  // permutation of the dof.

          il::Array2D<double> collocationPoints_;  //  collocation points coordinates

          int dimension_;      // spatial dimension
          int dof_dimension_;  // number of dof per nodes / collocation points

          bool isBuilt_;  // if the class instance is built
                          //  bool isLU_=false;

          // H-matrix parameters
          int max_leaf_size_;
          double eta_;
          double epsilon_aca_;
          //  double epsilon_lu_;
          std::string kernel_;

          public:
          //---------------------------------------------------------------------------

          Bigwhamio(){
              dimension_=0;dof_dimension_=0;
              isBuilt_= false;
              eta_=0.;
              epsilon_aca_=0.001;max_leaf_size_=1;
              kernel_="none";
          };
          ~Bigwhamio() = default;
          void set(const std::vector<double>& coor, const std::vector<int64_t>& conn,
                   const std::string& kernel, const std::vector<double>& properties,
                   const int max_leaf_size,const double eta,const double eps_aca)
                   {
                        // coor and conn are assumed to be passed in row-major storage format
                        kernel_ = kernel;

                        // switch depending on Kernels for mesh building
                        max_leaf_size_ = max_leaf_size;
                        eta_ = eta;
                        epsilon_aca_ = eps_aca;
                        std::cout << " Now setting things... " << kernel_ << "\n";
                        il::Timer tt;

                        // if on kernel name - separating 2D and 3D kernels,
                        // duplicating code for simplicity
                        if ((kernel_ == "2DP1") || (kernel_ == "S3DP0")) {
                          // step 1 - create the mesh object
                          dimension_ = 2;
                          dof_dimension_ = 2;

                          IL_ASSERT(coor.size() % dimension_ == 0);
                          IL_ASSERT(conn.size() % dimension_ == 0);

                          std::cout << "Number of nodes " << coor.size() / dimension_ << " .. mod"
                                    << (coor.size() % dimension_) << "\n";
                          std::cout << " Number of elts " << conn.size() / dimension_ << "\n";

                          il::int_t nvertex = coor.size() / dimension_;
                          il::int_t nelts = conn.size() / dimension_;
                          il::int_t nnodes_elts = dimension_;
                          il::Array2D<double> Coor{nvertex, dimension_,0.}; // columm major order
                          il::Array2D<il::int_t> Conn{nelts, nnodes_elts, 0};

                          // interpolation order
                          int p = 0;
                          if (kernel_ == "2DP1") {
                            p = 1;
                          }
                          std::cout << " interpolation order  " << p << "\n";
                          // populate mesh (loops could be optimized - passage row-major to col-major)
                          int index = 0;
                          for (il::int_t i = 0; i < Coor.size(0); i++) {
                            for (il::int_t j = 0; j < Coor.size(1); j++) {
                              Coor(i, j) = coor[index];
                              index++;
                            }
                          }

                          index = 0;
                          for (il::int_t i = 0; i < Conn.size(0); i++) {
                            for (il::int_t j = 0; j < Conn.size(1); j++) {
                              Conn(i, j) = conn[index];
                              index++;
                            }
                          }

                          bie::Mesh mesh2d(Coor, Conn, p);
                          std::cout << "... mesh done"<< "\n";
                          std::cout << "Number elts " << mesh2d.numberOfElts() <<"\n";

                          collocationPoints_ = mesh2d.getCollocationPoints();
                          std::cout << "Creating cluster tree - number of collocation pts" << collocationPoints_.size(0) <<"\n";

                          tt.Start();
                          const il::Cluster cluster =
                              il::cluster(max_leaf_size_, il::io, collocationPoints_);
                          tt.Stop();
                          std::cout << "Cluster tree creation time :  " << tt.time() << "\n";
                          tt.Reset();

                          permutation_ = cluster.permutation;
                          tt.Start();
                          std::cout << "Creating hmatrix  Tree - \n";
                          il::Tree<il::SubHMatrix, 4> hmatrix_tree =
                              il::hmatrixTree(collocationPoints_, cluster.partition, eta_);
                          tt.Stop();
                          std::cout << "hmatrix  tree creation time :  " << tt.time() << "\n";
                          tt.Reset();

                          // elastic properties
                          std::cout <<" properties vector size " << properties.size() <<"\n";

                          if (kernel_ == "2DP1") {
                            IL_ASSERT(properties.size() == 2);
                          } else if (kernel_ == "S3DP0") {
                            IL_ASSERT(properties.size() == 3);
                          }
                          bie::ElasticProperties elas(properties[0], properties[1]);

                          // now we can populate the h- matrix
                          tt.Start();
                          if (kernel_ == "2DP1")  // or maybe use switch
                          {
                            std::cout << "Kernel Isotropic ELasticity 2D P1 segment \n";
                            const bie::ElasticHMatrix2DP1<double> M{collocationPoints_,
                                                                    permutation_, mesh2d, elas};
                            h_ = il::toHMatrix(M, hmatrix_tree, epsilon_aca_);  //

                          } else if (kernel_ == "S3DP0") {
                            std::cout
                                << "Kernel Isotropic ELasticity Simplified_3D (2D) P0 segment \n";

                            const bie::ElasticHMatrix2DP0<double> M{
                                collocationPoints_, permutation_, mesh2d, elas, properties[2]};
                            h_ = il::toHMatrix(M, hmatrix_tree, epsilon_aca_);  //
                          }
                          tt.Stop();
                          std::cout << "H-mat time = :  " << tt.time() << "\n";
                          std::cout << "H mat set : CR = " << il::compressionRatio(h_)
                                    << " eps_aca " << epsilon_aca_ << " eta " << eta_ << "\n";
                          tt.Reset();
                        } else if (kernel_=="3DT6" || kernel_ == "3DR0_displ" || kernel_ == "3DR0_traction" || kernel_=="3DT0") {
                          // check this  NOTE 1: the 3D mesh uses points and connectivity matrix that are
                          // transposed w.r. to the 2D mesh

                          // step 1 - create the mesh object
                          dimension_ = 3;
                          dof_dimension_ = 3;

                          il::int_t nnodes_elts = 0; // n of nodes per element
                          int p = 0; // interpolation order
                          if (kernel_=="3DT6") {nnodes_elts = 3; p = 2;}
                          else if (kernel_=="3DT0") {nnodes_elts = 3; p = 0;}
                          else if (kernel_ == "3DR0_displ" || kernel_ == "3DR0_traction") {nnodes_elts = 4; p = 0;}
                          else {std::cout << "Invalid kernel name ---\n"; il::abort(); };

                          IL_ASSERT(conn.size() % nnodes_elts == 0);
                          IL_ASSERT(coor.size() % dimension_ == 0);

                          std::cout << " Number of nodes " << coor.size() / dimension_ << " .. mod "
                                    << (coor.size() % dimension_) << "\n";
                          std::cout << " Number of elts " << conn.size() / nnodes_elts << "\n";
                          std::cout << " Interpolation order  " << p << "\n";

                          il::int_t nelts = conn.size() / nnodes_elts;
                          il::int_t nvertex = coor.size() / dimension_;

                          il::Array2D<double> Coor{nvertex, dimension_,0.}; // columm major order
                          il::Array2D<il::int_t> Conn{nelts, nnodes_elts, 0};


                          // populate mesh (loops could be optimized - passage row-major to col-major)
                          int index = 0;
                          for (il::int_t i = 0; i < Coor.size(0); i++) {
                            for (il::int_t j = 0; j < Coor.size(1); j++) {
                              Coor(i, j) = coor[index];
                              index++;
                            }
                          }

                          index = 0;
                          for (il::int_t i = 0; i < Conn.size(0); i++) {
                            for (il::int_t j = 0; j < Conn.size(1); j++) {
                              Conn(i, j) = conn[index];
                              index++;
                            }
                          }

                          bie::Mesh3D mesh3d(Coor, Conn, p);
                          std::cout << "... mesh done"<< "\n";
                          std::cout << " Number elts " << mesh3d.numberOfElts() <<"\n";

                          collocationPoints_ = mesh3d.getCollocationPoints();

                          std::cout << " Coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";
                          std::cout << "Creating cluster tree - number of collocation pts: " << collocationPoints_.size(0) <<"\n";

                          tt.Start();
                          const il::Cluster cluster = il::cluster(max_leaf_size_, il::io, collocationPoints_);
                          tt.Stop();
                          std::cout << "Cluster tree creation time :  " << tt.time() << "\n";
                          tt.Reset();

                          permutation_ = cluster.permutation;
                          tt.Start();
                          std::cout << "Creating hmatrix  Tree - \n";
                          il::Tree<il::SubHMatrix, 4> hmatrix_tree =
                              il::hmatrixTree(collocationPoints_, cluster.partition, eta_);
                          tt.Stop();
                          std::cout << "hmatrix  tree creation time :  " << tt.time() << "\n";
                          tt.Reset();
                          std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";

                          // elastic properties
                          std::cout <<" properties vector size " << properties.size() <<"\n";

                          IL_ASSERT(properties.size() == 2);
                          bie::ElasticProperties elas(properties[0], properties[1]);

                          // now we can populate the h- matrix
                          tt.Start();
                          if (kernel_ == "3DT6")  // or maybe use switch
                          {
                            std::cout << "Kernel Isotropic ELasticity 3D T6 (quadratic) triangle \n";
                            std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";
                            const bie::ElasticHMatrix3DT6<double> M{
                                collocationPoints_, permutation_, mesh3d, elas, 0, 0};
                            h_ = il::toHMatrix(M, hmatrix_tree, epsilon_aca_);  //
                            std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";
                          }
                          else if (kernel_ == "3DT0")
                          {
                            std::cout << "Kernel Isotropic Elasticity 3D T0 (quadratic) triangle \n";
                            std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";
                            const bie::ElasticHMatrix3DT0<double> M{
                                      collocationPoints_, permutation_, mesh3d, elas, 0}; // local_global = 0 if local-local, 1 if global-global
                            h_ = il::toHMatrix(M, hmatrix_tree, epsilon_aca_);  //
                            std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";
                          }
                          else if (kernel_ == "3DR0_displ" || kernel_ == "3DR0_traction")
                          {
                            std::cout << "Kernel Isotropic ELasticity 3D R0 (constant) rectangle \n";
                            std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";

                            int I_want_DD_to_traction_kernel = 999;
                            if ( kernel_ == "3DR0_traction"){
                                // DD to traction HMAT
                                I_want_DD_to_traction_kernel = 1;}
                            else if (kernel_ == "3DR0_displ" ){
                                // DD to displacement HMAT
                                I_want_DD_to_traction_kernel = 0;}
                            const bie::ElasticHMatrix3DR0<double> M{collocationPoints_, permutation_, mesh3d, elas, 0, 0,I_want_DD_to_traction_kernel};
                            h_ = il::toHMatrix(M, hmatrix_tree, epsilon_aca_);
                            std::cout << "HMAT --> built \n";

                            std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";
                            std::cout << "H mat set : CR = " << il::compressionRatio(h_)
                                        << " eps_aca " << epsilon_aca_ << " eta " << eta_ << "\n";
                          }
                          tt.Stop();
                          std::cout << "H-mat time = :  " << tt.time() << "\n";
                          tt.Reset();
                          std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";

                        }
                        else {std::cout << "Invalid kernel name ---\n"; return;
                        };

                        tt.Start();
                        std::cout << "now saving the H-mat patterns...  ";
                        setHpattern();  // set Hpattern
                        tt.Stop();
                        std::cout << " in " << tt.time() << "\n";

                        if (h_.isBuilt()) {
                          isBuilt_ = true;
                        } else {
                          isBuilt_ = false;
                        }

                        std::cout << "end of set() of bigwhamio object \n";
                   }

         //---------------------------------------------------------------------------
          void setHpattern()
          {
                // store the h pattern in a il::Array2D<double> for future use

                IL_EXPECT_FAST(h_.isBuilt());
                il::Array2D<il::int_t> pattern = il::output_hmatPattern(h_);
                //  separate full rank and low rank blocks to speed up the hdot
                il::Array2D<il::int_t> lr_patt{pattern.size(0), 0};
                il::Array2D<il::int_t> fr_patt{pattern.size(0), 0};
                lr_patt.Reserve(3, pattern.size(1));
                fr_patt.Reserve(3, pattern.size(1));
                il::int_t nfb = 0;
                il::int_t nlb = 0;
                for (il::int_t i = 0; i < pattern.size(1); i++) {
                    il::spot_t s(pattern(0, i));
                    if (h_.isFullRank(s)) {
                    fr_patt.Resize(3, nfb + 1);
                    fr_patt(0, nfb) = pattern(0, i);
                    fr_patt(1, nfb) = pattern(1, i);
                    fr_patt(2, nfb) = pattern(2, i);
                    nfb++;
                    } else if (h_.isLowRank(s)) {
                    lr_patt.Resize(3, nlb + 1);
                    lr_patt(0, nlb) = pattern(0, i);
                    lr_patt(1, nlb) = pattern(1, i);
                    lr_patt(2, nlb) = pattern(2, i);
                    nlb++;
                    } else {
                    std::cout << "error in pattern !\n";
                    il::abort();
                    }
                }
                lr_pattern_ = lr_patt;
                fr_pattern_ = fr_patt;
          }

          bool isBuilt() {return isBuilt_;} ;

          //---------------------------------------------------------------------------
          //  get and other methods below
          std::vector<double> getCollocationPoints()
          {
            IL_EXPECT_FAST(isBuilt_);
            std::cout << "beginning of getCollocationPoints bigwham \n";
            std::cout << " spatial dim :" << dimension_ << " collocation dim size :" << collocationPoints_.size(1) << "\n";
            std::cout << " collocation npoints :" << collocationPoints_.size(0) << "\n";
            std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";

            IL_EXPECT_FAST(collocationPoints_.size(1) == dimension_);

            il::int_t npoints = collocationPoints_.size(0);

            std::vector<double> flat_col;
            flat_col.assign(npoints * dimension_, 0.);
            int index = 0;
            for (il::int_t i = 0; i < collocationPoints_.size(0); i++) {
              for (il::int_t j = 0; j < collocationPoints_.size(1); j++) {
                flat_col[index] = collocationPoints_(i, j);
                index++;
              }
            }
            std::cout << "end of getCollocationPoints bigwham \n";
            return flat_col;
          };

          //---------------------------------------------------------------------------
          std::vector<int> getPermutation()
          {
            IL_EXPECT_FAST(isBuilt_);
            std::cout << "coll points dim "<< collocationPoints_.size(0) << " - " << collocationPoints_.size(1) << "\n";

            std::vector<int> permut;
            permut.assign(permutation_.size(), 0);
            for (il::int_t i = 0; i < permutation_.size(); i++) {
              permut[i] = permutation_[i];
            }
            return permut;
          }
          //---------------------------------------------------------------------------

          //---------------------------------------------------------------------------
          double getCompressionRatio()
          {
            IL_EXPECT_FAST(isBuilt_);
              return il::compressionRatio(h_);
          }

          std::string getKernel()  {return  kernel_;}

          int getSpatialDimension() const  {return dimension_;}

          int matrixSize(int k) {return  h_.size(k);};

          //---------------------------------------------------------------------------
             std::vector<int> getHpattern()
          {
            // API function to output the hmatrix pattern
            //  as flattened list via a pointer
            //  the numberofblocks is also returned (by reference)
            //
            //  the pattern matrix is formatted as
            // row = 1 block : i_begin,j_begin, i_end,j_end,FLAG,entry_size
            // with FLAG=0 for full rank and FLAG=1 for low rank
            //
            // we output a flatten row-major order std::vector

            IL_EXPECT_FAST(isBuilt_);

            int fullRankPatternSize1, lowRankPatternSize1;

            lowRankPatternSize1 =  lr_pattern_.size(1);
            fullRankPatternSize1 = fr_pattern_.size(1);

            int numberofblocks = lowRankPatternSize1 + fullRankPatternSize1;
            int len = 6 * numberofblocks;
            std::cout << "number of blocks " << numberofblocks << "\n";

            std::vector<int> patternlist(len,0);

            int index = 0;
            //  starts with full rank

            for (il::int_t j = 0; j < fullRankPatternSize1; j++) {
              il::spot_t s(fr_pattern_(0, j));
              // check is low rank or not
              il::Array2DView<double> A = h_.asFullRank(s);
              //        std::cout << "block :" << i  << " | " << pat_SPOT(1,i) << "," <<
              //        pat_SPOT(2,i) <<
              //                  "/ " << pat_SPOT(1,i)+A.size(0)-1 << ","
              //                  <<pat_SPOT(2,i)+A.size(1)-1 << " -   "  << k << " - "
              //                  << A.size(0)*A.size(1) << "\n";
              patternlist[index++] = fr_pattern_(1, j);
              patternlist[index++] = fr_pattern_(2, j);
              patternlist[index++] = fr_pattern_(1, j) + A.size(0) - 1;
              patternlist[index++] = fr_pattern_(2, j) + A.size(1) - 1;
              patternlist[index++] = 0;
              patternlist[index++] = A.size(0) * A.size(1);
            }

            // then low ranks
                for (il::int_t j = 0; j < lowRankPatternSize1; j++) {
                    il::spot_t s(lr_pattern_(0, j));
                    il::Array2DView<double> A = h_.asLowRankA(s);
                    il::Array2DView<double> B = h_.asLowRankB(s);
                    //        std::cout << "block :" << i  << " | " << pat_SPOT(1,i) << "," <<
                    //        pat_SPOT(2,i) <<
                    //                  "/ " << pat_SPOT(1,i)+A.size(0)-1 << ","
                    //                  <<pat_SPOT(2,i)+B.size(0)-1
                    //                  << " - "  <<  k <<  " - "  <<
                    //                  A.size(0)*A.size(1)+B.size(0)*B.size(1) << "\n";
                    patternlist[index++] = lr_pattern_(1, j);
                    patternlist[index++] = lr_pattern_(2, j);
                    patternlist[index++] = lr_pattern_(1, j) + A.size(0) - 1;
                    patternlist[index++] = lr_pattern_(2, j) + B.size(0) - 1;
                    patternlist[index++] = 1;
                    patternlist[index++] = A.size(0) * A.size(1) + B.size(0) * B.size(1);
                }
            // return a row major flatten vector
            return patternlist;
          }

          //---------------------------------------------------------------------------
          void getFullBlocks(std::vector<double> & val_list,std::vector<int> & pos_list)
          {
            // return the full dense block entries of the hmat as
            // flattened lists
            // val_list(i) = H(pos_list(2*i),pos_list(2*i+1));

            IL_EXPECT_FAST(isBuilt_);

            il::int_t numberofunknowns = h_.size(1);
            il::Array2C<il::int_t> pos{0, 2};
            il::Array<double> val{};
            val.Reserve(numberofunknowns * 4);
            pos.Reserve(numberofunknowns * 4, 2);

            output_hmatFullBlocks(this->h_, val, pos);

            std::cout << "done Full Block: nval " << val.size() << " / " << pos.size(0)
                      << " n^2 " << numberofunknowns * numberofunknowns << "\n";

            IL_EXPECT_FAST((val.size()) == (pos.size(0)));
            IL_EXPECT_FAST(pos.size(1) == 2);

            // outputs
            pos_list.resize((pos.size(0)) * 2);
            val_list.resize(pos.size(0));

            int index = 0;
            for (il::int_t i = 0; i < pos.size(0); i++) {
              pos_list[index++] = pos(i, 0);
              pos_list[index++] = pos(i, 1);
            }
            index = 0;
            for (il::int_t i = 0; i < pos.size(0); i++) {
              val_list[index++] = val[i];
            }

          }


          // ---------------------------------------------------------------------------
          std::vector<double> hdotProduct(const std::vector<double>& x)
          {
            IL_EXPECT_FAST(this->isBuilt_);

                IL_EXPECT_FAST(h_.size(0) == h_.size(1));
                IL_EXPECT_FAST(h_.size(1) == x.size());

                il::Array<double> z{h_.size(1), 0.};

                // permutation of the dofs according to the re-ordering sue to clustering
                il::int_t numberofcollocationpoints = collocationPoints_.size(0);

                for (il::int_t i = 0; i < numberofcollocationpoints; i++) {
                    for (int j = 0; j < dof_dimension_; j++) {
                        z[dof_dimension_ * i + j] = x[dof_dimension_ * (permutation_[i]) + j];
                    }
                }

                z = il::dotwithpattern(h_, fr_pattern_, lr_pattern_, z);
                ////    z = il::dot(h_,z);

                std::vector<double> y;
                y.assign(z.size(), 0.);
                // permut back
                for (il::int_t i = 0; i < numberofcollocationpoints; i++) {
                    for (int j = 0; j < dof_dimension_; j++) {
                        y[dof_dimension_ * (this->permutation_[i]) + j] =
                                z[dof_dimension_ * i + j];
                    }
                }
                return y;
          }
          //---------------------------------------------------------------------------
          std::vector<double> hdotProductInPermutted(const std::vector<double> & x)
          {
                IL_EXPECT_FAST(this->isBuilt_);


                    IL_EXPECT_FAST(h_.size(0) == h_.size(1));
                    IL_EXPECT_FAST(h_.size(1) == x.size());

                    il::Array<double> z{h_.size(1), 0.};
                    for (il::int_t i = 0; i < h_.size(1); i++) {
                        z[i] = x[i];
                    }

                    z = il::dotwithpattern(h_,fr_pattern_,lr_pattern_,z);
                    std::vector<double> y=x;
                    for (il::int_t i = 0; i < h_.size(1); i++) {
                        y[i] = z[i];
                    }
                    return y;
          }

          int getNodesPerElem()
          {        int nnodes_elts =0;
                 if (kernel_=="3DT6" || kernel_ == "3DT0") {nnodes_elts = 3; }
                 else if (kernel_ == "3DR0_displ" || kernel_ == "3DR0_traction") {nnodes_elts = 4; }
                 else if (kernel_ == "2DP1") {nnodes_elts = 2;}
                 else if (kernel_ == "S3DP0") {nnodes_elts = 3;}
                 else {std::cout << "Invalid kernel name ---\n"; il::abort(); };
                 return nnodes_elts;
          }

          int getInterpOrder()
          {     int p =1000;
                if (kernel_=="3DT6") {p = 2;}
                else if (kernel_ == "3DT0") {p = 0;}
                else if (kernel_ == "3DR0_displ" || kernel_ == "3DR0_traction") {p = 0;}
                else if (kernel_ == "2DP1") {p = 1;}
                else if (kernel_ == "S3DP0") {p = 0;}
                else {std::cout << "Invalid kernel name ---\n"; il::abort(); };
                return p;
          }

          il::Array2D<il::int_t>getConn(const std::vector<int64_t>& conn){

                int nnodes_elts = getNodesPerElem();
                il::int_t nelts = conn.size() / nnodes_elts;

                il::Array2D<il::int_t> Conn{nelts, nnodes_elts, 0};
                int index = 0;
                for (il::int_t i = 0; i < Conn.size(0); i++) {
                    for (il::int_t j = 0; j < Conn.size(1); j++) {
                        Conn(i, j) = conn[index];
                        index++;
                    }
                }
                return Conn;
          }

          il::Array2D<double> getCoor(const std::vector<double>& coor){
            il::int_t nvertex = coor.size() / dimension_;
            il::Array2D<double> Coor{nvertex, dimension_, 0.}; // columm major order
            // populate mesh (loops could be optimized - passage row-major to col-major)
            int index = 0;
            for (il::int_t i = 0; i < Coor.size(0); i++) {
                for (il::int_t j = 0; j < Coor.size(1); j++) {
                    Coor(i, j) = coor[index];
                    index++;
                }
            }
            return Coor;
          }

          //---------------------------------------------------------------------------
          std::vector<double> computeStresses(const std::vector<double>& solution,
                                            const std::vector<double>& obsPts,
                                            const int npts, // do we really need it? we know dimension + size of inputs
                                            const std::vector<double>& properties,
                                            const std::vector<double>& coor,
                                            const std::vector<int64_t>& conn,
                                            const bool are_dd_global) {

                 /* BE CAREFUL 2D CASES NEVER TESTED! */

            // PURPOSE: compute stresses at list of points (of size npts )
            //          from a solution vector.
            // INPUT:   "solution" a flattened list containing the solution in terms of DDs
            //          "obsPts" a flattened list containing the observation points coordinates
            //          [ x(1), y(1), z(1), ... ,x(npts), y(npts), z(npts) ]
            //          "npts" the number of points
            //          "properties" is a vector of size 2 or 3 depending on the kernel
            //          It should contain at least
            //          - properties[0] = YoungModulus
            //          - properties[1] = PoissonRatio
            //          - properties[2] = fracture height (mandatory only for "S3DP0" kernel)
            //          "coor" are the coordinates of the vertices of the mesh
            //          "conn" is the connectivity matrix vertices to elements of the mesh
            // OUTPUT:  a flattened list containing the stress at each required point

            std::cout <<" Computing stress tensor ...\n";

            IL_EXPECT_FAST(this->isBuilt_);
            // note solution MUST be of length = number of dofs !

            il::Array2D<double> pts{npts, dimension_};
            il::int_t numberofunknowns = solution.size();
            il::Array<double> solu{numberofunknowns};
            il::Array2D<double> stress;
            bie::ElasticProperties elas(properties[0], properties[1]);

            int p = getInterpOrder(); // interpolation order

            il::Array2D<il::int_t> Conn = getConn(conn);
            il::Array2D<double> Coor = getCoor(coor);
            for (il::int_t i = 0; i < numberofunknowns; i++) {
                solu[i] = solution[i];
            }
            switch (dimension_) {
                case 2: {
                    std::cout << "\n WARNING: computeStresses never tested !!\n";

                    il::int_t index = 0;
                    for (il::int_t i = 0; i < npts; i++) {
                        pts(i, 0) = obsPts[index++];
                        pts(i, 1) = obsPts[index++];
                    }

                    std::cout << " compute stress - " << kernel_ <<"\n";
                    if (kernel_ == "2DP1") {
                        bie::Mesh mesh2d(Coor, Conn, p);
                        stress = bie::computeStresses2D(
                                pts, mesh2d, elas, solu, bie::point_stress_s2d_dp1_dd,0.);
                    } else if (kernel_ == "S3DP0") {
                        bie::Mesh mesh2d(Coor, Conn, p);
                        stress = bie::computeStresses2D(pts, mesh2d, elas, solu, bie::point_stress_s3d_dp0_dd, properties[2]);
                    }
                    break;

                }
                case 3: {
                    /*
                        not implemented yet for 3DT6
                    */
                    il::int_t index = 0;
                    for (il::int_t i = 0; i < npts; i++) {
                        pts(i, 0) = obsPts[index++];
                        pts(i, 1) = obsPts[index++];
                        pts(i, 2) = obsPts[index++];
                    }

                    std::cout << "\n compute stress - " << kernel_ <<"\n";
                    if (kernel_ == "3DT6") {
                        bie::Mesh3D mesh3d(Coor, Conn, p);
                        std::cout << "\n WARNING: not implemented !!\n";
                        il::abort();
                    } else if (kernel_ == "3DT0") {
                        bie::Mesh3D mesh3d(Coor, Conn, p, false);
                        stress = bie::computeStresses3D(pts, mesh3d, elas, solu, bie::point_stress_3DT0, are_dd_global);
                    } else if (kernel_ == "3DR0_traction") {
                        bie::Mesh3D mesh3d(Coor, Conn, p, false);
                        stress = bie::computeStresses3D(pts, mesh3d, elas, solu, bie::point_stress_3DR0, are_dd_global);
                    }
                    break;
                }
            }

            std::vector<double> stress_out(stress.size(0) * stress.size(1));

            il::int_t index = 0;
            for (il::int_t i= 0; i < stress.size(0); i++){
                for (il::int_t j=0; j < stress.size(1); j++){
                    stress_out[index]=stress(i,j);
                    index = index +1 ;
                } // loop on the stress components
            } // loop on the observation points

            return stress_out;
          }


        //---------------------------------------------------------------------------
        std::vector<double> computeDisplacements(std::vector<double>& solution,
                                                std::vector<double>& obsPts,
                                                int npts,
                                                const std::vector<double>& properties,
                                                const std::vector<double>& coor,
                                                const std::vector<int64_t>& conn,
                                                bool are_dd_global) {


            // PURPOSE: compute displacements at list of points (of size npts )
            //          from a solution vector.
            // INPUT:   "solution" a flattened list containing the solution in terms of DDs
            //          "obsPts" a flattened list containing the observation points coordinates
            //          [ x(1), y(1), z(1), ... ,x(npts), y(npts), z(npts) ]
            //          "npts" the number of points
            //          "properties" is a vector of size 2 or 3 depending on the kernel
            //          It should contain at least
            //          - properties[0] = YoungModulus
            //          - properties[1] = PoissonRatio
            //          - properties[2] = fracture height (mandatory only for "S3DP0" kernel)
            //          "coor" are the coordinates of the nodes of the mesh
            //          "conn" is the connectivity matrix node to elements of the mesh
            // OUTPUT:  a flattened list containing the displacements at each required point


            IL_EXPECT_FAST(this->isBuilt_);
            // note solution MUST be of length = number of dofs !

            il::Array2D<double> pts{npts, dimension_};
            il::int_t numberofunknowns = solution.size();
            il::Array<double> solu{numberofunknowns};
            il::Array2D<double> displacements;
            bie::ElasticProperties elas(properties[0], properties[1]);

            int p = getInterpOrder(); // interpolation order

            il::Array2D<il::int_t> Conn = getConn(conn);
            il::Array2D<double> Coor = getCoor(coor);
            for (il::int_t i = 0; i < numberofunknowns; i++) {
                solu[i] = solution[i];
            }
            switch (dimension_) {
                case 2: {

                    il::int_t index = 0;
                    for (il::int_t i = 0; i < npts; i++) {
                        pts(i, 0) = obsPts[index++];
                        pts(i, 1) = obsPts[index++];
                    }

                    std::cout << " compute displacement - " << kernel_ <<"\n";
                    if (kernel_ == "2DP1") {
                        bie::Mesh mesh2d(Coor, Conn, p);
                        std::cout << "\n WARNING: not implemented !!\n";
                        il::abort();
                    } else if (kernel_ == "S3DP0") {
                        bie::Mesh mesh2d(Coor, Conn, p);
                        std::cout << "\n WARNING: not implemented !!\n";
                        il::abort();
                    }
                    break;

                }
                case 3: {
                    /*
                        implemented only for constant DD over a rectangular element
                    */
                    il::int_t index = 0;
                    for (il::int_t i = 0; i < npts; i++) {
                        pts(i, 0) = obsPts[index++];
                        pts(i, 1) = obsPts[index++];
                        pts(i, 2) = obsPts[index++];
                    }

                    std::cout << "\n compute displacement - " << kernel_ <<"\n";
                    if (kernel_ == "3DT6") {
                        bie::Mesh3D mesh3d(Coor, Conn, p);
                        std::cout << "\n WARNING: not implemented !!\n";
                        il::abort();
                    } else if (kernel_ == "3DR0_displ") {
                        bie::Mesh3D mesh3d(Coor, Conn, p, false);
                        displacements = bie::computeDisplacements3D(pts, mesh3d, elas, solu, bie::point_displacement_3DR0, are_dd_global);
                    }
                    break;
                }
            }

            std::vector<double> displacements_out(displacements.size(0) * displacements.size(1));

            il::int_t index = 0;
            for (il::int_t i= 0; i < displacements.size(0); i++){
                for (il::int_t j=0; j < displacements.size(1); j++){
                    displacements_out[index]=displacements(i,j);
                    index = index +1 ;
                } // loop on the displacements components
            } // loop on the observation points

            return displacements_out;
        }

        std::vector<double> getInfluenceCoe(double x, double y, double  z, double  a, double  b, double G, double nu)
        {   std::vector<double> the_stress_out(18);
            il::StaticArray2D<double, 3, 6> Stress;
            Stress = bie::StressesKernelR0(x, y,  z,  a,  b, G, nu);

            il::int_t index = 0;
            for (il::int_t i= 0; i < Stress.size(0); i++){
                for (il::int_t j=0; j < Stress.size(1); j++){
                    the_stress_out[index]=Stress(i,j);
                    index = index +1 ;
                }
            }
            return the_stress_out;
        }
          //

          //  //---------------------------------------------------------------------------
          //  int h_IterativeSolve(double* y, bool permuted, int maxRestart, il::io_t,
          //                       int& nIts, double& rel_error, double* x0) {
          //    // api function for a GMRES iterative solver h_.x=y
          //    //  IMPORTANT: this GMRES  USE the  dot using  dotwithpattern
          //    //
          //    // inputs:
          //    // double * y : pointer to the RHS
          //    // bool permuted : true or false depending if y is already permutted
          //    // il::int_t maxRestart :: maximum number of restart its of the GMRES
          //    // il::io_t :: il flag - variables after are modifed by this function
          //    // nIts :: max number of allowed its of the GMRES, modified as the
          //    actual nb
          //    // of its
          //    // rel_error:: relative error of the GMRES, modifed as the actual
          //    // norm of current residuals
          //    // double * x0 : pointer to the initial guess,
          //    // modified as the actual solution out int :: 0 if status.Ok(),
          //    // output
          //    // int 0 for success, 1 for failure of gmres
          //    const il::int_t rows = h_.size(0);
          //    const il::int_t columns = h_.size(1);
          //    IL_EXPECT_FAST(columns == rows);
          //    IL_EXPECT_FAST(this->isHbuilt == 1);
          //    //
          //
          //    //
          //    il::Array<double> z{columns, 0.}, xini{columns, 0.};
          //
          //    // permut
          //
          //    // note we duplicate 2 vectors here....
          //    // if x is given as permutted or not
          //    if (!permuted) {
          //      for (il::int_t i = 0; i < numberofcollocationpoints; i++) {
          //        for (int j = 0; j < this->dof_dimension; j++) {
          //          z[dof_dimension * i + j] =
          //              y[dof_dimension * (this->permutation_[i]) + j];
          //          xini[dof_dimension * i + j] =
          //              x0[dof_dimension * (this->permutation_[i]) + j];
          //        }
          //      }
          //    } else {
          //      for (int j = 0; j < columns; j++) {
          //        z[j] = y[j];
          //        xini[j] = x0[j];
          //      }
          //    }
          //
          //    // ensure hpattern is set for hdot
          //    if (isHpattern==0){
          //      setHpattern();
          //    }
          //    // prepare call to GMRES
          //    // Instantiation of Matrix type object for GMRES
          //    bie::Matrix<double> A(this->h_,this->fr_pattern_,this->lr_pattern_);
          //    il::Status status{};
          //
          //    // Jacobi preconditioner only so far
          //    bie::DiagPreconditioner<double> Prec(this->diagonal);
          //
          //    il::Gmres<double> gmres{A, Prec, maxRestart};
          //    // solution of the systme via GMRES
          //    gmres.Solve(z.view(), rel_error, nIts, il::io, xini.Edit(), status);
          //    std::cout << "end of Gmres solve, its # " << gmres.nbIterations()
          //              << " norm of residuals: " << gmres.normResidual() << "\n";
          //    // return norm and nits
          //    rel_error = gmres.normResidual();
          //    nIts = gmres.nbIterations();
          //
          //    // return solution in x0
          //    if (!permuted) {
          //      // permut back
          //      for (il::int_t i = 0; i < numberofcollocationpoints; i++) {
          //        for (int j = 0; j < this->dof_dimension; j++) {
          //          x0[dof_dimension * (this->permutation_[i]) + j] =
          //              xini[dof_dimension * i + j];
          //        }
          //      }
          //    } else {
          //      for (int j = 0; j < columns; j++) {
          //        x0[j] = xini[j];
          //      }
          //    }
          //
          //    if (status.Ok()) {
          //      return 0;
          //    } else {
          //      return 1;
          //    }
          //  }

        };  // end class bigwhamio
