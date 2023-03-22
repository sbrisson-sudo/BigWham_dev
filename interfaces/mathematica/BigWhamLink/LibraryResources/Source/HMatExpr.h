//
// Created by Brice Lecampion on 10.11.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Nov. 21 2020


#include <LTemplate.h>

#include<iostream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <sstream>

#include <BigWhamIO.h>

// A class for hierarchical matrices interface to mma
class HMatExpr {
   
    Bigwhamio bie_obj_;

    int dimension_;

    std::string kernel_;

 public:
    
    HMatExpr() {mma::print(" constructor called "); }
    ~HMatExpr() { mma::print(" destructor called ");   }

    // Set the value of the vector
    void set(const mma::RealMatrixRef coor,const mma::IntMatrixRef conn,const char* ker,const mma::RealTensorRef properties,const mint max_leaf_size,
             const double eta,const double eps_aca
             ) {
        //
        std::string kernel = ker;
        kernel_ = kernel;
        dimension_ = coor.cols(); //spatial dimension
        // flatten the tensors ....
        std::vector<double> flat_coor;
        flat_coor.assign(coor.rows()*coor.cols(), 0.);
        for (int i=0;i<coor.rows();i++){
            for (int j=0; j<dimension_;j++){
                flat_coor[dimension_*i+j]=coor(i,j); // row major format
            }
        }
        
        int nvert_elts= conn.cols();
        std::vector<int> flat_con;
        flat_con.assign(conn.rows()*nvert_elts, 0);
        for (int i=0;i<conn.rows();i++){
            for (int j=0; j<nvert_elts;j++){
                flat_con[nvert_elts*i+j]=conn(i,j); // row major format
            }
        }
        std::vector<double> prop;
        prop.assign(properties.begin(),properties.end());
        
        std::cout <<"Now setting up and constructing the HMatrix ...\n";
        mma::print("Now setting up and constructing the HMatrix ...");

        bie_obj_.set(flat_coor,flat_con,kernel,prop,max_leaf_size,eta,eps_aca);
        
        int dim_aux=bie_obj_.getSpatialDimension();
        
        massert(dimension_==dim_aux);
        
        mma::print(" HMatrix is now set !");

    };

    bool isBuilt()  {return bie_obj_.isBuilt(); }
    
    const char * getKernel()  {return kernel_.c_str();}
    
    // get compression
    double getCompressionRatio() {return bie_obj_.getCompressionRatio();}

    // get spatial dimension
    int getSpatialDimension() {return bie_obj_.getSpatialDimension();}

    // get problem dimension (number of unknowns per nodes)
    int getProblemDimension() {return bie_obj_.getProblemDimension();}

    // get matrix size
    mma::IntTensorRef getSize(){
        int nr=bie_obj_.matrixSize(0);
        int nc=bie_obj_.matrixSize(1);
        std::vector<int>  Msize={nr,nc};
        return mma::makeVector<mint>(Msize.size(), Msize.data());
    }
    
    // get permutation of points
    mma::IntTensorRef getPermutation(){
        std::vector<long> permutation=bie_obj_.getPermutation();
        return mma::makeVector<mint>(permutation.size(), permutation.data());
    }
    
    // get collocation points
    mma::RealTensorRef get_collocation_points(){
        std::vector<double> coll_pts=bie_obj_.get_collocation_points();
        int sp_dimension =bie_obj_.getSpatialDimension();
        int npts=coll_pts.size()/sp_dimension;
        return mma::makeMatrix<double>(npts, sp_dimension, coll_pts.data());
    }

    // get hPattern
    mma::IntTensorRef getHpattern(){
          std::vector<long> pattern=bie_obj_.getHpattern();
          int nblocks = pattern.size()/6;
          return mma::makeMatrix<mint>(nblocks, 6, pattern.data());
    }

//   // get hPattern
//    mma::IntTensorRef getHpattern2(){
//      std::vector<int> pattern=bie_obj_.getHpattern2();
//      int nblocks = pattern.size()/6;
//      return mma::makeMatrix<mint>(nblocks, 6, pattern.data());
//    }

    // hdot (in the original non-permutted state)
    mma::RealTensorRef hdot(mma::RealTensorRef xv)   {
        std::vector<double> x;
        x.assign(xv.begin(), xv.end());
        // needs to put here some checks on size to avoid kernel crash
        std::vector<double> y=bie_obj_.matvect(x);
        return mma::makeVector<double>(y.size(), y.data());
    }

//    // hdot tbb recursive (in non-permutted state)
//    mma::RealTensorRef hdotRec(mma::RealTensorRef xv)   {
//      std::vector<double> x;
//      x.assign(xv.begin(), xv.end());
//      // needs to put here some checks on size to avoid kernel crash
//      std::vector<double> y=bie_obj_.hdotProductRec(x);
//      return mma::makeVector<double>(y.size(), y.data());
//    }

//    // compute stresses
//    mma::RealTensorRef computeStresses(
//            const mma::RealTensorRef sol,
//            const mma::RealMatrixRef obsPts,
//            const mint nPts,
//            const mma::RealTensorRef properties,
//            const mma::RealMatrixRef coor,
//            const mma::IntMatrixRef conn,
//            const mbool areDDglobal)   {
//
//        dimension_ = coor.cols(); //spatial dimension
//
//        // solution vector
//        std::vector<double> solu;
//        solu.assign(sol.begin(), sol.end());
//
//        // properties vector
//        std::vector<double> prop;
//        prop.assign(properties.begin(),properties.end());
//
//        // flatten the tensors ....
//
//        // observation points
//        std::vector<double> flat_obsPts;
//        flat_obsPts.assign(obsPts.rows()*obsPts.cols(), 0.);
//        for (int i=0;i<obsPts.rows();i++){
//            for (int j=0; j<dimension_;j++){
//                flat_obsPts[dimension_*i+j]=obsPts(i,j); // row major format
//            }
//        }
//
//        // vertices(nodes) coordinates
//        std::vector<double> flat_coor;
//        flat_coor.assign(coor.rows()*coor.cols(), 0.);
//        for (int i=0;i<coor.rows();i++){
//            for (int j=0; j<dimension_;j++){
//                flat_coor[dimension_*i+j]=coor(i,j); // row major format
//            }
//        }
//
//        // connectivity
//        int nvert_elts= conn.cols();
//        std::vector<int> flat_con;
//        flat_con.assign(conn.rows()*nvert_elts, 0);
//        for (int i=0;i<conn.rows();i++){
//            for (int j=0; j<nvert_elts;j++){
//                flat_con[nvert_elts*i+j]=conn(i,j); // row major format
//            }
//        }
//
//        mma::print(" calling compute stress function ... ");
//
//        // calling function
//        std::vector<double> stressTensor = bie_obj_.computeStresses(solu,flat_obsPts,nPts,prop,flat_coor,flat_con,areDDglobal);
//
//        mma::print(" finished! ");
//
//        int nCompStress = stressTensor.size()/nPts;
//
//        return mma::makeMatrix<double>(nPts, nCompStress, stressTensor.data());
//
//    }
//
//
//    // compute displacements
//    mma::RealTensorRef computeDisplacements(
//            const mma::RealTensorRef sol,
//            const mma::RealMatrixRef obsPts,
//            const mint nPts,
//            const mma::RealTensorRef properties,
//            const mma::RealMatrixRef coor,
//            const mma::IntMatrixRef conn,
//            const mbool areDDglobal)   {
//
//        dimension_ = coor.cols(); //spatial dimension
//
//        // solution vector
//        std::vector<double> solu;
//        solu.assign(sol.begin(), sol.end());
//
//        // properties vector
//        std::vector<double> prop;
//        prop.assign(properties.begin(),properties.end());
//
//        // flatten the tensors ....
//
//        // observation points
//        std::vector<double> flat_obsPts;
//        flat_obsPts.assign(obsPts.rows()*obsPts.cols(), 0.);
//        for (int i=0;i<obsPts.rows();i++){
//            for (int j=0; j<dimension_;j++){
//                flat_obsPts[dimension_*i+j]=obsPts(i,j); // row major format
//            }
//        }
//
//        // vertices(nodes) coordinates
//        std::vector<double> flat_coor;
//        flat_coor.assign(coor.rows()*coor.cols(), 0.);
//        for (int i=0;i<coor.rows();i++){
//            for (int j=0; j<dimension_;j++){
//                flat_coor[dimension_*i+j]=coor(i,j); // row major format
//            }
//        }
//
//        // connectivity
//        int nvert_elts= conn.cols();
//        std::vector<int> flat_con;
//        flat_con.assign(conn.rows()*nvert_elts, 0);
//        for (int i=0;i<conn.rows();i++){
//            for (int j=0; j<nvert_elts;j++){
//                flat_con[nvert_elts*i+j]=conn(i,j); // row major format
//            }
//        }
//
//        mma::print(" calling compute displacement function ... ");
//
//        // calling function
//        std::vector<double> dispVector = bie_obj_.computeDisplacements(solu,flat_obsPts,nPts,prop,flat_coor,flat_con,areDDglobal);
//
//        mma::print(" finished computing displacement ! ");
//
//        int nCompStress = dispVector.size()/nPts;
//
//        return mma::makeMatrix<double>(nPts, nCompStress, dispVector.data());
//
//    }
//
  // get fullBlocks
    mma::SparseMatrixRef<double> getFullBlocks(){
        
        std::vector<int> pos_list;
        std::vector<double> val_list;
        
        std::cout << " calling getFullBlocks \n";
        mma::print("calling getFullBlocks");

        bie_obj_.getFullBlocks(val_list,pos_list);
        std::cout << " n entries: " <<  (val_list.size()) << "\n";
        std::cout << " preparing the Mtensors \n";
        mma::print(" preparing the Mtensors ");
        // be careful makeSparseMatrix require 1-indexing for position !
        for (int i=0;i<pos_list.size();i++){
            pos_list[i]=pos_list[i]+1;
        }
        
        auto vals = mma::makeVector<double>(val_list.size(),val_list.data());
        auto pos  = mma::makeMatrix<mint>(val_list.size(),2,pos_list.data());
        
        std::cout << " preparing the sparse matrix \n";
        mma::print(" preparing the sparse matrix  ");
        int naux = bie_obj_.matrixSize(0);
        
        mma::SparseMatrixRef<double> sm = mma::makeSparseMatrix(pos, vals, naux, naux);
       
        pos.free();
        vals.free();
        
        return sm;
        
    }
    
};
