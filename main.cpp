//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 22.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <iostream>
#include <il/Array2D.h>
#include <src/core/FaceData.cpp>
#include <src/core/FaceData.h>
#include <BigWham.h>



#include <string>
#include <random>



#include <Hmat-lib/cluster/cluster.h>
#include <Hmat-lib/compression/toHMatrix.h>
#include <Hmat-lib/hmatrix/HMatrix.h>
#include <Hmat-lib/hmatrix/HMatrixUtils.h>
#include <Hmat-lib/linearAlgebra/blas/hdot.h>
#include <Hmat-lib/linearAlgebra/factorization/luDecomposition.h>
#include <elasticity/2d/ElasticHMatrix2DP0.h>
#include <elasticity/2d/ElasticHMatrix2DP1.h>
#include <elasticity/PostProcessDDM_2d.h>
#include <elasticity/2d/FullMatrixAssembly2D.h>
#include <src/core/ElasticProperties.h>
//#include <src/solvers/HIterativeSolverUtilities.h>
#include <src/_test/elastic3DR0_element_benchmark.h>

int test2DP1(){

  std::cout << "-------------- test2DP1 ---------------------\n";
  // this routine is for testing different component of the library
  //  here 2DP1 kernel.

// simple 1D mesh
  il::int_t ne=40;

  il::Array2D<double> nodes{ne+1,2,0.};
  il::Array2D<il::int_t> conn{ne,2};
  for (il::int_t i=0;i<ne+1;i++){
    nodes(i,0)=-1.0+2.0*i/ne;
    nodes(i,1)=0.;
  }
 for (il::int_t i=0;i<ne;i++){
   conn(i,0)=i;
   conn(i,1)=i+1;
 }

  bie::Mesh test2(nodes,conn,1);
  il::Array2D<double> coll_points = test2.getCollocationPoints();

  il::int_t leaf_size=32; //test->max_leaf_size;
  il::Timer tt;
  tt.Start();
  const il::Cluster cluster = il::cluster(leaf_size, il::io, coll_points);

  const il::Tree<il::SubHMatrix, 4> hmatrix_tree =
      il::hmatrixTree(coll_points, cluster.partition, 0.0);
  tt.Stop();

  std::cout << "Time for cluster construction " << tt.time() <<"\n";
  tt.Reset();
  std::cout << cluster.partition.depth() <<"\n";

  bie::ElasticProperties elas_aux(1.0,0.2);

  il::HMatrix<double> h_ ;
  const bie::ElasticHMatrix2DP1<double> M{coll_points, cluster.permutation,
                                            test2, elas_aux};

  std::cout << " create h mat " << M.size(0) <<"\n";

  h_= il::toHMatrix(M, hmatrix_tree, 0.001);
  std::cout << " create h mat ended " << h_.isBuilt() <<"\n";
  std::cout << " compression ratio " << il::compressionRatio(h_)<<"\n";

  double Ep=1.0/(1.0-0.2*0.2);
  double sig = 1.0;
  double a=1.0;
  double coef = 4.0 * sig / (Ep);
  // at collocation points
  il::Array<double> wsol_coll{coll_points.size(0), 0.};
  for (int i = 0; i < coll_points.size(0); ++i) {
    if (std::abs(coll_points(i,0)) < a) {
      wsol_coll[i] = coef * sqrt(pow(a, 2) - pow(coll_points(i,0), 2));
    }
  }
  //at corresponding nodes - works here due to simple mesh....
  il::Array<double> wsol_nodes{coll_points.size(0), 0.};
  double aux_x=0.; int j=0;
  for (int i = 0; i < conn.size(0); ++i) {

    aux_x=nodes(conn(i,0),0);
    if (std::abs(aux_x) < a) {
      wsol_nodes[j] = coef * sqrt(pow(a, 2) - pow(aux_x, 2));
    }
    j++;
    aux_x=nodes(conn(i,1),0);
    if (std::abs(aux_x) < a) {
      wsol_nodes[j] = coef * sqrt(pow(a, 2) - pow(aux_x, 2));
    }
    j++;
  }

  il::Array<double> xx{coll_points.size(0)*2};

  for (il::int_t i=0;i<xx.size()/2;i++){
     xx[2*i]=0.0;
     xx[2*i+1]=wsol_nodes[i];
  }

 il::Array< double> y=il::dot(h_,xx);

  // now testing the Bigwhamio class...
  Bigwhamio testbie;
  //Bigwhamio *test = new Bigwhamio();


  std::vector<double> f_coor;
  f_coor.assign(2*test2.numberOfNodes(),0.);
  for (il::int_t i=0;i<test2.numberOfNodes();i++){
    f_coor[2*i]=test2.coordinates(i,0);
    f_coor[2*i+1]=test2.coordinates(i,1);
  }

  std::vector<int64_t> f_conn;
  f_conn.assign(2*test2.numberOfElts(),0);
  for (il::int_t i=0;i<test2.numberOfElts();i++){
    f_conn[2*i]=test2.connectivity(i,0);
    f_conn[2*i+1]=test2.connectivity(i,1);
  }

  std::vector<double> f_prop;
  f_prop.assign(2,0);
  f_prop[0]=elas_aux.getE();
  f_prop[1]=elas_aux.getNu();

  std::string ker="2DP1";

  std::cout << " now setting things in bigwhamio obj \n";

  testbie.set(f_coor,f_conn,ker,f_prop,leaf_size,0.,0.001);


  std::cout        <<  " C R :"<< testbie.getCompressionRatio() <<"\n";

  std::vector<double> x;
  x.assign(xx.size(),0.);
  for (il::int_t i=0;i<xx.size();i++){
    x[i]=xx[i];
  }

  std::vector<double> y2=testbie.hdotProduct(x);
  std::cout << " elastic dot product solution:: \n";
  for (il::int_t i=0;i<coll_points.size(0);i++){
    std::cout <<" y " << i  <<" " << y[2*i] << " - with obj - " << y2[2*i]  <<"\n";
    std::cout <<" y " << i  <<" " << y[2*i+1] << " - with obj - " << y2[2*i+1] <<"\n";
  }

  std::cout << "-------------- End of 2DP1 - tests ------------- " << "\n";

  return  0;

}
///////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
int testS3DP0(){
// testing the S3DP0 kernel implementation on griffith crack
  std::cout << "-------------- testS3DP0 ---------------------\n";

  std::cout << " testing the S3DP0 kernel implementation on griffith crack" << "\n";
// simple 1D mesh
  il::int_t ne=140;

  il::Array2D<double> nodes{ne+1,2,0.};
  il::Array2D<il::int_t> conn{ne,2};

  for (il::int_t i=0;i<ne+1;i++){
    nodes(i,0)=-1.0+2.0*i/ne;
    nodes(i,1)=0.;
  }

  for (il::int_t i=0;i<ne;i++){
    conn(i,0)=i;
    conn(i,1)=i+1;
  }
  bie::Mesh Mesh0(nodes,conn,0);

  il::Array2D<double> coll_points = Mesh0.getCollocationPoints();

  il::int_t leaf_size=12;
  il::Timer tt;
  tt.Start();
  const il::Cluster cluster = il::cluster(leaf_size, il::io, coll_points);

  const il::Tree<il::SubHMatrix, 4> hmatrix_tree =
      il::hmatrixTree(coll_points, cluster.partition, 10.0);
  tt.Stop();
  std::cout << "Time for cluster construction " << tt.time() <<"\n";
  tt.Reset();
  std::cout << cluster.partition.depth() <<"\n";

  bie::ElasticProperties elas_aux(1.0,0.2);

  il::HMatrix<double> h_ ;
  const bie::ElasticHMatrix2DP0<double> M{coll_points, cluster.permutation,
                                          Mesh0, elas_aux, 10000.};

//  std::cout << " in set h mat 2 " << cluster.permutation.size() << " e aca " << test->epsilon_aca <<"\n";

  std::cout << " create h mat " << M.size(0) <<"\n";
  h_= il::toHMatrix(M, hmatrix_tree, 0.001);
  std::cout << " create h mat ended " << h_.isBuilt() <<"\n";
  std::cout << " compression ratio " << il::compressionRatio(h_)<<"\n";

  double Ep=1.0/(1.0-0.2*0.2);
  double sig = 1.0;
  double a=1.0;

  double coef = 4.0 * sig / (Ep);

  // compute analyticla solution at collocation points
  il::Array<double> wsol_coll{coll_points.size(0), 0.};

  for (int i = 0; i < coll_points.size(0); ++i) {
    if (std::abs(coll_points(i,0)) < a) {
      wsol_coll[i] = coef * sqrt(pow(a, 2) - pow(coll_points(i,0), 2));
    }
  }
  //at corresponding nodes - works here due to simple mesh....

  il::Array<double> xx{coll_points.size(0)*2};

  for (il::int_t i=0;i<xx.size()/2;i++){
    xx[2*i]=0.0;
    xx[2*i+1]=wsol_coll[i];
  }

  il::Array< double> y=il::dot(h_,xx);

  std::cout << " Now creating similar BigWham instance \n";
  std::string kernelname = "S3DP0";

  // now testing the Bigwhamio class...
  Bigwhamio testbie;
  //Bigwhamio *test = new Bigwhamio();


  std::vector<double> f_coor;
  f_coor.assign(2*Mesh0.numberOfNodes(),0.);
  for (il::int_t i=0;i<Mesh0.numberOfNodes();i++){
    f_coor[2*i]=Mesh0.coordinates(i,0);
    f_coor[2*i+1]=Mesh0.coordinates(i,1);
  }

  std::vector<int64_t> f_conn;
  f_conn.assign(2*Mesh0.numberOfElts(),0);
  for (il::int_t i=0;i<Mesh0.numberOfElts();i++){
    f_conn[2*i]=Mesh0.connectivity(i,0);
    f_conn[2*i+1]=Mesh0.connectivity(i,1);
  }

  std::vector<double> f_prop;
  f_prop.assign(3,0);
  f_prop[0]=elas_aux.getE();
  f_prop[1]=elas_aux.getNu();
  f_prop[2]=1000.;

  std::cout << " now setting things in bigwhamio obj \n";

  testbie.set(f_coor,f_conn,kernelname,f_prop,22,0.4,0.001);

  std::cout        <<  " C R :"<< testbie.getCompressionRatio() <<"\n";

  std::vector<double> x;
  x.assign(xx.size(),0.);
  for (il::int_t i=0;i<xx.size();i++){
    x[i]=xx[i];
  }

  std::vector<double> y2=testbie.hdotProduct(x);
  std::cout << " elastic dot product solution:: \n";
  for (il::int_t i=0;i<coll_points.size(0);i++){
    std::cout <<" y " << i  <<" " << y[2*i] << " - with obj - " << y2[2*i]  <<"\n";
    std::cout <<" y " << i  <<" " << y[2*i+1] << " - with obj - " << y2[2*i+1] <<"\n";
  }


  std::cout << " build ?" << testbie.isBuilt()
            <<  " CR :"<< testbie.getCompressionRatio() << " CR   " <<"\n";


  std::cout << " TEst h - pattern \n";

  std::string hpatfilename = "/Users/bricelecampion/ClionProjects/BigWham/cmake-build-debug/Hpat.csv";

  il::output_hmatPatternF(h_,hpatfilename);

  il::Array2D<il::int_t> pat_SPOT = il::output_hmatPattern(h_);

  std::cout << "n blocks "  << " / " << pat_SPOT.size(1) << "\n";

  int k=0;
  for (il::int_t i=0;i<pat_SPOT.size(1);i++){
    il::spot_t s(pat_SPOT(0,i));
    k=2;
    if (h_.isFullRank(s)) {
      k=0;
      il::Array2DView<double > A =h_.asFullRank(s);
       std::cout << "block :" << i  << " | " << pat_SPOT(1,i) << "," << pat_SPOT(2,i) <<
       "/ " << pat_SPOT(1,i)+A.size(0)-1 << "," <<pat_SPOT(2,i)+A.size(1)-1 <<
       " -   "  << k << " - "  << A.size(0)*A.size(1) << "\n";
    } else if(h_.isLowRank(s)) {
      k=1;
      il::Array2DView<double > A = h_.asLowRankA(s);
      il::Array2DView<double> B = h_.asLowRankB(s);
      std::cout << "block :" << i  << " | " << pat_SPOT(1,i) << "," << pat_SPOT(2,i) <<
                "/ " << pat_SPOT(1,i)+A.size(0)-1 << "," <<pat_SPOT(2,i)+B.size(0)-1
                << " - "  <<  k <<  " - "  << A.size(0)*A.size(1)+B.size(0)*B.size(1) << "\n";
    }
      std::cout <<" block " << i << " :: " << pat_SPOT(0,i) << " - " << h_.isHierarchical(s) <<"\n";
  }

  std::vector<double> val_list;
  std::vector<int> pos_list;

  testbie.getFullBlocks(val_list, pos_list);

  std::cout << "n full block entry" << val_list.size() <<"\n";
  std::cout << " hmat size " << testbie.matrixSize(0) <<"\n";

  il::Array2D<double> mypts{3,2,0.};
  mypts(0,0)=1.5;
  mypts(1,0)=2.;
  mypts(2,0)=2.5;

  il::Array2D<double> stress=bie::computeStresses2D(mypts,Mesh0,elas_aux,xx,bie::point_stress_s3d_dp0_dd,1000.);
  for (il::int_t i=0;i<3;i++){
    std::cout << "stresses " << stress(i,0) << " - " << stress(i,1) <<  " - " << stress(i,2) <<"\n";

  }
  std::cout << "----------end of S3DP0 test ---------------------\n";

  return 0;
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
int testFullMat()
{
  std::cout << "--------------- testFullMat ---------------------\n";
// simple 1D mesh
  il::int_t ne=100;

  il::Array2D<double> nodes{ne+1,2,0.};
  il::Array2D<il::int_t> conn{ne,2};
  for (il::int_t i=0;i<ne+1;i++){
    nodes(i,0)=-1.0+2.0*i/ne;
    nodes(i,1)=0.;
  }
  for (il::int_t i=0;i<ne;i++){
    conn(i,0)=i;
    conn(i,1)=i+1;
  }

  bie::Mesh mesh(nodes,conn,1);

  bie::ElasticProperties elas(1.0,0.2);

  il::Timer tt;
  tt.Start();
  std::cout << "start of assembly \n " ;
  il::Array2D<double> Mfull = serialFullMatrix2d(mesh,elas,bie::normal_shear_stress_kernel_dp1_dd_nodal,0.);
  tt.Stop();
  std::cout << "time of assembly " << tt.time() <<"\n";
  tt.Reset();

  il::Array<il::int_t> permut{mesh.numberCollocationPoints(),0};
  for (il::int_t i=0;i<mesh.numberCollocationPoints();i++){
    permut[i]=i;
  }
  tt.Start();
  il::Array2D<double> Mfull2 =
      parallelFullMatrix2d(mesh, elas, permut, 0.,
                           bie::normal_shear_stress_kernel_dp1_dd_nodal);
  tt.Stop();
  std::cout << "time of parallel assembly  " << tt.time();

  /// dot-product of linear system and checks
  double Ep=1.0/(1.0-0.2*0.2);
  double sig = 1.0;
  double a=1.0;
  double coef = 4.0 * sig / (Ep);
  // at collocation points
  il::Array<double> wsol_coll{mesh.numberCollocationPoints(), 0.};
  auto coll_points=mesh.getCollocationPoints();
  for (int i = 0; i < coll_points.size(0); ++i) {
    if (std::abs(coll_points(i,0)) < a) {
      wsol_coll[i] = coef * sqrt(pow(a, 2) - pow(coll_points(i,0), 2));
    }
  }
  //at corresponding nodes - works here due to simple mesh....
  il::Array<double> wsol_nodes{coll_points.size(0), 0.};
  double aux_x=0.; int j=0;
  for (int i = 0; i < conn.size(0); ++i) {

    aux_x=nodes(conn(i,0),0);
    if (std::abs(aux_x) < a) {
      wsol_nodes[j] = coef * sqrt(pow(a, 2) - pow(aux_x, 2));
    }
    j++;
    aux_x=nodes(conn(i,1),0);
    if (std::abs(aux_x) < a) {
      wsol_nodes[j] = coef * sqrt(pow(a, 2) - pow(aux_x, 2));
    }
    j++;
  }

  il::Array<double> xx{coll_points.size(0)*2};

  for (il::int_t i=0;i<xx.size()/2;i++){
    xx[2*i]=0.0;
    xx[2*i+1]=wsol_nodes[i];
//    std::cout <<" x_s " << i <<" "<< xx[2*i] <<"\n";
//    std::cout <<" x_n " << i <<" "<< xx[2*i+1] <<"\n";
  }

  il::Array< double> y1=il::dot(Mfull,xx);
  il::Array< double> y2=il::dot(Mfull2,xx);

  for (il::int_t i=0;i<coll_points.size(0);i++){
    std::cout <<" y (1) " << i  <<" / " << y1[2*i] <<" y (2) " << i  << " /  " << y2[2*i] <<"\n";
    std::cout <<" y (1) " << i  <<" / " << y1[2*i+1] <<" y (2) " << i  <<" / " << y2[2*i+1] <<"\n";
  }

  std::cout << "----------end of testFullMat ---------------------\n";

  return 0;

}
///////////////////////////////////////////////////////////////////////////////

int testHdot() {
  std::cout << "--------------- testHdot ---------------------\n";

  // star cracks mesh - crack length unity
  il::int_t nfracs=6;
  il::int_t ne_per_frac=800;
  il::Array<double> rad{ne_per_frac+1,0.};
  il::Array<double> angle{nfracs,0.};
//
  for (il::int_t i=0;i<=ne_per_frac;i++){
    rad[i]=1.*i/ne_per_frac;
  }
  for (il::int_t i=0;i<nfracs;i++){
    angle[i]=2*(il::pi)*i/nfracs;
  }
//
  il::int_t ne=ne_per_frac*nfracs;
  il::Array2D<double> nodes{ne+1,2,0.};
  il::Array2D<il::int_t> conn{ne,2};
  il::int_t index=0;
   nodes(index,0)=0.;
   nodes(index,1)=0.;
  for (il::int_t i=0;i<nfracs;i++){
    for (il::int_t j=1;j<ne_per_frac+1;j++){
      index++;
      nodes(index,0)=rad[j]*cos(angle[i]);
      nodes(index,1)=rad[j]*sin(angle[i]);
    }
  }

  for (il::int_t i=0;i<ne;i++){
    conn(i,0)=i;
    conn(i,1)=i+1;
  }
  for (il::int_t j=0;j<nfracs;j++){
    // first node of each fract is the first nod of the mesh.
    conn(j*ne_per_frac,0)=0;
  }

  std::cout << " N unknowns: " << 4*ne <<"\n";
  bie::Mesh mesh(nodes,conn,1);

  bie::ElasticProperties elas(1.0,0.2);
  il::Array2D<double> coll_points = mesh.getCollocationPoints();

  il::int_t leaf_size=32;
  il::Timer tt;
  tt.Start();
  const il::Cluster cluster = il::cluster(leaf_size, il::io, coll_points);

  tt.Stop();
  std::cout << "Time for  cluster tree construction " << tt.time() <<"\n";
  std::cout << " cluster depth ..." <<  cluster.partition.depth() <<"\n";
  std::cout << " cluster - part " << cluster.permutation.size() << "\n";

  tt.Reset();
//  std::cout << " press enter to continue ...\n";
//  std::cin.ignore(); // pause while user do not enter return

  tt.Start();
  const il::Tree<il::SubHMatrix, 4> hmatrix_tree =
      il::hmatrixTree(coll_points, cluster.partition, 3.0);
  tt.Stop();
  std::cout << "Time for binary cluster tree construction " << tt.time() <<"\n";
  std::cout << " binary cluster depth ..." << hmatrix_tree.depth() << "\n";
  std::cout << " root - " << hmatrix_tree.root().index << "\n";

  tt.Reset();
 // std::cout << " press enter to continue ...\n";
  //std::cin.ignore(); // pause while user do not enter return


  il::HMatrix<double> h_ ;
  const bie::ElasticHMatrix2DP1<double> M{coll_points, cluster.permutation,
                                          mesh, elas};

  std::cout << " create h mat - size " << M.size(0) << " * " << M.size(1) <<"\n";
  tt.Start();
  h_= il::toHMatrix(M, hmatrix_tree, 0.0001);
  tt.Stop();

  std::cout << " create h mat ended " << h_.isBuilt() <<  " in " << tt.time() <<  "\n";
  std::cout << " compression ratio " << il::compressionRatio(h_)<<"\n";
  std::cout << " Compressed memory " << (il::compressionRatio(h_)*8*(4*ne*4*ne)) <<"\n";
  std::cout << " dense case memory " << 8*(4*ne*4*ne) << "\n";

 // std::cout << " press enter to continue ...\n";
  //std::cin.ignore(); // pause while user do not enter return
  tt.Reset();
  tt.Start();
  il::Array2D<il::int_t> pattern=output_hmatPattern(h_);
  tt.Stop();
  std::cout << " time for getting pattern "  << tt.time() <<  " number of blocks "  << pattern.size(1) << "\n";
  tt.Reset();
  /// dot-product of linear system and checks
  double Ep=1.0/(1.0-0.2*0.2);
  double sig = 1.0;
  double a=1.0;
  double coef = 4.0 * sig / (Ep);
  // at collocation points
  il::Array<double> wsol_coll{mesh.numberCollocationPoints(), 0.};
  for (int i = 0; i < coll_points.size(0); ++i) {
    if (std::abs(coll_points(i,0)) < a) {
      wsol_coll[i] = coef * sqrt(pow(a, 2) - pow(coll_points(i,0), 2));
    }
  }
  //at corresponding nodes - works here due to simple mesh....
  il::Array<double> wsol_nodes{coll_points.size(0), 0.};
  double aux_x=0.; int j=0;
  for (int i = 0; i < conn.size(0); ++i) {
    aux_x=nodes(conn(i,0),0);
    if (std::abs(aux_x) < a) {
      wsol_nodes[j] = coef * sqrt(pow(a, 2) - pow(aux_x, 2));
    }
    j++;
    aux_x=nodes(conn(i,1),0);
    if (std::abs(aux_x) < a) {
      wsol_nodes[j] = coef * sqrt(pow(a, 2) - pow(aux_x, 2));
    }
    j++;
  }

  il::Array<double> xx{coll_points.size(0)*2};

  for (il::int_t i=0;i<xx.size()/2;i++){
    xx[2*i]=0.0;
    xx[2*i+1]=wsol_nodes[i];
  }

  tt.Start();
  il::Array< double> y1=il::dot(h_,xx);
  tt.Stop();
  std::cout << " time for recursive hdot " << tt.time() <<"\n";
  tt.Reset();
  tt.Start();
  il::Array< double> y2= il::dotwithpattern_serial(h_, pattern, xx);
  tt.Stop();
  std::cout << " time for Non-recursive hdot " << tt.time() <<"\n";

  il::Array2D<il::int_t> lr_patt{pattern.size(0),0};
  il::Array2D<il::int_t> fr_patt{pattern.size(0),0};
  lr_patt.Reserve(3,pattern.size(1));
  fr_patt.Reserve(3,pattern.size(1));
  il::int_t nfb=0;il::int_t nlb=0;
  for (il::int_t i=0;i<pattern.size(1);i++){
    il::spot_t s(pattern(0,i));
    if (h_.isFullRank(s)){

      fr_patt.Resize(3,nfb+1);
      fr_patt(0,nfb)=pattern(0,i);
      fr_patt(1,nfb)=pattern(1,i);
      fr_patt(2,nfb)=pattern(2,i);
      nfb++;
    } else if (h_.isLowRank(s)){

      lr_patt.Resize(3,nlb+1);
      lr_patt(0,nlb)=pattern(0,i);
      lr_patt(1,nlb)=pattern(1,i);
      lr_patt(2,nlb)=pattern(2,i);
      nlb++;
    } else {
      std::cout <<"error in pattern !\n" ;
      il::abort();
    }
  }

  tt.Reset();
  tt.Start();

  il::Array< double> y3= il::dotwithpattern(h_, fr_patt, lr_patt, xx);
  tt.Stop();
  std::cout << " time for Non-recursive hdot v2 " << tt.time() <<"\n";
  std::cout << " norm y1 " << il::norm(y1,il::Norm::L2) << " y2 " << il::norm(y2,il::Norm::L2)
  <<" y3 " << il::norm(y3,il::Norm::L2) << "\n";

  std::cout << " n FR blocks " << fr_patt.size(1) << "  n LR block " << lr_patt.size(1) <<"\n";

  std::cout << "----------end of test hdot  ---------------------\n";

  return 0;

}

int test3DR0() {
    std::cout << "-------------- test3DR0 ---------------------\n";

    // coordinates
    const std::vector<double> coor={-1.,-1.,0.,
                                    1.,-1.,0.,
                                    1.,1.,0.,
                                    -1.,1.,0.,
                                    -1.,2.,0.,
                                    1.,2.,0.};
    // connectivity
    const std::vector<int64_t> conn = {0,1,2,3,
                                       3,2,5,4};

    const std::vector<double> properties = {100, 0.2}; // Young Modulus , Poisson's ratio
    const int max_leaf_size = 1000;
    const double eta = 0.;
    const double eps_aca = 0.001;

    // create displacement HMAT
    const std::string displacementKernel = "3DR0_displ";
    Bigwhamio displacementHMAT;
    displacementHMAT.set(coor,conn,displacementKernel,properties,
                         max_leaf_size, eta, eps_aca);


    const std::string tractionKernel = "3DR0_traction";
    Bigwhamio tractionHMAT;
    tractionHMAT.set(coor,conn,tractionKernel,properties,
                     max_leaf_size, eta, eps_aca);

    // use the Hdot product
    const std::vector<double>  xx = {1.,2.,3.,4.,5.,6.};
    std::vector<double> x;
    std::cout << "Traction HMAT dot product \n" ;
    x = tractionHMAT.hdotProduct(xx);
    for(int i=0; i<x.size(); ++i) std::cout << x[i] << " ";
    std::cout << "\n" ;
    std::cout << "Displacement HMAT dot product \n" ;
    x = displacementHMAT.hdotProduct(xx);
    for(int i=0; i<x.size(); ++i) std::cout << x[i] << ' ';

    // compute stresses at a set of observation points
    std::vector<double> coorobsp={-10.,-10.,0.,
                                   20.,-20.,0.};
    std::vector<double> mysol = {1.,1.,1.,1.,1.,1.};
    std::vector<double> bb = tractionHMAT.computeStresses(mysol, coorobsp, 2, properties, coor, conn, true);

    for(int i=0; i<bb.size()/6; ++i) {
        std::cout << "\n stress at point #" << i << "\n ";
        for(int j=0; j<6; ++j){
            std::cout << bb[i*6+j] << " ";
        }
    }

    // compute displacements at a set of observation points

    std::vector<double> cc = displacementHMAT.computeDisplacements(mysol, coorobsp, 2, properties, coor, conn, true);

    for(int i=0; i<cc.size()/3; ++i) {
        std::cout << "\n displacement at point #" << i << "\n ";
        for(int j=0; j<3; ++j){
            std::cout << cc[i*3+j] << " ";
        }
    }

    std::cout << "\n----------end of test 3DR0  ---------------------\n";
    return 0;
}


int main() {

  std::cout << "++++++++++++++++++++\n";
  test3DR0();

  //test2DP1();

  //testS3DP0();

//  testFullMat();
//  testHdot();

//  test3DT6Mesh();
//    std::string vertices_file = "/home/alexis/bigwham/vertices5000.csv";
//    std::string connectivity_file = "/home/alexis/bigwham/connectivity5000.csv";
//  test3DT6PennyShaped(vertices_file,connectivity_file);

  std::cout << "\n End of BigWham - exe " << "\n";

}

