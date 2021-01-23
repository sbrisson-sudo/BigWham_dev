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
#include <src/elasticity/3d/Elastic3DT0_element.h>

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

int test3DT0() {

    std::cout << "--------------- test3DT0 ---------------------\n";

    il::StaticArray<double, 3> x;
    x[0] = 0.3;
    x[1] = 0.3;
    x[2] = 0.33;
    il::StaticArray2D<double, 3, 3> xv;
    xv(0,0) = 0.0;
    xv(0,1) = 0.0;
    xv(0,2) = 0.0;
    xv(1,0) = 1.0;
    xv(1,1) = 0.0;
    xv(1,2) = 0.0;
    xv(2,0) = 0.0;
    xv(2,1) = 1.0;
    xv(2,2) = 0.0;

    // get triangle vertices coordinates separated
    il::StaticArray<double, 3> y1, y2, y3;
    for (il::int_t i = 0; i < 3; i++) {
        y1[i] = xv(0, i);
        y2[i] = xv(1, i);
        y3[i] = xv(2, i);
    }
    // subtractions of coordinates vectors
    il::StaticArray<double, 3> y31, y21, y1x;
    for (il::int_t i = 0; i < 3; i++) {
        y31[i] = y3[i] - y1[i];
        y21[i] = y2[i] - y1[i];
        y1x[i] = y1[i] - x[i];
    }

    // local reference system (e1,e2,e3)
    il::StaticArray<double, 3> e1, e2, e3;

    // e1
    e1[0] = y21[0] / il::norm(y21, il::Norm::L2);
    e1[1] = y21[1] / il::norm(y21, il::Norm::L2);
    e1[2] = y21[2] / il::norm(y21, il::Norm::L2);

    // e2
    double sigma = il::dot(y31,y21) / pow(il::norm(y21, il::Norm::L2),2.0);
    il::StaticArray<double, 3> e2vector; // not normalized
    e2vector[0] = y31[0] - sigma * y21[0];
    e2vector[1] = y31[1] - sigma * y21[1];
    e2vector[2] = y31[2] - sigma * y21[2];
    e2[0] = e2vector[0] / il::norm(e2vector, il::Norm::L2);
    e2[1] = e2vector[1] / il::norm(e2vector, il::Norm::L2);
    e2[2] = e2vector[2] / il::norm(e2vector, il::Norm::L2);

    // e3
    e3 = il::cross(e1,e2);

    std::cout << " e1: " << e1[0] << " " << e1[1] << " " << e1[2] <<"\n";
    std::cout << " e2: " << e2[0] << " " << e2[1] << " " << e2[2] <<"\n";
    std::cout << " e3: " << e3[0] << " " << e3[1] << " " << e3[2] <<"\n";

    // geometrical parameters
    double a, b, c;
    a = il::dot(y31,e2);
    b = il::dot(y21,e1);
    c = il::dot(y31,e1);
    double theta1, theta2, theta3;
    theta1 = acos(c / sqrt( pow(c,2.0) + pow(a,2.0) ) );
    theta2 = acos((b - c) / sqrt( pow(b-c,2.0) + pow(a,2.0) ) );
    theta3 = il::pi - (theta1 + theta2);

    std::cout << " a, b, c: " << a << " " << b << " " << c <<"\n";
    std::cout << " theta1, theta2, theta3: " << theta1 << " " << theta2 << " " << theta3 <<"\n";

    // endpoints' coordinates of each triangle edge for the local coordinate systems of each edge

    // vertex 1 coordinates
    double xi1, zeta1, eta;
    xi1 = il::dot(y1x,e1);
    zeta1 = il::dot(y1x,e2);
    eta = il::dot(y1x,e3);

    // auxiliary angles
    double alpha1, alpha2, alpha3;
    alpha1 = 0.;
    alpha2 = il::pi - theta2;
    alpha3 = il::pi + theta1;

    // endpoints' coordinates edge L1
    double p11, p12, q1;
    p11 = xi1;
    p12 = b + xi1;
    q1 = zeta1;

    // endpoints' coordinates edge L2
    double p22, p23, q2;
    p22 = (b + xi1) * cos(alpha2) + zeta1 * sin(alpha2);
    p23 = (c + xi1) * cos(alpha2) + (a + zeta1) * sin(alpha2);
    q2 = -(b + xi1) * sin(alpha2) + zeta1 * cos(alpha2);

    // endpoints' coordinates edge L3
    double p33, p31, q3;
    p33 = (c + xi1) * cos(alpha3) + (a + zeta1) * sin(alpha3);
    p31 = xi1 * cos(alpha3) + zeta1 * sin(alpha3);
    q3 = -xi1 * sin(alpha3) + zeta1 * cos(alpha3);

    // generic recursive functions

    // previous definitions

    // p
    il::StaticArray2D<double, 3, 3> p;
    p(0,0) = p11;
    p(1,1) = p22;
    p(2,2) = p33;
    p(0,1) = p12;
    p(1,2) = p23;
    p(2,0) = p31;

    // q
    il::StaticArray<double, 3> q;
    q[0] = q1;
    q[1] = q2;
    q[2] = q3;

    std::cout << " p11, p12, q1: " << p(0,0) << " " << p(0,1) << " " << q[0] <<"\n";
    std::cout << " p22, p23, q2: " << p(1,1) << " " << p(1,2) << " " << q[1] <<"\n";
    std::cout << " p33, p31, q3: " << p(2,2) << " " << p(2,0) << " " << q[2] <<"\n";

    // rho - distance between source/receiver point 'x' and i-th vertex of triangle
    il::StaticArray<double, 3> rho;
    for (int i = 0; i < 3; i++) {
        rho[i] = sqrt( pow(p(i,i),2.0) + pow(q[i],2.0) + pow(eta,2.0) );
    }

    std::cout << " rho: " << rho[0] << " " << rho[1] << " " << rho[2] <<"\n";

    // definition of cases - where the projected source/receiver point lies
    // also, computation of theta_0

    // rho_plane - distance between projected source/receiver point 'x' and i-th vertex of triangle
    il::StaticArray<double, 3> rho_plane;
    for (int i = 0; i < 3; i++) {
        rho_plane[i] = sqrt( pow(p(i,i),2.0) + pow(q[i],2.0) );
    }

    std::cout << " rho_plane: " << rho_plane[0] << " " << rho_plane[1] << " " << rho_plane[2] <<"\n";

    double eps_tol = 2.221e-016; // 1.0e-015 parameter used for "if conditions" involving inequalities due to numerical precision

    int id;
    double theta_0;
    if( q[0] > eps_tol || q[1] > eps_tol || q[2] > eps_tol ){
        // if projected point lies strictly outside of the triangle
        id = -1;
        theta_0 = 0.0;
    } else if( q[0] < -eps_tol && q[1] < -eps_tol && q[2] < -eps_tol ){
        // if projected point lies strictly inside of the triangle
        id = 0;
        theta_0 = 2*il::pi;
    } else if( rho_plane[0] > eps_tol && rho_plane[1] > eps_tol && rho_plane[2] > eps_tol){
        // if projected point lies on an edge of the triangle but not on a vertex
        id = 4;
        theta_0 = il::pi;
    } else if( rho_plane[0] < eps_tol ){
        // if projected point lies exactly on vertex 1
        id = 1;
        theta_0 = theta1;
    } else if( rho_plane[1] < eps_tol ){
        // if projected point lies exactly on vertex 2
        id = 2;
        theta_0 = theta2;
    } else if( rho_plane[2] < eps_tol ){
        // if projected point lies exactly on vertex 3
        id = 3;
        theta_0 = theta3;
    }

    std::cout << " id: " << id << "\n";
    std::cout << " theta_0 for unit triangle: " << theta_0 << "\n";
    double teta2 = atan(1.0);
    double teta1 = il::pi/2;
    double teta3 = il::pi-teta1-teta2;

    std::cout << " tetas: " << teta1 << " " << teta2 << " " << teta3 <<"\n";

    // check if summation algorithm is working

    double sum = 0.0;
    for (int i = 0; i < 3; i++) {
        sum += rho_plane[i];
    }

    std::cout << " sum: " << sum << "\n";

    // now the generic recursive functions

    int i1; // i+1, will be used in all recursive generic functions

    // d
    il::StaticArray<double, 3> d;
    for (int i = 0; i < 3; i++) {
        d[i] = pow(q[i],2.0) + pow(eta, 2.0);
    }

    // rho tilde
    il::StaticArray<double, 3> rho_t;
    for (int i = 0; i < 3; i++) {
        i1 = i + 1;
        if (i1 > 2) { i1 = 0; } // if index is equal to 3, replace by 0
        rho_t[i] = rho[i] - rho[i1];
    }

    // phi
    il::StaticArray<double, 3> phi;
    for (int i = 0; i < 3; i++) {
        i1 = i + 1;
        if (i1 > 2) { i1 = 0; } // if index is equal to 3, replace by 0
        phi[i] = p(i, i) * rho[i] - p(i, i1) * rho[i1];
    }

    // chi
    // Note: in Nintcheu Fata (2009, 2011) there are no comments for the case
    // in which the log arguments are negative or zero and this can happen! Brice
    // regularized this in his mathematica notebook but I don't know yet how (Alexis).
    // I use here Brice's solution
    il::StaticArray<double, 3> chi;
    for (int i = 0; i < 3; i++) {
        i1 = i + 1;
        if (i1 > 2) { i1 = 0; } // if index is equal to 3, replace by 0
        if ((p(i, i) + rho[i]) > eps_tol && (p(i, i1) + rho[i1]) > eps_tol) {
            // if the log arguments are strictly positive
            chi[i] = log(p(i, i) + rho[i]) - log(p(i, i1) + rho[i1]);
        } else if ((p(i, i) + rho[i]) < -eps_tol && (p(i, i1) + rho[i1]) < -eps_tol) {
            // if the log arguments are strictly negative
            double d = pow(q[i], 2.0) + pow(eta, 2.0);
            double p1 = p(i, i) + rho[i];
            double p2 = p(i, i1) + rho[i1];
            double z1 = (d / p1) / p1;
            double z2 = (d / p2) / p2;
            chi[i] = log(
                    abs(p2 / p1) * ((1.0 - 0.25 * z1 * (1.0 - 0.5 * z1)) / (1.0 - 0.25 * z2 * (1.0 - 0.5 * z2))));
        } else {
            // if any of the log arguments is zero
            chi[i] = 0.0;
        }
    }

    // delta
    il::StaticArray<double, 3> delta;
    for (int i = 0; i < 3; i++) {
        i1 = i + 1;
        if (i1 > 2) { i1 = 0; } // if index is equal to 3, replace by 0
        delta[i] = p(i, i) / rho[i] - p(i, i1) / rho[i1];
    }

    // L
    il::StaticArray<double, 3> L;
    for (int i = 0; i < 3; i++) {
        i1 = i + 1;
        if (i1 > 2) { i1 = 0; } // if index is equal to 3, replace by 0
        L[i] = 1.0 / rho[i] - 1.0 / rho[i1];
    }

    // D
    il::StaticArray<double, 3> D;
    for (int i = 0; i < 3; i++) {
        i1 = i + 1;
        if (i1 > 2) { i1 = 0; } // if index is equal to 3, replace by 0
        D[i] = p(i, i) / pow(rho[i], 3.0) - p(i, i1) / pow(rho[i1], 3.0);
    }

    // Lambda
    il::StaticArray<double, 3> Lambda;
    for (int i = 0; i < 3; i++) {
        i1 = i + 1;
        if (i1 > 2) { i1 = 0; } // if index is equal to 3, replace by 0
        Lambda[i] = 1.0 / pow(rho[i], 3.0) - 1.0 / pow(rho[i1], 3.0);
    }

    // gamma

    // Note: in Nintcheu Fata (2009, 2011) there are no comments for the case
    // in which the arctan arguments are not defined (in this case, 0/0) and this can happen! Brice
    // regularized this in his mathematica notebook by taking the limits, I understand the case of
    // the edge, but not the cases for the vertices. I just coded up Brice's solution for now, but
    // I didn't include the case inside the element and the case right at the vertices
    il::StaticArray<double, 3> gamma;

    if (id == 1) {
        // if the projected source/receiver point lies on vertex 1
        gamma[0] = 0.0;
        gamma[2] = 0.0;
        int i = 1;
        i1 = 2;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
    } else if (id == 2) {
        // if the projected source/receiver point lies on vertex 2
        gamma[0] = 0.0;
        gamma[1] = 0.0;
        int i = 2;
        i1 = 0;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
    } else if (id == 3) {
        // if the projected source/receiver point lies on vertex 3
        gamma[1] = 0.0;
        gamma[2] = 0.0;
        int i = 0;
        i1 = 1;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
    } else if (abs(q[0]) < eps_tol) {
        // if the projected source/receiver point lies on the edge 1 or its prolongation (in both directions)
        // the limit is needed only for the case eta = 0. However, for eta =! 0, the value of gamma1 is the same
        // as the limit (=zero), so that's why I don't make the difference
        gamma[0] = 0.0;
        int i = 1;
        i1 = 2;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
        i = 2;
        i1 = 0;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
    } else if (abs(q[1]) < eps_tol) {
        // if the projected source/receiver point lies on the edge 2 or its prolongation (in both directions)
        // the limit is needed only for the case eta = 0. However, for eta =! 0, the value of gamma2 is the same
        // as the limit (=zero), so that's why I don't make the difference
        gamma[1] = 0.0;
        int i = 0;
        i1 = 1;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
        i = 2;
        i1 = 0;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
    } else if (abs(q[2]) < eps_tol) {
        // if the projected source/receiver point lies on the edge 3 or its prolongation (in both directions)
        // the limit is needed only for the case eta = 0. However, for eta =! 0, the value of gamma3 is the same
        // as the limit (=zero), so that's why I don't make the difference
        gamma[2] = 0.0;
        int i = 0;
        i1 = 1;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
        i = 1;
        i1 = 2;
        gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                        (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                   - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                          (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
    } else {
        // any other case
        for (int i = 0; i < 3; i++) {
            i1 = i + 1;
            if (i1 > 2) { i1 = 0; } // if index is equal to 3, replace by 0
            gamma[i] = atan((-2.0 * p(i, i) * q[i] * eta * rho[i]) /
                            (pow(q[i], 2.0) * pow(rho[i], 2.0) - pow(p(i, i), 2.0) * pow(eta, 2.0)))
                       - atan((-2.0 * p(i, i1) * q[i] * eta * rho[i1]) /
                              (pow(q[i], 2.0) * pow(rho[i1], 2.0) - pow(p(i, i1), 2.0) * pow(eta, 2.0)));
        }
    }

    // theta for displacements and tractions at the boundary of the domain
    double theta;
    double thetaSum = 0.0;
    for (int i = 0; i < 3; i++) {
        thetaSum += gamma[i];
    }
    if (eta < eps_tol) {
        // if eta is negative or zero
        theta = 0.5 * thetaSum - theta_0;
    } else {
        // if eta is strictly positive
        theta = 0.5 * thetaSum + theta_0;
    }

    std::cout << " Checking parameters sensitive to cases... " << "\n";

    std::cout << " theta: " << theta << "\n";
    std::cout << " eta: " << eta << "\n";
    std::cout << " gamma 1: " << gamma[0] << " gamma 2: " << gamma[1] << " gamma 3: " << gamma[2] << "\n";
    std::cout << " theta_0: " << theta_0 << "\n";

    // sin(alpha[i]) and cos(alpha[i])

    // alpha
    il::StaticArray<double, 3> alpha;
    alpha[0] = alpha1;
    alpha[1] = alpha2;
    alpha[2] = alpha3;

    il::StaticArray<double, 3> sinAlpha;
    il::StaticArray<double, 3> cosAlpha;
    for (int i = 0; i < 3; i++) {
        sinAlpha[i] = sin(alpha[i]);
        cosAlpha[i] = cos(alpha[i]);
    }

    // check generic integrals for DD1

    double I5_Xi = bie::i5_Xi(delta,d,sinAlpha);
    double I7_Xi_Xi_Xi = bie::i7_Xi_Xi_Xi(q,d,sinAlpha,cosAlpha,D,Lambda,delta);
    double I7_Xi_Zeta_Zeta = bie::i7_Xi_Zeta_Zeta(q,d,sinAlpha,cosAlpha,D,Lambda,delta);
    double I7_Xi = bie::i7_Xi(d,D,delta,sinAlpha);
    double I7_Xi_Xi_Zeta = bie::i7_Xi_Xi_Zeta(q,d,sinAlpha,cosAlpha,D,Lambda,delta);
    double I5_Zeta = bie::i5_Zeta(delta,d,cosAlpha);
    double I7_Xi_Xi_Aux = bie::i7_Xi_Xi_Aux(eta,cosAlpha,Lambda,sinAlpha,q,d,D,delta);
    double I7_Xi_Xi_By_Aux = I7_Xi_Xi_Aux/pow(eta,2.0)+theta/(15.0*pow(eta,3.0)); // to verify when eta =! 0 against mma
    double I5_Zeta_Zeta_Aux = bie::i5_Zeta_Zeta_Aux(L,sinAlpha,q,d,delta,cosAlpha);
    double I5_Zeta_Zeta_By_Aux = I5_Zeta_Zeta_Aux + theta/(3*eta); // to verify when eta =! 0 against mma
    double I7_Xi_Zeta = bie::i7_Xi_Zeta(sinAlpha,Lambda,cosAlpha,q,d,D,delta);
    double I5_Xi_Zeta = bie::i5_Xi_Zeta(L,sinAlpha,q,d,delta,cosAlpha);

    std::cout << " Checking generic integrals for DD1 ... " << "\n";

    std::cout << " I5_Xi: " << I5_Xi << "\n";
    std::cout << " I7_Xi_Xi_Xi: " << I7_Xi_Xi_Xi << "\n";
    std::cout << " I7_Xi_Zeta_Zeta: " << I7_Xi_Zeta_Zeta << "\n";
    std::cout << " I7_Xi: " << I7_Xi << "\n";
    std::cout << " I7_Xi_Xi_Zeta: " << I7_Xi_Xi_Zeta << "\n";
    std::cout << " I5_Zeta: " << I5_Zeta << "\n";
//    std::cout << " I7_Xi_Xi_Aux: " << I7_Xi_Xi_Aux << "\n";
    std::cout << " I7_Xi_Xi_By_Aux: " << I7_Xi_Xi_By_Aux << "\n";
//    std::cout << " I5_Zeta_Zeta_Aux: " << I5_Zeta_Zeta_Aux << "\n";
    std::cout << " I5_Zeta_Zeta_By_Aux: " << I5_Zeta_Zeta_By_Aux << "\n";
    std::cout << " I7_Xi_Zeta: " << I7_Xi_Zeta << "\n";
    std::cout << " I5_Xi_Zeta: " << I5_Xi_Zeta << "\n";

    // check generic integrals for DD2

    double I7_Zeta_Zeta_Zeta = bie::i7_Zeta_Zeta_Zeta(q,d,sinAlpha,cosAlpha,D,Lambda,delta);
    double I7_Zeta = bie::i7_Zeta(d,D,delta,cosAlpha);
    double I7_Zeta_Zeta_Aux = bie::i7_Zeta_Zeta_Aux(eta,cosAlpha,Lambda,sinAlpha,q,d,D,delta);
    double I7_Zeta_Zeta_By_Aux = I7_Zeta_Zeta_Aux/pow(eta,2.0)+theta/(15.0*pow(eta,3.0)); // to verify when eta =! 0 against mma
    double I5_Xi_Xi_Aux = bie::i5_Xi_Xi_Aux(L,sinAlpha,q,d,delta,cosAlpha);
    double I5_Xi_Xi_By_Aux = I5_Xi_Xi_Aux + theta/(3*eta); // to verify when eta =! 0 against mma

    std::cout << " Checking generic integrals for DD2 ... " << "\n";

    std::cout << " I7_Zeta_Zeta_Zeta: " << I7_Zeta_Zeta_Zeta << "\n";
    std::cout << " I7_Zeta: " << I7_Zeta << "\n";
//    std::cout << " I7_Zeta_Zeta_Aux: " << I7_Zeta_Zeta_Aux << "\n";
    std::cout << " I7_Zeta_Zeta_By_Aux: " << I7_Zeta_Zeta_By_Aux << "\n";
//    std::cout << " I5_Xi_Xi_Aux: " << I5_Xi_Xi_Aux << "\n";
    std::cout << " I5_Xi_Xi_By_Aux: " << I5_Xi_Xi_By_Aux << "\n";

    // check generic integrals for DD2

    double I5_Aux = bie::i5_Aux(q,d,delta);
    double I5_By_Aux = I5_Aux/pow(eta,2.0) + theta/(3.0*pow(eta,3.0)); // to verify when eta =! 0 against mma
    double I7_Aux = bie::i7_Aux(eta,q,d,D,delta);
    double I7_By_Aux = I7_Aux/pow(eta,4.0) + theta/(5.0*pow(eta,5.0)); // to verify when eta =! 0 against mma

    std::cout << " Checking generic integrals for DD3 ... " << "\n";

    //    std::cout << " I5_Aux: " << I5_Aux << "\n";
    std::cout << " I5_By_Aux: " << I5_By_Aux << "\n";
//    std::cout << " I7_Aux: " << I7_Aux << "\n";
    std::cout << " I7_By_Aux: " << I7_By_Aux << "\n";

    std::cout << "----------end of test 3DT0  ---------------------\n";

    return 0;

}

int main() {

  std::cout << "++++++++++++++++++++\n";
//  test3DR0();

  //test2DP1();

  //testS3DP0();

//  testFullMat();
//  testHdot();

//  test3DT6Mesh();
//    std::string vertices_file = "/home/alexis/bigwham/vertices5000.csv";
//    std::string connectivity_file = "/home/alexis/bigwham/connectivity5000.csv";
//  test3DT6PennyShaped(vertices_file,connectivity_file);

test3DT0();

  std::cout << "\n End of BigWham - exe " << "\n";

}

