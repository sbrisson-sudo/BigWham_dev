//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 22.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2019.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <iostream>
#include <string>
#include <random>

#include <il/Array2D.h>

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
#include <src/solvers/HIterativeSolverUtilities.h>

#include <BigWham.h>

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

//  std::cout << ne << " elt " << conn.size(0) << " nodes " << nodes.size(0) <<"\n";
//  for (il::int_t i=0;i<nodes.size(0);i++){
//    std::cout << " nodes points " << nodes(i,0) <<" / " <<nodes(i,1) <<"\n";
//  }
//
//  for (il::int_t i=0;i<conn.size(0);i++){
//    std::cout << " conn   " << conn(i,0) <<" / " <<conn(i,1) <<"\n";
//  }

  bie::Mesh test2(nodes,conn,1);
  il::Array2D<double> coll_points = test2.getCollocationPoints();

  il::int_t leaf_size=32;//test->max_leaf_size;
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
//  std::cout << " in set h mat 2 " << cluster.permutation.size() << " e aca " << test2->epsilon_aca <<"\n";

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
//    std::cout <<" x_s " << i <<" "<< xx[2*i] <<"\n";
//    std::cout <<" x_n " << i <<" "<< xx[2*i+1] <<"\n";
  }

 il::Array< double> y=il::dot(h_,xx);
  std::cout << " Elas times solution:: \n";
  for (il::int_t i=0;i<coll_points.size(0);i++){
    std::cout <<" y " << i  <<" " << y[2*i] <<"\n";
    std::cout <<" y " << i  <<" " << y[2*i+1] <<"\n";
  }

  // now testing the Bigwhamio class...
  Bigwhamio *test = new Bigwhamio();
//
  test->numberofnodes =nodes.size(0);
  test->numberofelements=conn.size(0);
  test->nodesperelement =2;
  test->dimension=2;

  int len = test->numberofnodes *2;
  double * pointlist = new double[len];
  int iter=0;
  for ( int i = 0; i < nodes.size(0); i++) {
    for (int j=0;j<nodes.size(1);j++){
       pointlist[iter++] = nodes(i,j);
    }
  }

  len = test->numberofelements*2;
  int * connlist = new int[len];
   iter=0;
  for ( int i = 0; i < conn.size(0); i++) {
    for (int j=0;j<conn.size(1);j++){
       connlist[iter++] = conn(i,j);
    }
  }
//
  test->setMesh(pointlist,connlist,1);
//  test->setCollocationPoints();
//
  test->eta=1.0;
  test->max_leaf_size=22;
  test->epsilon_aca=1.e-3;

  test->setElasticProperties(1.,0.2);

  test->setHmat();
  test->setDiag();

  std::cout << " build ?" << test->isHbuilt
            <<  " C R :"<<  test->compressionRatio() <<"\n";

  auto *x_rhs = (new double[coll_points.size(0)*2]) ;

  for (il::int_t i=0;i<coll_points.size(0)*2;i++){
    x_rhs[i]=xx[i];
  }

  double *sol = test->hdotProduct(x_rhs);

  for (il::int_t i=0;i<coll_points.size(0);i++){
    std::cout <<" y " << 2*i  <<" : " << y[2*i] <<" ----- ";
    std::cout << " y_b " << 2*i << " : " << sol[2*i] <<"\n";
    std::cout <<" y " <<  2*i+1   << " : " << y[2*i+1] <<"----- ";
    std::cout << " y_b " << 2*i+1 << " : " << sol[2*i+1] <<"\n";
    std::cout << " -- \n";
  }
 // il::Array2D<double> coll2=test2.getCollocationPoints();
//  std::cout << " coll points " << coll_points(0,0) <<" / " <<coll_points(1,0) <<"\n";

  std::cout << test->numberofelements << "=" << test->numberofnodes << " - "<< test->numberofcollocationpoints <<"\n";

  test->kernel="2DP1";

  std::cout << test->kernel <<"\n";
//  std::cout <<  " comp " << std::strcmp(test->kernel,s2) <<"\n";

//
  double * collpointlist=test->getCollocationPoints();

  std::cout << " Interpolation order " << test2.interpolationOrder() <<"\n";
  std::cout << " coll points " << coll_points(0,0) <<" / " <<coll_points(2*ne-1,0) <<"\n";
  std::cout << " coll points " <<  collpointlist[0] << " /// "<<  " /// "<<  collpointlist[(2*ne-1)*2]<<"\n";

  std::cout << " build ? " << test->isHbuilt
            <<  " C R :"<<  test->compressionRatio() <<"\n";
//    std::cout << " nu  " << test->nu <<  " yme :"<<  test->yme <<"\n";
  std::cout << " coll points " <<  collpointlist[2] << " /// "<<  " /// "<<  collpointlist[(2*ne-1)*2]<<"\n";

  auto * ff = (new  double [coll_points.size(0)*2]);
  auto * x0 = (new  double [coll_points.size(0)*2]);
  for (int i=0;i<coll_points.size(0);i++ ){
    ff[2*i]=0.;
    ff[2*i+1]=1.;
    x0[2*i]=0.;x0[2*i+1]=0.;
  }
  // gmres via bigwham object.

  int nIts = 1000; double tolerance=1.e-9;
  int stat= test->h_IterativeSolve(ff,false,1000,il::io, nIts,tolerance, x0);

  std::cout << " Gmress success (0 ?) " << stat <<"\n";

  // compute relative error
  il::Array<double> err{coll_points.size(0)-2,0.};
  for (int i=1;i<coll_points.size(0)-1;i++ ){
    err[i-1]=il::abs(x0[2*i+1]-wsol_nodes[i])/wsol_nodes[i];
  }
  std::cout << " rel err L2 norm  " << il::norm(err,il::Norm::L2) <<"\n";
  // hlu
  il::HMatrix<double> h2=test->h_;
  il::luDecomposition(0.001,il::io,h2);
  std::cout << "Compression ratio of HLU-decomposition:  "
            << il::compressionRatio(h2) << std::endl;
  double compress_ratio_lu = il::compressionRatio(h2);
  il::Array<double> ff2{coll_points.size(0)*2,1.};
  il::solve(h2, il::io, ff2.Edit());
  std::cout << ff2[10] << "\n";
  // compute relative error
  il::Array<double> errlu{coll_points.size(0)-2,0.};
  for (int i=1;i<coll_points.size(0)-1;i++ ){
    errlu[i-1]=il::abs(ff2[2*i+1]-wsol_nodes[i])/wsol_nodes[i];
  }
  std::cout << " rel err L2 norm  " << il::norm(errlu,il::Norm::L2) <<"\n";


//  delete test2;
  std::cout << "-------------- End of 2DP1 - tests ------------- " << "\n";

  return  0;

}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
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
//    std::cout <<" x_s " << i <<" "<< xx[2*i] <<"\n";
//    std::cout <<" x_n " << i <<" "<< xx[2*i+1] <<"\n";
  }

  il::Array< double> y=il::dot(h_,xx);

  std::cout << " Now creating similar BigWham instance \n";
  std::string kernelname = "S3DP0";

  auto *test = new Bigwhamio();
//
  test->numberofnodes =nodes.size(0);
  test->numberofelements=conn.size(0);
  test->nodesperelement =2;

  int len = test->numberofnodes *2;
  double* pointlist = new double[len];
  int iter=0;
  for ( int i = 0; i < nodes.size(0); i++) {
    for (int j=0;j<nodes.size(1);j++){
       pointlist[iter++] = nodes(i,j);
    }
  }
  test->dimension=2;

  int lenc = test->numberofelements*2;
  int * connlist = new int[lenc];
  iter=0;
  for ( int i = 0; i < conn.size(0); i++) {
    for (int j=0;j<conn.size(1);j++){
      connlist[iter++] = conn(i,j);
    }
  }
 test->setMesh(pointlist,connlist,0);

//  test->setCollocationPoints();
//
  std::cout << " preparing bigwham instance for hmat \n";
//
  std::cout << "coll points set" <<"\n";

  test->eta=4.0;
  test->max_leaf_size=22;
  test->epsilon_aca=1.e-3;
  test->kernel = kernelname;

  test->setElasticProperties(1.,0.2);

  test->setHmat();
  test->setDiag();

  std::cout << " build ?" << test->isHbuilt
            <<  " CR :"<<  test->compressionRatio() << " CR   " <<"\n";

  auto *x_rhs = (new double[coll_points.size(0)*2]) ;

  for (il::int_t i=0;i<coll_points.size(0)*2;i++){
    x_rhs[i]=xx[i];
  }

  //////////// dot product.
  double *sol = test->hdotProduct(x_rhs);

//  for (il::int_t i=0;i<coll_points.size(0);i++){
//    std::cout <<" y " << 2*i  <<" : " << y[2*i] <<" ----- ";
//    std::cout << " y_b " << 2*i << " : " << sol[2*i] <<"\n";
//    std::cout <<" y " <<  2*i+1   << " : " << y[2*i+1] <<"----- ";
//    std::cout << " y_b " << 2*i+1 << " : " << sol[2*i+1] <<"\n";
//    std::cout << " -- \n";
//  }

  //
  std::cout << " TEst h - pattern \n";

//

  std::string hpatfilename = "/Users/bricelecampion/ClionProjects/BigWham/cmake-build-debug/Hpat.csv";

  il::output_hmatPatternF(h_,hpatfilename);

//  il::Array2D<il::int_t> pat = il::output_hmatPattern(test->h_);
//  std::cout << " out from hmat pattern\n";
//
//  for (il::int_t j=0;j<pat.size(1);j++){
//    for (il::int_t i=0;i<pat.size(0);i++){
//      std::cout << pat(i,j) << " - " ;
//    }
//    std::cout << std::endl;
//  }

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

  il::SparseMatrixCSR<il::int_t, double>  sparsy=test->getFullBlocks();

  double * val_list = (double *) nullptr ;
  int * pos_list = (int *) nullptr;
  int nentries=0;

  std::cout << "test with pointers " <<"\n";
  test->getFullBlocks(val_list, pos_list,nentries);
  std::cout << val_list[0] <<"\n";

  il::spot_t ss = (test->h_.root());

   std::cout << " root " << ss.index << " sign " << " valid ? " << ss.isValid() <<"\n";


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

int test3DT6Mesh() {
// testing the 3DT6 kernel implementation

    std::cout << " testing the 3DT6 kernel implementation" << "\n";
// simple 2D mesh

// transposed!

    il::Array2D<double> nodes{5, 3, 0.};
    il::Array2D<il::int_t> conn{4, 3};

    nodes(0, 0) = 0;
    nodes(0, 1) = 0;
    nodes(0, 2) = 0;

    nodes(1, 0) = 2;
    nodes(1, 1) = 0;
    nodes(1, 2) = 0;

    nodes(2, 0) = 0;
    nodes(2, 1) = 2;
    nodes(2, 2) = 0;

    nodes(3, 0) = -2;
    nodes(3, 1) = 0;
    nodes(3, 2) = 0;

    nodes(4, 0) = 0;
    nodes(4, 1) = -2;
    nodes(4, 2) = 0;

    conn(0, 0) = 0;
    conn(0, 1) = 1;
    conn(0, 2) = 2;

    conn(1, 0) = 0;
    conn(1, 1) = 2;
    conn(1, 2) = 3;

    conn(2, 0) = 0;
    conn(2, 1) = 3;
    conn(2, 2) = 4;

    conn(3, 0) = 0;
    conn(3, 1) = 4;
    conn(3, 2) = 1;

    bie::Mesh3D Mesh2(nodes, conn, 2);

    il::Array2D<double> coll_points = Mesh2.getCollocationPoints();
    il::int_t el = 1; // element number to be tested
    std::cout << "collocation points element " << el << "\n";
    std::cout << "size = " << coll_points.size(0) << " x " << coll_points.size(1) << "\n";

    std::cout << coll_points(0 + 6 * el,0) << " " << coll_points(0 + 6 * el,1) << " " << coll_points(0 + 6 * el,2) <<  "\n";
    std::cout << coll_points(1 + 6 * el,0) << " " << coll_points(1 + 6 * el,1) << " " << coll_points(1 + 6 * el,2) <<  "\n";
    std::cout << coll_points(2 + 6 * el,0) << " " << coll_points(2 + 6 * el,1) << " " << coll_points(2 + 6 * el,2) <<  "\n";
    std::cout << coll_points(3 + 6 * el,0) << " " << coll_points(3 + 6 * el,1) << " " << coll_points(3 + 6 * el,2) <<  "\n";
    std::cout << coll_points(4 + 6 * el,0) << " " << coll_points(4 + 6 * el,1) << " " << coll_points(4 + 6 * el,2) <<  "\n";
    std::cout << coll_points(5 + 6 * el,0) << " " << coll_points(5 + 6 * el,1) << " " << coll_points(5 + 6 * el,2) <<  "\n";

    il::StaticArray2D<double,3,3> vert_elt = Mesh2.getVerticesElt(el);

    std::cout << "Vertices element " << el << "\n";
    std::cout << vert_elt(0,0) << " " << vert_elt(0,1) << " " << vert_elt(0,2) <<  "\n";
    std::cout << vert_elt(1,0) << " " << vert_elt(1,1) << " " << vert_elt(1,2) <<  "\n";
    std::cout << vert_elt(2,0) << " " << vert_elt(2,1) << " " << vert_elt(2,2) <<  "\n";

    // testing normal and tangent unit vector computations with new constructor

    // unit normal to element in global system of coordinates
    il::StaticArray<double, 3> n_;
    // unit tangent 1 to element in global system of coordinates
    il::StaticArray<double, 2> s_;
    // unit tangent 2 to element in global system of coordinates
    il::StaticArray<double, 2> t_;
    // centroid in global system of coordinates
    il::StaticArray<double, 2> Xc_;

    // prepare vectors for cross product
    // vec01 goes from vertex 0 to vertex 1
    // vec02 goes from vertex 0 to vertex 2
    il::StaticArray2D<double,3,3> Xs = vert_elt;
    il::StaticArray<double,3> vec01, vec02;
    vec01[0] = Xs(1,0) - Xs(0,0);
    vec01[1] = Xs(1,1) - Xs(0,1);
    vec01[2] = Xs(1,2) - Xs(0,2);
    vec02[0] = Xs(2,0) - Xs(0,0);
    vec02[1] = Xs(2,1) - Xs(0,1);
    vec02[2] = Xs(2,2) - Xs(0,2);

    // cross product to get normal vector
    il::StaticArray<double, 3> n;
    n[0] = vec01[1] * vec02[2] - vec01[2] * vec02[1];
    n[1] = vec01[2] * vec02[0] - vec01[0] * vec02[2];
    n[2] = vec01[0] * vec02[1] - vec01[1] * vec02[0];

    // normalize normal vector
    double norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] = n[0]/norm;
    n[1] = n[1]/norm;
    n[2] = n[2]/norm;

    std::cout << "unit normal vector " << "\n";
    std::cout << n[0] << " " << n[1] << " " << n[2] <<  "\n";

    // compute unit tangent vector s_ that goes from vertex 0 to vertex 1
    il::StaticArray<double,3> s;
    s[0] = Xs(1,0) - Xs(0,0);
    s[1] = Xs(1,1) - Xs(0,1);
    s[2] = Xs(1,2) - Xs(0,2);

    // normalize tangent vector s_
    norm = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
    s[0] = s[0]/norm;
    s[1] = s[1]/norm;
    s[2] = s[2]/norm;

    std::cout << "unit tangent vector s_ " << "\n";
    std::cout << s[0] << " " << s[1] << " " << s[2] <<  "\n";

    // compute unit tangent vector t_ orthogonal to n_ and s_

    // cross product between n_ and s_
    il::StaticArray<double, 3> t;
    t[0] = n[1] * s[2] - n[2] * s[1];
    t[1] = n[2] * s[0] - n[0] * s[2];
    t[2] = n[0] * s[1] - n[1] * s[0];

    // normalize tangent vector t_
    norm = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
    t[0] = t[0]/norm;
    t[1] = t[1]/norm;
    t[2] = t[2]/norm;

    // compute centroid of the element
    il::StaticArray<double, 3> xc;
    xc[0] = ( Xs(0,0) + Xs(1,0) + Xs(2,0) ) / 3.;
    xc[1] = ( Xs(0,1) + Xs(1,1) + Xs(2,1) ) / 3.;
    xc[2] = ( Xs(0,2) + Xs(1,2) + Xs(2,2) ) / 3.;

    std::cout << "centroid " << "\n";
    std::cout << xc[0] << " " << xc[1] << " " << xc[2] <<  "\n";

    std::cout << "unit tangent vector t_ " << "\n";
    std::cout << t[0] << " " << t[1] << " " << t[2] <<  "\n";

    std::cout << "testing nodes and collocation points for T0, T3 and T6 " << "\n";

    std::cout << "nodes T0" << "\n";
    bie::TriangularElementData elt0(vert_elt, 0);
    il::Array2D<double> n_elt0 = elt0.getNodes();
    std::cout << "Size = " << n_elt0.size(0) << " x " << n_elt0.size(1) <<  "\n";
    std::cout << n_elt0(0,0) << " " << n_elt0(0,1) << " " << n_elt0(0,2) <<  "\n";

    std::cout << "CPs T0" << "\n";
    il::Array2D<double> cp_elt0 = elt0.getCollocationPoints();
    std::cout << "Size = " << cp_elt0.size(0) << " x " << cp_elt0.size(1) <<  "\n";
    std::cout << cp_elt0(0,0) << " " << cp_elt0(0,1) << " " << cp_elt0(0,2) <<  "\n";

    std::cout << "Rotation Matrix" << "\n";
    il::StaticArray2D<double, 3, 3> rM = elt0.rotationMatrix();
    std::cout << rM(0,0) << " " << rM(0,1) << " " << rM(0,2) <<  "\n";
    std::cout << rM(1,0) << " " << rM(1,1) << " " << rM(1,2) <<  "\n";
    std::cout << rM(2,0) << " " << rM(2,1) << " " << rM(2,2) <<  "\n";

    std::cout << "nodes T3" << "\n";
    bie::TriangularElementData elt1(vert_elt, 1);
    il::Array2D<double> n_elt1 = elt1.getNodes();
    std::cout << "Size = " << n_elt1.size(0) << " x " << n_elt1.size(1) <<  "\n";
    std::cout << n_elt1(0,0) << " " << n_elt1(0,1) << " " << n_elt1(0,2) <<  "\n";
    std::cout << n_elt1(1,0) << " " << n_elt1(1,1) << " " << n_elt1(1,2) <<  "\n";
    std::cout << n_elt1(2,0) << " " << n_elt1(2,1) << " " << n_elt1(2,2) <<  "\n";

    std::cout << "CPs T3" << "\n";
    il::Array2D<double> cp_elt1 = elt1.getCollocationPoints();
    std::cout << "Size = " << cp_elt1.size(0) << " x " << cp_elt1.size(1) <<  "\n";
    std::cout << cp_elt1(0,0) << " " << cp_elt1(0,1) << " " << cp_elt1(0,2) <<  "\n";
    std::cout << cp_elt1(1,0) << " " << cp_elt1(1,1) << " " << cp_elt1(1,2) <<  "\n";
    std::cout << cp_elt1(2,0) << " " << cp_elt1(2,1) << " " << cp_elt1(2,2) <<  "\n";

    std::cout << "nodes T6" << "\n";
    bie::TriangularElementData elt2(vert_elt, 2);
    il::Array2D<double> n_elt2 = elt2.getNodes();
    std::cout << "Size = " << n_elt2.size(0) << " x " << n_elt2.size(1) <<  "\n";
    std::cout << n_elt2(0,0) << " " << n_elt2(0,1) << " " << n_elt2(0,2) <<  "\n";
    std::cout << n_elt2(1,0) << " " << n_elt2(1,1) << " " << n_elt2(1,2) <<  "\n";
    std::cout << n_elt2(2,0) << " " << n_elt2(2,1) << " " << n_elt2(2,2) <<  "\n";
    std::cout << n_elt2(3,0) << " " << n_elt2(3,1) << " " << n_elt2(3,2) <<  "\n";
    std::cout << n_elt2(4,0) << " " << n_elt2(4,1) << " " << n_elt2(4,2) <<  "\n";
    std::cout << n_elt2(5,0) << " " << n_elt2(5,1) << " " << n_elt2(5,2) <<  "\n";

    std::cout << "CPs T6" << "\n";
    il::Array2D<double> cp_elt2 = elt2.getCollocationPoints();
    std::cout << "Size = " << cp_elt2.size(0) << " x " << cp_elt2.size(1) <<  "\n";
    std::cout << cp_elt2(0,0) << " " << cp_elt2(0,1) << " " << cp_elt2(0,2) <<  "\n";
    std::cout << cp_elt2(1,0) << " " << cp_elt2(1,1) << " " << cp_elt2(1,2) <<  "\n";
    std::cout << cp_elt2(2,0) << " " << cp_elt2(2,1) << " " << cp_elt2(2,2) <<  "\n";
    std::cout << cp_elt2(3,0) << " " << cp_elt2(3,1) << " " << cp_elt2(3,2) <<  "\n";
    std::cout << cp_elt2(4,0) << " " << cp_elt2(4,1) << " " << cp_elt2(4,2) <<  "\n";
    std::cout << cp_elt2(5,0) << " " << cp_elt2(5,1) << " " << cp_elt2(5,2) <<  "\n";

    std::cout << " Creating BigWham instance \n";
    std::string kernelname = "3DT6";

    auto *test = new Bigwhamio();

    test->numberofnodes =nodes.size(0); // these are actually vertices in 3D
    test->numberofelements=conn.size(0);
    test->nodesperelement =3;
    test->dimension =3;

    int len = test->numberofnodes * test->dimension;
    double* pointlist = new double[len];
    int iter=0;
    for ( int i = 0; i < test->numberofnodes; i++) {
        for (int j = 0; j < test->dimension; j++){
            pointlist[iter++] = nodes(i,j);
        }
    }

    int lenc = test->numberofelements * test->nodesperelement;
    int * connlist = new int[lenc];
    iter=0;
    for ( int i = 0; i < test->numberofelements; i++) {
        for (int j = 0; j < test->nodesperelement; j++){
            connlist[iter++] = conn(i,j);
        }
    }

    test->setMesh(pointlist,connlist,2);

//    std::cout << "Preparing bigwham instance for hmat \n";
//
//    test->eta=25;
//    test->max_leaf_size=25;
//    test->epsilon_aca=1.e-4;
//    test->kernel = kernelname;
//
//    test->setElasticProperties(1.,0.2);
//
//
//    test->setHmat();
//
//    test->setDiag();
//
//    std::cout << " build ?" << test->isHbuilt
//              <<  " CR :"<<  test->compressionRatio() << " CR   " <<"\n";

    return 0;
}

// I (alexis) just copy/paste these two functions written by Carlo, to read and load the connectivity and coordinate matrices from a CSV file

il::Array2D<double> read_coord_CSV(std::string& inputfileCOORD){

    il::Array2D<double> coordinates_table;

    std::string line = "", subline = "";
    char delim = ',';

    std::ifstream nf; // node coordinate file stream
    nf.open(inputfileCOORD.c_str());
    if(nf.is_open()) {
        // counting the array size
        il::int_t n_row = 0, n_col = 0;
        while(!nf.eof()) {
            std::getline(nf, line);
            if (line.length() == 0) {
                std::cout << "Empty row" << std::endl;
            } else {
                std::stringstream linestream(line);
                ++n_row;
                il::int_t n_col_t = 0;
                do {
                    subline = "";
                    std::getline(linestream, subline, delim);
                    if (subline.length() > 0){
                        ++n_col_t;
                    }
                } while(subline.length() > 0);
                if (n_col_t != n_col && n_col > 0) {
                    std::cout << "Row size is not constant" << std::endl;
                }
                else n_col = n_col_t;
            }
        }
        nf.clear();
        nf.seekg(0);
        // resizing the output array (node coordinates)
        coordinates_table.Resize(n_row,n_col);
        // importing the array
        for (il::int_t k = 0; k < coordinates_table.size(0); ++k) {
            std::getline(nf, line);
            std::stringstream linestream(line);
            for (il::int_t j = 0; j < coordinates_table.size(1); ++j) {
                subline = "";
                std::getline(linestream, subline, delim);
                if (subline.length() > 0){
                    double value = -1.0;
                    std::stringstream sublinestream(subline);
                    sublinestream >> value;
                    coordinates_table(k,j) = value;
                }
            }
        }
        nf.close();
        return coordinates_table;
    }
    else {
        std::cout << "Can't open the file" << std::endl;
    }
}

il::Array2D<il::int_t> read_conn_CSV(std::string& inputfileCONN){

    il::Array2D<il::int_t> connectivity_table;

    std::string line = "", subline = "";
    char delim = ',';
    std::ifstream cf; // connectivity file stream
    cf.open(inputfileCONN.c_str());
    if(cf.is_open()) {
        // counting the array size
        il::int_t n_row = 0, n_col = 0;
        while(!cf.eof()) {
            std::getline(cf, line);
            if (line.length() == 0) {
                std::cout << "Empty row" << std::endl;
            } else {
                std::stringstream linestream(line);
                ++n_row;
                il::int_t n_col_t = 0;
                do {
                    subline = "";
                    std::getline(linestream, subline, delim);
                    if (subline.length() > 0){
                        ++n_col_t;
                    }
                } while(subline.length() > 0);
                if (n_col_t != n_col && n_col > 0) {
                    std::cout << "Row size is not constant" << std::endl;
                }
                else n_col = n_col_t;
            }
        }
        cf.clear();
        cf.seekg(0);
        // resizing the output array (mesh connectivity)
        connectivity_table.Resize(n_row,n_col);
        // importing the array
        for (il::int_t k = 0; k < connectivity_table.size(0); ++k) {
            std::getline(cf, line);
            std::stringstream linestream(line);
            for (il::int_t j = 0; j < connectivity_table.size(1); ++j) {
                subline = "";
                std::getline(linestream, subline, delim);
                if (subline.length() > 0){
                    il::int_t value = -1;
                    std::stringstream sublinestream(subline);
                    sublinestream >> value;
                    connectivity_table(k, j) = value;
                }
            }
        }

        cf.close();
        return connectivity_table;
    }
    else {
        std::cout << "Can't open the file" << std::endl;
    }

}

int test3DT6PennyShaped(std::string& vertices_file, std::string& connectivity_file) {

// testing the 3DT6 kernel implementation

    std::cout << " testing the 3DT6 kernel implementation with penny-shaped crack" << "\n";

    il::Array2D<double> nodes = read_coord_CSV(vertices_file);
//    std::cout << nodes(0,0) << " " << nodes(0,1) << " " << nodes(0,2) <<  "\n";
//    std::cout << nodes(1,0) << " " << nodes(1,1) << " " << nodes(1,2) <<  "\n";
//    std::cout << nodes(2,0) << " " << nodes(2,1) << " " << nodes(2,2) <<  "\n";
//    std::cout << nodes(3,0) << " " << nodes(3,1) << " " << nodes(3,2) <<  "\n";

    il::Array2D<il::int_t> conn = read_conn_CSV(connectivity_file);
//    std::cout << conn(0,0) << " " << conn(0,1) << " " << conn(0,2) <<  "\n";
//    std::cout << conn(1,0) << " " << conn(1,1) << " " << conn(1,2) <<  "\n";
//    std::cout << conn(2,0) << " " << conn(2,1) << " " << conn(2,2) <<  "\n";
//    std::cout << conn(3,0) << " " << conn(3,1) << " " << conn(3,2) <<  "\n";

    bie::Mesh3D Mesh(nodes, conn, 2);

    std::cout << " Creating BigWham instance \n";
    auto *test = new Bigwhamio();
    std::string kernelname = "3DT6";
    test->kernel = kernelname;
    test->setElasticProperties(1.,0.2);

    test->numberofnodes =nodes.size(0); // these are actually vertices in 3D
    test->numberofelements=conn.size(0);
    test->nodesperelement =3;
    test->dimension =3;

    int len = test->numberofnodes * test->dimension;
    double* pointlist = new double[len];
    int iter=0;
    for ( int i = 0; i < test->numberofnodes; i++) {
        for (int j = 0; j < test->dimension; j++){
            pointlist[iter++] = nodes(i,j);
        }
    }

    int lenc = test->numberofelements * test->nodesperelement;
    int * connlist = new int[lenc];
    iter=0;
    for ( int i = 0; i < test->numberofelements; i++) {
        for (int j = 0; j < test->nodesperelement; j++){
            connlist[iter++] = conn(i,j);
        }
    }

    test->setMesh(pointlist,connlist,2);

    std::cout << "Preparing bigwham instance for hmat \n";

    test->eta=10;
    test->max_leaf_size=10;
    test->epsilon_aca=1.e-4;

    test->setHmat();

    std::cout << " build ? " << test->isHbuilt
              <<  " CR :"<<  test->compressionRatio() <<"\n";

    // analytical solution at the nodes

    // compute radius of every node

//    double *coordList = test -> getNodesCoordinates();
//
//    int numRealNodes = test -> numberofelements * 6;
//
//    il::Array2D<double> coordMatrix{numRealNodes, 3, 0};
//
//    int index = 0;
//    for (il::int_t i = 0; i < numRealNodes; i++) {
//        coordMatrix(i, 0) = coordList[index++];
//        coordMatrix(i, 1) = coordList[index++];
//        coordMatrix(i, 2) = coordList[index++];
//    }
//
//    double* radius = new double[numRealNodes];
//    for ( int i = 0; i < numRealNodes; i++) {
//        radius[i] = sqrt(pow(coordMatrix(i,0), 2) + pow(coordMatrix(i,1), 2));
//    }
//
//    std::cout << radius[200] << " " <<  radius[500] <<  " " << radius[1000] << "\n";
//    std::cout << coordMatrix(200,0) << " " <<  coordMatrix(200,1) << " " <<  coordMatrix(200,2) << "\n";
//
//    // compute RHS for mode I
//    auto* ddsol_nodes = new double[numRealNodes * 3];
//    for ( int i = 0; i < numRealNodes; i++) {
//        ddsol_nodes[3*i] = 0.0;
//        ddsol_nodes[3*i+1] = 0.0;
//        ddsol_nodes[3*i+2] = ( 8.0 * 1.0 / ( il::pi * (1.0 / (1-pow(0.2,2)) ) ) ) * sqrt(std::abs(pow(5.0, 2) - pow(radius[i], 2)));
//    }
//
//    for ( int i = 0; i < numRealNodes * 3; i++) {
//        std::cout << ddsol_nodes[3*i] << " " <<  ddsol_nodes[3*i+1] <<  " " << ddsol_nodes[3*i+2] << "\n";
//    }
//
//    il::Timer tt;
//    tt.Start();
//    std::cout << "Calling Hdot -------------- \n";
//
//    double *sol;
////    sol = test->hdotProduct(ddsol_nodes);
//    sol = test->hdotProductInPermutted(ddsol_nodes);
//
//    tt.Stop();
//    std::cout << "H-dot time = :  " << tt.time() << "\n";
//
//    for ( int i = 0; i < numRealNodes * 3; i++) {
//        std::cout << sol[3*i] << " " <<  sol[3*i+1] <<  " " << sol[3*i+2] << "\n";
//    }

//
//    for (int i = 0; i < coll_points.size(0); ++i) {
//        if (std::abs(coll_points(i,0)) < a) {
//            wsol_coll[i] = coef * sqrt(pow(a, 2) - pow(coll_points(i,0), 2));
//        }
//    }


//    double Ep=1.0/(1.0-0.2*0.2);
//    double sig = 1.0;
//    double a=1.0;
//
//    double coef = 4.0 * sig / (Ep);
//
//    il::Array<double> wsol_coll{coll_points.size(0), 0.};
//
//    for (int i = 0; i < coll_points.size(0); ++i) {
//        if (std::abs(coll_points(i,0)) < a) {
//            wsol_coll[i] = coef * sqrt(pow(a, 2) - pow(coll_points(i,0), 2));
//        }
//    }
    //at corresponding nodes - works here due to simple mesh....

//    il::Array<double> xx{coll_points.size(0)*2};
//
//    for (il::int_t i=0;i<xx.size()/2;i++){
//        xx[2*i]=0.0;
//        xx[2*i+1]=wsol_coll[i];
//
//
//
//    auto *x_rhs = (new double[coll_points.size(0)*2]) ;
//
//    for (il::int_t i=0;i<coll_points.size(0)*2;i++){
//        x_rhs[i]=xx[i];
//    }

//    double *sol = test->hdotProduct(x_rhs);

    return 0;

}

int main() {

  std::cout << "++++++++++++++++++++\n";

  test2DP1();
//
  testS3DP0();
//  testFullMat();
//  testHdot();

//  test3DT6Mesh();
//    std::string vertices_file = "/home/alexis/bigwham/vertices5000.csv";
//    std::string connectivity_file = "/home/alexis/bigwham/connectivity5000.csv";
//  test3DT6PennyShaped(vertices_file,connectivity_file);

  std::cout << " End of BigWham - exe " << "\n";

}

