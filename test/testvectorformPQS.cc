// testvectorformPQS.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// testing the vector form-factors partially quenched staggered
// as expressed in lowest order masses

#include<iostream>
#include<complex>
#include "vectorformPQ.h"
#include "vectorformPQS.h"

int main(void){
  using namespace std;


  double m11,m22,m33,m44,m55,m66,F0,mu,mp,mk,qsq,l9r;
  double DS,DV,DA,DT,dV,dA;// taste splittings and hairpin vertices
  mu = 0.77;
  mp = 0.135;
  mk = 0.495;
  F0 = 0.0923;
  l9r = 0.;
  m11 = mp*mp;
  m22 = 0.53*m11;
  m33 = 2.*mk*mk-m11;
  m44 = 1.1*m11;
  m55 = 1.11*m11;
  m66 = 1.2*m33;
  qsq = -0.25;
  DS = 0.01; DV = 0.015; DA = 0.03; DT = 0.08;
  dV = 0.045; dA = 0.06;


  cout << "Kl3 F_+ formfactors: partially quenched staggered\n";
  cout << "DS DV DA DT dV dA : "
       <<DS<<' '<< DV<<' '<< DA<<' '<< DT<<' '<< dV<<' '<< dA<<'\n';

  cout << "\nRpart only here\n";

  cout <<"2 valence 2 sea masses:\n";
  quarkmassnf qmass22({m11/2.,m33/2.,m44/2.,m66/2.},F0,mu);
  cout << qmass22;
  cout << "PQ staggered : "<<fvpv2s2nf3p4RS(qsq,qmass22,DS,DV,DA,DT,dV,dA) <<'\n';
  cout << "DS,dV,... = 0 reduces to PQ:\n";
  cout << "PQ v2s2      : "<<fvpv2s2nf3p4R(qsq,qmass22) <<'\n';
  cout << "PQ staggered : "<<fvpv2s2nf3p4RS(qsq,qmass22,0.,0.,0.,0.,0.,0.) <<'\n';

  cout <<"\n2 valence 3 sea masses:\n";
  quarkmassnf qmass23({m11/2.,m33/2.,m44/2.,m55/2.,m66/2.},F0,mu);
  cout << qmass23;
  cout << "PQ staggered : "<<fvpv2s3nf3p4RS(qsq,qmass23,DS,DV,DA,DT,dV,dA) <<'\n';
  cout << "DS,dV,... = 0 reduces to PQ:\n";
  cout << "PQ v2s3      : "<<fvpv2s3nf3p4R(qsq,qmass23) <<'\n';
  cout << "PQ staggered : "<<fvpv2s3nf3p4RS(qsq,qmass23,0.,0.,0.,0.,0.,0.) <<'\n';

  cout << "\n3 sea reduces to 2 sea\n";
  double eps = 1e-2, eps1 = 1.+eps,eps2=1.-eps;
  quarkmassnf qmass1({m11/2.,m33/2.,m44/2.,m44/2.*eps1,m66/2.},F0,0.77);
  quarkmassnf qmass2({m11/2.,m33/2.,m44/2.,m44/2.*eps2,m66/2.},F0,0.77);
  dcomplex fv1  = fvpv2s3nf3p4RS(qsq,qmass1,DS,DV,DA,DT,dV,dA);
  dcomplex fv2  = fvpv2s3nf3p4RS(qsq,qmass2,DS,DV,DA,DT,dV,dA);
  cout << "two points close to 2 sea mass case:\n"
       << fv1 <<' ' << fv2 <<'\n'
       << "average and comparison with 2 sea :\n"
       << (fv1+fv2)/2.<<' '<<fvpv2s2nf3p4RS(qsq,qmass22,DS,DV,DA,DT,dV,dA)<<'\n';

  cout << "\nNo dependence on the spectator mass\n";
  quarkmassnf qmass32({m11/2.,m22/2.,m33/2.,m44/2.,m66/2.},F0,mu);
  quarkmassnf qmass33({m11/2.,m22/2.,m33/2.,m44/2.,m55/2.,m66/2.},F0,mu);
  cout << "PQ v2s2 : "<<fvpv2s2nf3p4RS(qsq,qmass22,DS,DV,DA,DT,dV,dA) <<'\n';
  cout << "PQ v3s2 : "<<fvpv3s2nf3p4RS(qsq,qmass32,DS,DV,DA,DT,dV,dA) <<'\n';
  cout << "PQ v2s3 : "<<fvpv2s3nf3p4RS(qsq,qmass23,DS,DV,DA,DT,dV,dA) <<'\n';
  cout << "PQ v3s3 : "<<fvpv3s3nf3p4RS(qsq,qmass33,DS,DV,DA,DT,dV,dA) <<'\n';




  return 0;
}
