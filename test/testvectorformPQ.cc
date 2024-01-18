// testvectorformPQ.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// testing the vector form-factors as expressed in lowest order masses

#include<iostream>
#include<complex>
#include "vectorformPQ.h"
#include "vectorformlo.h"

int main(void){
  using namespace std;
  double m1 = 0.13956995;
  double mk = 0.493677;
  double F0 = 0.0924;
  double mu = 0.77;

  lomass lomass1(m1,mk,F0,mu);
  double m11 = pow(m1,2);
  double m33 = 2.*mk*mk-m11;

  Li Li1;
  Li1.setmu(mu);
  Li1.setli(9,0.0069);
  Linf Lin;
  Lin.setmu(mu);
  Lin.setlinf(9,0.0069);

  quarkmassnf qmassn({m11/2.,m33/2.},F0,0.77);

  double qsq = -0.1;
  cout << "Kl3 F_+ formfactors\n\n";
  cout << "L_9^r dependent parts\n";
  cout << "standard : "<<fvpkpip4Llo(qsq,lomass1,Li1)<<'\n';
  cout << "PQ       : "<<fvpnf3p4L(qsq,qmassn,Lin)<<'\n';

  cout << "\nRpart 2 valence and 2 sea quark masses\n";
  quarkmassnf qmass22({m11/2.,m33/2.,m11/1.5,m33/1.08},F0,0.77);
  cout << "example point f_+ : "<<fvpv2s2nf3p4R(qsq,qmass22) <<'\n';
  cout << "example point f_+ : "<<fvpv2s2nf3p4Rp(qsq,qmass22) 
       <<" other program\n";


  double eps = 1e-2, eps1 = 1.+eps,eps2=1.-eps;
  quarkmassnf qmass1({m11/2.,m33/2.,m11/2.*eps1,m33/2.*eps1},F0,0.77);
  quarkmassnf qmass2({m11/2.,m33/2.,m11/2.*eps2,m33/2.*eps1},F0,0.77);
  quarkmassnf qmass3({m11/2.,m33/2.,m11/2.*eps1,m33/2.*eps2},F0,0.77);
  quarkmassnf qmass4({m11/2.,m33/2.,m11/2.*eps2,m33/2.*eps2},F0,0.77);
  dcomplex fv1  = fvpv2s2nf3p4R(qsq,qmass1);
  dcomplex fv2  = fvpv2s2nf3p4R(qsq,qmass2);
  dcomplex fv3  = fvpv2s2nf3p4R(qsq,qmass3);
  dcomplex fv4  = fvpv2s2nf3p4R(qsq,qmass4);
  cout << "four points close to unquenched f_+ :\n"
       << fv1 <<' ' << fv2 <<' ' << fv3 <<' ' << fv4<<'\n'
       << "average and comparison with unquenched :\n"
       << (fv1+fv2+fv3+fv4)/4. <<' '<<fvpkpip4Rlo(qsq,lomass1)<<'\n';

  cout << "\nRpart 2 valence and 3 sea quark masses\n";
  quarkmassnf qmass23({m11/2.,m33/2.,m11/1.5,m11/1.2,m33/1.08},F0,0.77);
  cout << "example point f_+ : "<<fvpv2s3nf3p4R(qsq,qmass23) <<'\n';
  cout << "example point f_+ : "<<fvpv2s3nf3p4Rp(qsq,qmass23) 
       <<" other program\n";

  quarkmassnf qmass32({m11/2.,0.11,m33/2.,m11/1.5,m33/1.08},F0,0.77);
  quarkmassnf qmass33({m11/2.,0.11,m33/2.,m11/1.5,m11/1.2,m33/1.08},F0,0.77);
  cout << "\nThe spectator mass does not matter\n";
  cout << "example point f_+ 22: "<<fvpv2s2nf3p4R(qsq,qmass22) <<'\n';
  cout << "example point f_+ 32: "<<fvpv3s2nf3p4R(qsq,qmass32) <<'\n';
  cout << "example point f_+ 23: "<<fvpv2s3nf3p4R(qsq,qmass23) <<'\n';
  cout << "example point f_+ 32: "<<fvpv3s3nf3p4R(qsq,qmass33) <<'\n';

  cout << "\nThree sea mass case reduces to two sea mass case\n";
  quarkmassnf qmass231({m11/2.,m33/2.,m11/1.5,m11/1.5*eps1,m33/1.08},F0,0.77);
  quarkmassnf qmass232({m11/2.,m33/2.,m11/1.5,m11/1.5*eps2,m33/1.08},F0,0.77);
  cout << "Two points around 2 sea mass case\n";
  fv1 = fvpv2s3nf3p4R(qsq,qmass231);
  fv2 = fvpv2s3nf3p4R(qsq,qmass232);
  cout << fv1<<' '<<fv2<<'\n';
  cout <<"average and two sea mass case\n"
       <<(fv1+fv2)/2.<<' '<<fvpv2s2nf3p4R(qsq,qmass22)<<'\n';


  return 0;
}
