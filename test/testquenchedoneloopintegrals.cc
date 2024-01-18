// testquenchedoneloopintegrals.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include<cmath>
#include<iostream>
#include "quenchedoneloopintegrals.h"
// note that that also includes oneloopintegrals.h


int main(void){
  using namespace std;
  double m1sq = 0.11;
  double m2sq = 0.7;
  double qsq = -0.13;
  double mu2 = pow(0.77,2);
  
  cout<< "nprop = 1 done correctly\n";
  cout<<"Bb  : "<<Bb(1,m1sq,m2sq,qsq,mu2)<<' '<<Bb(m1sq,m2sq,qsq,mu2)<<'\n';
  cout<<"B1b : "<<B1b(1,m1sq,m2sq,qsq,mu2)<<' '<<B1b(m1sq,m2sq,qsq,mu2)<<'\n';
  cout<<"B21b: "<<B21b(1,m1sq,m2sq,qsq,mu2)<<' '<<B21b(m1sq,m2sq,qsq,mu2)<<'\n';
  cout<<"B22b: "<<B22b(1,m1sq,m2sq,qsq,mu2)<<' '<<B22b(m1sq,m2sq,qsq,mu2)<<'\n';
  cout<<"B31b: "<<B31b(1,m1sq,m2sq,qsq,mu2)<<' '<<B31b(m1sq,m2sq,qsq,mu2)<<'\n';
  cout<<"B22b: "<<B32b(1,m1sq,m2sq,qsq,mu2)<<' '<<B32b(m1sq,m2sq,qsq,mu2)<<'\n';


  dcomplex deriv;
  double eps = 1e-3;
  cout << "\nnprop = 2 versus numerical derivative\n";
  deriv = (Bb(m1sq*(1.+eps),m2sq,qsq,mu2)-Bb(m1sq*(1.-eps),m2sq,qsq,mu2))
    /(2.*eps*m1sq);
  cout <<"Bb  : "<<Bb(2,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B1b(m1sq*(1.+eps),m2sq,qsq,mu2)-B1b(m1sq*(1.-eps),m2sq,qsq,mu2))
    /(2.*eps*m1sq);
  cout <<"B1b : "<<B1b(2,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B21b(m1sq*(1.+eps),m2sq,qsq,mu2)-B21b(m1sq*(1.-eps),m2sq,qsq,mu2))
    /(2.*eps*m1sq);
  cout <<"B21b: "<<B21b(2,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B22b(m1sq*(1.+eps),m2sq,qsq,mu2)-B22b(m1sq*(1.-eps),m2sq,qsq,mu2))
    /(2.*eps*m1sq);
  cout <<"B22b: "<<B22b(2,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B31b(m1sq*(1.+eps),m2sq,qsq,mu2)-B31b(m1sq*(1.-eps),m2sq,qsq,mu2))
    /(2.*eps*m1sq);
  cout <<"B31b: "<<B31b(2,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B32b(m1sq*(1.+eps),m2sq,qsq,mu2)-B32b(m1sq*(1.-eps),m2sq,qsq,mu2))
    /(2.*eps*m1sq);
  cout <<"B32b: "<<B32b(2,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';

  cout << "\nnprop = 3 versus numerical derivative\n";
  deriv = (Bb(m1sq,m2sq*(1.+eps),qsq,mu2)-Bb(m1sq,m2sq*(1.-eps),qsq,mu2))
    /(2.*eps*m2sq);
  cout <<"Bb  : "<<Bb(3,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B1b(m1sq,m2sq*(1.+eps),qsq,mu2)-B1b(m1sq,m2sq*(1.-eps),qsq,mu2))
    /(2.*eps*m2sq);
  cout <<"B1b : "<<B1b(3,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B21b(m1sq,m2sq*(1.+eps),qsq,mu2)-B21b(m1sq,m2sq*(1.-eps),qsq,mu2))
    /(2.*eps*m2sq);
  cout <<"B21b: "<<B21b(3,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B22b(m1sq,m2sq*(1.+eps),qsq,mu2)-B22b(m1sq,m2sq*(1.-eps),qsq,mu2))
    /(2.*eps*m2sq);
  cout <<"B22b: "<<B22b(3,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B31b(m1sq,m2sq*(1.+eps),qsq,mu2)-B31b(m1sq,m2sq*(1.-eps),qsq,mu2))
    /(2.*eps*m2sq);
  cout <<"B31b: "<<B31b(3,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';
  deriv = (B32b(m1sq,m2sq*(1.+eps),qsq,mu2)-B32b(m1sq,m2sq*(1.-eps),qsq,mu2))
    /(2.*eps*m2sq);
  cout <<"B32b: "<<B32b(3,m1sq,m2sq,qsq,mu2)<<' '<<deriv<<'\n';


  return 0;
}
