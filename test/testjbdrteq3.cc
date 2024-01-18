// testjbdrteq3.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include<iostream>
#include "jbnumlib.h"

// roots of the equation x^3+r x^2+ s x+t=0
// d <=0 real roots
// d>0 2 complex roots, one real
using namespace std;
int main(void){

  double r=1.,s=1.,t=1.,x[3]={0.,0.,0.},d;

  r= -2.; s = -5.; t = 6.;// roots 1. -2. 3.
  jbdrteq3(r,s,t,x,d);
  cout << "roots: 1,-2,3 and the discriminant(<0)\n"; 
  cout << x[0] << ' ' << x[1] << ' '<<x[2] << ' '<<d << endl;

  r= 1.; s = -7.; t = -15.;
  jbdrteq3(r,s,t,x,d);
  cout << "roots: 1,-2+i,-2-i and the discriminant(>0)\n"; 
  cout << x[0] << ' ' << x[1] << ' '<<x[2] << ' '<<d << endl;

  r= 9.; s = 27.; t = 27.;
  cout << "Triple root: -3 and the discriminant (=0)\n";
  jbdrteq3(r,s,t,x,d);
  cout << x[0] << ' ' << x[1] << ' '<<x[2] << ' '<<d << endl;

  r= 9.; s = 27.; t = 27.;
  cout << "Triple root: -3 and the discriminant (=0)\n";
  jbdrteq3(r,s,t,x,d);
  cout << x[0] << ' ' << x[1] << ' '<<x[2] << ' '<<d << endl;

  r= 3.; s = 0.; t = -4.;
  cout << "Double+single root: 1,-2,-2 and the discriminant (=0)\n";
  jbdrteq3(r,s,t,x,d);
  cout << x[0] << ' ' << x[1] << ' '<<x[2] << ' '<<d << endl;
}
