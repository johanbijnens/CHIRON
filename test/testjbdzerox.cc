// testjbdzerox.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include<iostream>
#include "jbnumlib.h"

// small testing function
double ff(const double x){
  return cos(x);
}

int main(void){
  using namespace std;
  double result = jbdzerox(ff,1.,4.,1e-10,1000,0);
  cout << result <<' '<<result*2./M_PI-1. << endl;
  result = jbdzerox(ff,1.,4.,1e-10,1000,1);
  cout << result <<' '<<result*2./M_PI-1. << endl;
  return 0;
}
