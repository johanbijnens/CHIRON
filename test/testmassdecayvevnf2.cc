// testmassdecayvevnf2.cc is part of the CHIRON ChPT at two loops program
// collection
// Copyright (C) 2019 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include <iostream>
#include <iomanip>
#include <string>
#include <complex>
#include "massdecayvevnf2.h"

int main(void){
  using namespace std;
  const double mpi = 0.13957061;
  const physmassnf2 stdmass(mpi,0.0922,0.77);
  libar lb1(-0.4,4.3,2.9,4.4);
  li li1(lb1);
  cout << stdmass;
  cout << lb1;
  cout << li1;
  cout <<"mpi4 via lb1 : "<<mpi4nf2(stdmass,lb1)<<'\n';
  cout <<"mpi4 via li1 : "<<mpi4nf2(stdmass,li1)<<'\n';
  cout <<"mpi4r        : "<<mpi4rnf2(stdmass)<<'\n';
  cout <<"mpi4l        : "<<mpi4lnf2(stdmass,li1)<<'\n';
  cout <<"fpi4 via lb1 : "<<fpi4nf2(stdmass,lb1)<<'\n';
  cout <<"fpi4 via li1 : "<<fpi4nf2(stdmass,li1)<<'\n';
  cout <<"fpi4r        : "<<fpi4rnf2(stdmass)<<'\n';
  cout <<"fpi4l        : "<<fpi4lnf2(stdmass,li1)<<'\n';
  cout <<"qqup4 via lb1: "<<qqup4nf2(stdmass,lb1)<<'\n';
  cout <<"qqup4 via li1: "<<qqup4nf2(stdmass,li1)<<'\n';
  cout <<"qqup4r       : "<<qqup4rnf2(stdmass)<<'\n';
  cout <<"qqup4l       : "<<qqup4lnf2(stdmass,li1)<<'\n';
  return 0;
}
