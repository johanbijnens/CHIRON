// testmassdecayvevnf2V.cc is part of the CHIRON ChPT at two loops program
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
#include "massdecayvevnf2V.h"
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"
const double hbarc = 0.1973269718;// hbarc in units of GeV fm
int main(void){
  using namespace std;
  const double mpi = 0.13957061;
  const double fpi = 0.0922;
  const double mu = 0.77;
  const double mu2 = mu*mu;
  const double mp2 = mpi*mpi;
  const physmassnf2 stdmass(mpi,fpi,mu);
  libar lb1(-0.4,4.3,2.9,4.4);
  li li1(lb1);
  cout << stdmass;
  cout << lb1;
  cout << li1;
  // precision one and two loop integrals
  double tprecision = 1e-9;
  int besselsum = 200;
  setprecisionfinitevolumeoneloopt(tprecision,1e-9,true);
  setprecisionfinitevolumeoneloopb(besselsum,1e-9,true);
  setprecisionfinitevolumesunsetb(200,40,1e-13,1e-12,true);
  setprecisionfinitevolumesunsett(1e-13,1e-12,true);
  double L = 4./mpi;
  cout <<"hhVb "<<hhVb(mp2,mp2,mp2,mp2,L,mu2)<<'\n';
  cout <<"hhVt "<<hhVt(mp2,mp2,mp2,mp2,L,mu2)<<'\n';
  cout <<"hh21Vb "<<hh21Vb(mp2,mp2,mp2,mp2,L,mu2)<<'\n';
  cout <<"hh21Vt "<<hh21Vt(mp2,mp2,mp2,mp2,L,mu2)<<'\n';
  cout <<"hh27Vb "<<hh27Vb(mp2,mp2,mp2,mp2,L,mu2)<<'\n';
  cout <<"hh27Vt "<<hh27Vt(mp2,mp2,mp2,mp2,L,mu2)<<'\n';
  cout <<"mpi^2 phys          : "<<mpi*mpi<<'\n';
  cout <<"mpi^2*mpi4nf2 via li1 : "<<mpi*mpi*mpi4nf2(stdmass,li1)<<'\n';
  cout <<"mpi4nf2V Bessel     : "<<mpi4nf2Vb(stdmass,L)<<'\n';
  cout <<"mpi4nf2V theta      : "<<mpi4nf2Vt(stdmass,L)<<'\n';
  cout <<"mpi6lnf2V Bessel    : "<<mpi6lnf2Vb(stdmass,li1,L)<<'\n';
  cout <<"mpi6lnf2V theta     : "<<mpi6lnf2Vt(stdmass,li1,L)<<'\n';
  cout <<"mpi6nf2V Bessel     : "<<mpi6nf2Vb(stdmass,li1,L)<<'\n';
  cout <<"mpi6nf2V theta      : "<<mpi6nf2Vt(stdmass,li1,L)<<'\n';
  cout <<"fpi phys            : "<<stdmass.getfpi()<<'\n';
  cout <<"fpi*fpinf24 via li1 : "<<fpi4nf2(stdmass,li1)*stdmass.getfpi()<<'\n';
  cout <<"fpi4nf2V Bessel     : "<<fpi4nf2Vb(stdmass,L)<<'\n';
  cout <<"fpi4nf2V theta      : "<<fpi4nf2Vt(stdmass,L)<<'\n';
  cout <<"fpi6lnf2V Bessel    : "<<fpi6lnf2Vb(stdmass,li1,L)<<'\n';
  cout <<"fpi6lnf2V theta     : "<<fpi6lnf2Vt(stdmass,li1,L)<<'\n';
  cout <<"fpi6nf2V Bessel     : "<<fpi6nf2Vb(stdmass,li1,L)<<'\n';
  cout <<"fpi6nf2V theta      : "<<fpi6nf2Vt(stdmass,li1,L)<<'\n';
  cout <<"qqup4nf2 via li1    : "<<qqup4nf2(stdmass,li1)<<'\n';
  cout <<"qqup4nf2V Bessel    : "<<qqup4nf2Vb(stdmass,L)<<'\n';
  cout <<"qqup4nf2V theta     : "<<qqup4nf2Vt(stdmass,L)<<'\n';
 return 0;
}
