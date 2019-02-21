// vectorformlo.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.



#include<complex>
#include "oneloopintegrals.h"
#include "vectorformlo.h"

typedef std::complex<double> dcomplex;

// only used internally here
dcomplex hhem(const double msq, const double t,const double mu2){
  return 0.5*Ab(msq,mu2)-B22b(msq,t,mu2);
}
dcomplex hhkpi(const double m1sq, const double m2sq, const double t,
	       const double mu2){
  return 0.25*(Ab(m1sq,mu2)+Ab(m2sq,mu2))-B22b(m1sq,m2sq,t,mu2);
}

//***************************************************************************
// electromagnetic form-factors
dcomplex fvpipp4lo(const double t, const lomass mass,const Li Liin){
  return fvpipp4Llo(t,mass,Liin)+fvpipp4Rlo(t,mass);
}

dcomplex fvpipp4Rlo(const double t, const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double mu2 = pow(mass.getmu(),2);
  return (2.*hhem(mp2,t,mu2)+hhem(mk2,t,mu2))/pow(mass.getf0(),2);
}

double fvpipp4Llo(const double t, const lomass mass, const Li Liin){
  return 2.*t*Liin.out(9)/pow(mass.getf0(),2);
}

dcomplex fvkpp4lo(const double t, const lomass mass,const Li Liin){
  return fvkpp4Llo(t,mass,Liin)+fvkpp4Rlo(t,mass);
}

dcomplex fvkpp4Rlo(const double t, const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double mu2 = pow(mass.getmu(),2);
  return (hhem(mp2,t,mu2)+2.*hhem(mk2,t,mu2))/pow(mass.getf0(),2);
}

double fvkpp4Llo(const double t, const lomass mass, const Li Liin){
  return 2.*t*Liin.out(9)/pow(mass.getf0(),2);
}

dcomplex fvkop4lo(const double t, const lomass mass,const Li Liin){
  return fvkop4Rlo(t,mass);
}

dcomplex fvkop4Rlo(const double t, const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double mu2 = pow(mass.getmu(),2);
  return (-hhem(mp2,t,mu2)+hhem(mk2,t,mu2))/pow(mass.getf0(),2);
}

double fvkop4Llo(const double t, const lomass mass, const Li Liin){
  return 0.;
}

//****************************************************************************
// Kpi or Kl3 formfactors
// f+
dcomplex fvpkpip4lo(const double t, const lomass mass,const Li Liin){
  return fvpkpip4Llo(t,mass,Liin)+fvpkpip4Rlo(t,mass);
}

double fvpkpip4Llo(const double t, const lomass mass, const Li Liin){
  return 2.*t*Liin.out(9)/pow(mass.getf0(),2);
}

dcomplex fvpkpip4Rlo(const double t, const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  return 1.5*(hhkpi(mp2,mk2,t,mu2)+hhkpi(mk2,me2,t,mu2))/pow(mass.getf0(),2);
}

// f_-
dcomplex fvmkpip4lo(const double t, const lomass mass,const Li Liin){
  return fvmkpip4Llo(t,mass,Liin)+fvmkpip4Rlo(t,mass);
}

double fvmkpip4Llo(const double t, const lomass mass, const Li Liin){
  double L5r = Liin.out(5);
  double L9r = Liin.out(9);
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  return (mk2-mp2)*(4.*L5r-2.*L9r)/pow(mass.getf0(),2);
}

// in the simplified form from when working on twisted boundary conditions
dcomplex fvmkpip4Rlo(const double t, const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  dcomplex result = -0.5*mp2*Bb(mp2,mk2,t,mu2)-0.5*mk2*Bb(mk2,me2,t,mu2)
    +(2.*mk2-mp2)*(B1b(mp2,mk2,t,mu2)+B1b(mk2,me2,t,mu2))
    -3./2.*(mk2-mp2)*(B21b(mp2,mk2,t,mu2)+B21b(mk2,me2,t,mu2));
  return result/pow(mass.getf0(),2);
}

// f_- R in the form of Bijnens Talavera 2003
// difference with above is using Passarino-Veltman relation
// q^2 B21 + B22 = ...
dcomplex fvmkpip4RloBT(const double qsq, const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double xmu2 = pow(mass.getmu(),2);
  double abmp2 = Ab(mp2,xmu2);
  double abmk2 = Ab(mk2,xmu2);
  double abme2 = Ab(me2,xmu2);
  dcomplex FM4462 =
       + Bb(mp2,mk2,qsq,xmu2) * ( 5./12.*qsq - 5./12.*mk2 - 1./12.*
         mp2 );

      FM4462 +=  + Bb(mk2,me2,qsq,xmu2) * ( 1./4.*qsq - 7./12.*mk2
          + 1./12.*mp2 );

      FM4462 +=  + B1b(mp2,mk2,qsq,xmu2) * (  - 5./12.*qsq + 19./12.*
         mk2 - 7./12.*mp2 );

      FM4462 +=  + B1b(mk2,me2,qsq,xmu2) * (  - 1./4.*qsq + 23./12.*
         mk2 - 11./12.*mp2 );

      FM4462 +=  + B21b(mp2,mk2,qsq,xmu2) * (  - 5./6.*qsq - 3./2.*
         mk2 + 3./2.*mp2 );

      FM4462 +=  + B21b(mk2,me2,qsq,xmu2) * (  - 1./2.*qsq - 3./2.*
         mk2 + 3./2.*mp2 );

      FM4462 +=  + B22b(mp2,mk2,qsq,xmu2) * (  - 5./6. );

      FM4462 +=  + B22b(mk2,me2,qsq,xmu2) * (  - 1./2. );

      FM4462 +=  + 7./12.*abmk2 - 5./12.*abmp2 + 1./2.*abme2;

  return FM4462/pow(mass.getf0(),2);
}
