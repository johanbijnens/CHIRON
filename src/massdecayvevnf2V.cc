// massdecayvevnf2V.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2019 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains the isospin limit expressions for the finite volume
// corrections for masses and decay constants
// in terms of the infinite volume physical masses and Fpi 

// by setting the compile flag, it either produces the
// Besselfunction (CHIRONBESSEL) or the Theta function (CHIRONTHETA) output
// typically: large mL Bessel better
//            intermediate or small mL theta function version better
//   however jbdadmul as used in theta version has quite often problems
#include <iostream>
#include <cmath>

#include "inputsnf2.h"
#include "linf2.h"
#include "oneloopintegrals.h"
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"
#include "massdecayvevnf2V.h"

const double pi16 = 1./(16.*M_PI*M_PI);
const double pi162 = pi16*pi16;

// the code that produces results with Bessel method or theta method
#ifdef CHIRONTHETA
#define AbV AbVt
#define BbV BbVt
#define A23bV A23bVt
#define hhV hhVt
#define hh1V hh1Vt
#define hh21V hh21Vt
#define hh27V hh27Vt
#define hhdV hhdVt
#define hh1dV hh1dVt
#define hh21dV hh21dVt
#define hh27dV hh27dVt
#define mpi4nf2V mpi4nf2Vt
#define mpi6nf2V mpi6nf2Vt
#define mpi6lnf2V mpi6lnf2Vt
#define mpi6rnf2V mpi6rnf2Vt
#define fpi4nf2V fpi4nf2Vt
#define fpi6nf2V fpi6nf2Vt
#define fpi6lnf2V fpi6lnf2Vt
#define fpi6rnf2V fpi6rnf2Vt
#define qqup4nf2V qqup4nf2Vt
#else
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef CHIRONBESSEL
#define AbV AbVb
#define BbV BbVb
#define A23bV A23bVb
#define hhV hhVb
#define hh1V hh1Vb
#define hh21V hh21Vb
#define hh27V hh27Vb
#define hhdV hhdVb
#define hh1dV hh1dVb
#define hh21dV hh21dVb
#define hh27dV hh27dVb
#define mpi4nf2V mpi4nf2Vb
#define mpi6nf2V mpi6nf2Vb
#define mpi6lnf2V mpi6lnf2Vb
#define mpi6rnf2V mpi6rnf2Vb
#define fpi4nf2V fpi4nf2Vb
#define fpi6nf2V fpi6nf2Vb
#define fpi6lnf2V fpi6lnf2Vb
#define fpi6rnf2V fpi6rnf2Vb
#define qqup4nf2V qqup4nf2Vb
#else
// just some garbage to produce an error
x = 1./0.0.;
#endif
#endif
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double mpi4nf2V(const physmassnf2 mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double M4V = (  - 1./2.*AbV(mp2,xl)*mp2);
  return M4V/pow(mass.getfpi(),2);
}

double mpi6nf2V(const physmassnf2 mass, const li liin, const double xl){
  return  mpi6lnf2V(mass,liin,xl)+mpi6rnf2V(mass,xl);
}

double mpi6lnf2V(const physmassnf2 mass, const li liin, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double l1r,l2r,l3r,l4r,l5r,l6r,l7r,mu;
  liin.out(l1r,l2r,l3r,l4r,l5r,l6r,l7r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6lnf2V\n";}

  double massp6LV =
       + AbV(mp2,xl) * (  - l4r*pow(mp2,2) + 5*l3r*pow(mp2,2) + 8*l2r*
         pow(mp2,2) + 14*l1r*pow(mp2,2) );

      massp6LV +=  + A23bV(mp2,xl) * (  - 12*l2r*mp2 - 6*l1r*mp2 );

      return massp6LV/pow(mass.getfpi(),4);
}

double mpi6rnf2V(const physmassnf2 mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mu2 = pow(mass.getmu(),2);

  double massp6RV =
       + pow(AbV(mp2,xl),2) * (  - 3./8.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*BbV(mp2,xl) * ( 1./4.*pow(mp2,2) );

      massp6RV +=  + AbV(mp2,xl) * ( 13./12.*pi16*pow(mp2,2) - 7./4.*
         Ab(mp2,mu2)*mp2 );

      massp6RV +=  + hhV(mp2,mp2,mp2,mp2,xl,mu2) * ( 5./6.*pow(mp2,2)
          );

      massp6RV +=  + hh21V(mp2,mp2,mp2,mp2,xl,mu2) * ( 3*pow(mp2,2) );

      massp6RV +=  + hh27V(mp2,mp2,mp2,mp2,xl,mu2) * (  - 3*mp2 );

      return massp6RV/pow(mass.getfpi(),4);
}



//**********************************************************************
double fpi4nf2V(const physmassnf2 mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double F4V = (  AbV(mp2,xl));
  return F4V/pow(mass.getfpi(),1);
}

double fpi6nf2V(const physmassnf2 mass, const li liin, const double xl){
  return fpi6lnf2V(mass,liin,xl)+fpi6rnf2V(mass,xl);
}

double fpi6lnf2V(const physmassnf2 mass, const li liin, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double l1r,l2r,l3r,l4r,l5r,l6r,l7r,mu;
  liin.out(l1r,l2r,l3r,l4r,l5r,l6r,l7r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi6lnf2V\n";}

  double    decayp6LV =
       + AbV(mp2,xl) * ( 5./2.*l4r*mp2 - 4*l2r*mp2 - 7*l1r*mp2 );

      decayp6LV +=  + A23bV(mp2,xl) * ( 6*l2r + 3*l1r );

      return decayp6LV/pow(mass.getfpi(),3);
}

double fpi6rnf2V(const physmassnf2 mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mu2 = pow(mass.getmu(),2);

  double decayp6RV =
       + AbV(mp2,xl)*BbV(mp2,xl) * (  - 1./2.*mp2 );

      decayp6RV +=  + AbV(mp2,xl) * (  - 1./3.*pi16*mp2 + 3./2.*Ab(mp2,mu2));

      decayp6RV +=  + hhV(mp2,mp2,mp2,mp2,xl,mu2) * (  - 1./2.*mp2 );

      decayp6RV +=  + hh27V(mp2,mp2,mp2,mp2,xl,mu2) * ( 3./2. );

      decayp6RV +=  + hhV(mp2,mp2,mp2,mp2,xl,mu2) * ( 5./12.*pow(
         mp2,2) );

      decayp6RV +=  + hh21V(mp2,mp2,mp2,mp2,xl,mu2) * ( 3./2.*pow(
         mp2,2) );

      decayp6RV +=  + hh27dV(mp2,mp2,mp2,mp2,xl,mu2) * (  - 3./2.*mp2 )
         ;

      return decayp6RV/pow(mass.getfpi(),3);
}

//**********************************************************************
double qqup4nf2V(const physmassnf2 mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double F4V = (1.5*  AbV(mp2,xl));
  return F4V/pow(mass.getfpi(),2);
}


