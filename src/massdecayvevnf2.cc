// massdecayvevnf2.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains the isospin limit expressions for masses, decay
// constants and vacuum expectation values 
// in terms of physical masses and Fpi for the two-flavour case

#include <iostream>
#include <iomanip>
#include <string>

#include "inputsnf2.h"
#include "linf2.h"
//#include "ci.h"
#include "oneloopintegrals.h"
//#include "sunsetintegrals.h"
#include "massdecayvevnf2.h"

const double pi = M_PI;
const double pi16 = 1./(16.*pi*pi);

// mpi ******************************************************************

double mpi4nf2(const physmassnf2 mass, const libar libarin){
  return -0.5*pi16*pow(mass.getmpi()/mass.getfpi(),2)*libarin.out(3);
}

double mpi4nf2(const physmassnf2 mass, const li liin){
  return mpi4lnf2(mass,liin)+mpi4rnf2(mass);
}

double mpi4lnf2(const physmassnf2 mass, const li liin){
  double l3r = liin.out(3);
  if (mass.getmu() != liin.getmu()){
    cout << "ERROR: differing mu's in mpi4lnf2\n";}
  return pow(mass.getmpi()/mass.getfpi(),2)*2.*l3r;
}

double mpi4rnf2(const physmassnf2 mass){
  double mp2 = pow(mass.getmpi(),2);
  double mu2 = pow(mass.getmu(),2);
  return Ab(mp2,mu2) * (  - 1./2.)/pow(mass.getfpi(),2);
}
// fpi ******************************************************************

double fpi4nf2(const physmassnf2 mass, const libar libarin){
  return pi16*pow(mass.getmpi()/mass.getfpi(),2)*libarin.out(4);
}

double fpi4nf2(const physmassnf2 mass, const li liin){
  return fpi4lnf2(mass,liin)+fpi4rnf2(mass);
}

double fpi4lnf2(const physmassnf2 mass, const li liin){
  double l4r = liin.out(4);
  if (mass.getmu() != liin.getmu()){
    cout << "ERROR: differing mu's in fpi4lnf2\n";}
  return pow(mass.getmpi()/mass.getfpi(),2)*l4r;
}

double fpi4rnf2(const physmassnf2 mass){
  double mp2 = pow(mass.getmpi(),2);
  double mu2 = pow(mass.getmu(),2);
  return Ab(mp2,mu2)/pow(mass.getfpi(),2);
}
// qbarq ******************************************************************

double qqup4nf2(const physmassnf2 mass, const libar libarin){
  return pi16*pow(mass.getmpi()/mass.getfpi(),2)*
    (2.*libarin.out(8)-0.5*libarin.out(3));
}

double qqup4nf2(const physmassnf2 mass, const li liin){
  return qqup4lnf2(mass,liin)+qqup4rnf2(mass);
}

double qqup4lnf2(const physmassnf2 mass, const li liin){
  double mu = liin.getmu();
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi4lnf2\n";}
  return pow(mass.getmpi()/mass.getfpi(),2)*2.*(liin.out(3)+liin.out(8));
}
double qqup4rnf2(const physmassnf2 mass){
  double mp2 = pow(mass.getmpi(),2);
  double mu2 = pow(mass.getmu(),2);
  return  Ab(mp2,mu2) * ( 3./2. )/pow(mass.getfpi(),2);
}
