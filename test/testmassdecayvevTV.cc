// testmassdecayvevTV.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// implements the formulas for masses etc. for twisted boundary conditions of
//  J.~Bijnens and J.~Relefors,
//  ``Masses, Decay Constants and Electromagnetic Form-factors with
//  Twisted Boundary Conditions,''
//  JHEP {\bf 1405} (2014) 015
//  doi:10.1007/JHEP05(2014)015
//  [arXiv:1402.1385 [hep-lat]].


#include<iostream>
#include<iomanip>
#include<cmath>
#include "massdecayvevTV.h"
#include "massdecayvevloV.h"
#include "finitevolumeonelooptwist.h"
#include "finitevolumeoneloopintegrals.h"

int main(void){
  // setting precision finitevolumeintegrals
  // Note only first one relevant here, and you might want to play with it
  // a bit
  setprecisionfinitevolumeonelooptwistt(1e-10,1e-10,true);
  setprecisionfinitevolumeoneloopt(1e-10,1e-10,true);

  // inputs used in Bijnens-Relefors
  // remeber the plots in that paper are with RELATIVE corrections
  double mp0 = 0.1395;
  double mk0 = 0.495;
  double F0 = 0.0922;
  double mu = 0.77;
  lomass mass1(mp0,mk0,F0,mu);
  fourvector pp,thetau,thetad,thetas;

  double mp2 = mp0*mp0;
  double mk2 = mp0*mp0;
  double me2 = 4./3.*mk2-1./3.*mp2;
 
  cout << "inputs used\n";
  cout << mass1;

  //**************************************************************************
  cout << "\nNo twist and center of mass: reducing to older result (relative correction)\n";
  cout << setprecision(7);

  double mpiL = 2.;
  double L = mpiL/mp0;
  double dmpi = mpi4loVt(mass1,L);
  double dmpip = mpipp4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "particle  mp0*L  older result   programs with twist\n";
  cout << "pi+         "<< mpiL <<' '
       <<setw(15)<<dmpi/mp2 << "    " << dmpip/mp2 << '\n';
  mpiL = 4.;
  L = mpiL/mp0;
  dmpi = mpi4loVt(mass1,L);
  dmpip = mpipp4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "pi+         "<< mpiL <<' '
       <<setw(15)<<dmpi/mp2 << "    " << dmpip/mp2 << '\n';

  mpiL = 2.;
  L = mpiL/mp0;
  dmpi = mpi4loVt(mass1,L);
  dmpip = mpiop4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "pi0         "<< mpiL <<' '
       <<setw(15)<<dmpi/mp2 << "    " << dmpip/mp2 << '\n';
  mpiL = 4.;
  L = mpiL/mp0;
  dmpi = mpi4loVt(mass1,L);
  dmpip = mpiop4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "pi0         "<< mpiL <<' '
       <<setw(15)<<dmpi/mp2 << "    " << dmpip/mp2 << '\n';

  mpiL = 2.;
  L = mpiL/mp0;
  dmpi = mk4loVt(mass1,L);
  dmpip = mkpp4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "K+          "<< mpiL <<' '
       <<setw(15)<<dmpi/mk2 << "    " << dmpip/mk2 << '\n';
  mpiL = 4.;
  L = mpiL/mp0;
  dmpi = mk4loVt(mass1,L);
  dmpip = mkpp4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "K+          "<< mpiL <<' '
       <<setw(15)<<dmpi/mk2 << "    " << dmpip/mk2 << '\n';

  mpiL = 2.;
  L = mpiL/mp0;
  dmpi = mk4loVt(mass1,L);
  dmpip = mkop4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "K0          "<< mpiL <<' '
       <<setw(15)<<dmpi/mk2 << "    " << dmpip/mk2 << '\n';
  mpiL = 4.;
  L = mpiL/mp0;
  dmpi = mk4loVt(mass1,L);
  dmpip = mkop4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "K0          "<< mpiL <<' '
       <<setw(15)<<dmpi/mk2 << "    " << dmpip/mk2 << '\n';

  mpiL = 2.;
  L = mpiL/mp0;
  dmpi = meta4loVt(mass1,L);
  dmpip = metap4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "eta         "<< mpiL <<' '
       <<setw(15)<<dmpi/me2 << "    " << dmpip/me2 << '\n';
  mpiL = 4.;
  L = mpiL/mp0;
  dmpi = meta4loVt(mass1,L);
  dmpip = metap4TV(mass1,L,pp,thetau,thetad,thetas);
  cout << "eta         "<< mpiL <<' '
       <<setw(15)<<dmpi/me2 << "    " << dmpip/me2 << '\n';

  //************************************************************************
  cout << "\nNo twist but moving frames\n";
  fourvector pp1,pp2,pp3,pp4;
  pp1 = fourvector(0.,0.,0.,0.);
  pp2 = fourvector(0.,2.*M_PI/L,0.,0.);
  pp3 = fourvector(0.,2.*M_PI/L,2.*M_PI/L,0.);
  pp4 = fourvector(0.,4.*M_PI/L,0.,0.);
  cout << "particle  mp0*L      p_i=0       p_x = 2pi/L     p_x=p_y=2pi/L   p_x=4pi/L\n";
  double dm1,dm2,dm3,dm4,dm5;
  mpiL = 2.;
  L = mpiL/mp0;
  dm1 = mpipp4TV(mass1,L,pp1,thetau,thetad,thetas);
  dm2 = mpipp4TV(mass1,L,pp2,thetau,thetad,thetas);
  dm3 = mpipp4TV(mass1,L,pp3,thetau,thetad,thetas);
  dm4 = mpipp4TV(mass1,L,pp4,thetau,thetad,thetas);
  cout << "pi+         "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3<<setw(16)<<dm4<<'\n';
  dm1 = mpiop4TV(mass1,L,pp1,thetau,thetad,thetas);
  dm2 = mpiop4TV(mass1,L,pp2,thetau,thetad,thetas);
  dm3 = mpiop4TV(mass1,L,pp3,thetau,thetad,thetas);
  dm4 = mpiop4TV(mass1,L,pp4,thetau,thetad,thetas);
  cout << "pi0         "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3<<setw(16)<<dm4<<'\n';
  dm1 = mkpp4TV(mass1,L,pp1,thetau,thetad,thetas);
  dm2 = mkpp4TV(mass1,L,pp2,thetau,thetad,thetas);
  dm3 = mkpp4TV(mass1,L,pp3,thetau,thetad,thetas);
  dm4 = mkpp4TV(mass1,L,pp4,thetau,thetad,thetas);
  cout << "k+          "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3<<setw(16)<<dm4<<'\n';
  dm1 = mkop4TV(mass1,L,pp1,thetau,thetad,thetas);
  dm2 = mkop4TV(mass1,L,pp2,thetau,thetad,thetas);
  dm3 = mkop4TV(mass1,L,pp3,thetau,thetad,thetas);
  dm4 = mkop4TV(mass1,L,pp4,thetau,thetad,thetas);
  cout << "k0          "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3<<setw(16)<<dm4<<'\n';
  dm1 = metap4TV(mass1,L,pp1,thetau,thetad,thetas);
  dm2 = metap4TV(mass1,L,pp2,thetau,thetad,thetas);
  dm3 = metap4TV(mass1,L,pp3,thetau,thetad,thetas);
  dm4 = metap4TV(mass1,L,pp4,thetau,thetad,thetas);
  cout << "eta         "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3<<setw(16)<<dm4<<'\n';

  //************************************************************************
  cout << "\nTwist and extra momentum dp_i (added to needed twist part)\n";
  thetau = fourvector(0.,M_PI/4.,0.,0.);
  cout << "thetau = " << thetau << "\nthetad=thetas=0\n";
  cout << "particle  mp0*L     dp_i=0       dp_x = 2pi/L    dp_x=p_y=2pi/L   dp_x=4pi/L     dp_x = -2pi/L\n";
  mpiL = 2.;
  L = mpiL/mp0;
  dm1 = mpipp4TV(mass1,L,pp1+thetau/L,thetau,thetad,thetas);
  dm2 = mpipp4TV(mass1,L,pp2+thetau/L,thetau,thetad,thetas);
  dm3 = mpipp4TV(mass1,L,pp3+thetau/L,thetau,thetad,thetas);
  dm4 = mpipp4TV(mass1,L,pp4+thetau/L,thetau,thetad,thetas);
  dm5 = mpipp4TV(mass1,L,-pp2+thetau/L,thetau,thetad,thetas);
  cout << "pi+         "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3
       <<setw(16)<<dm4<<setw(16)<<dm5<<'\n';
  dm1 = mpiop4TV(mass1,L,pp1,thetau,thetad,thetas);
  dm2 = mpiop4TV(mass1,L,pp2,thetau,thetad,thetas);
  dm3 = mpiop4TV(mass1,L,pp3,thetau,thetad,thetas);
  dm4 = mpiop4TV(mass1,L,pp4,thetau,thetad,thetas);
  dm5 = mpiop4TV(mass1,L,-pp2,thetau,thetad,thetas);
  cout << "pi0         "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3
       <<setw(16)<<dm4<<setw(16)<<dm5<<'\n';
  dm1 = mkpp4TV(mass1,L,pp1+thetau/L,thetau,thetad,thetas);
  dm2 = mkpp4TV(mass1,L,pp2+thetau/L,thetau,thetad,thetas);
  dm3 = mkpp4TV(mass1,L,pp3+thetau/L,thetau,thetad,thetas);
  dm4 = mkpp4TV(mass1,L,pp4+thetau/L,thetau,thetad,thetas);
  dm5 = mkpp4TV(mass1,L,-pp2+thetau/L,thetau,thetad,thetas);
  cout << "k+          "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3
       <<setw(16)<<dm4<<setw(16)<<dm5<<'\n';
  dm1 = mkop4TV(mass1,L,pp1,thetau,thetad,thetas);
  dm2 = mkop4TV(mass1,L,pp2,thetau,thetad,thetas);
  dm3 = mkop4TV(mass1,L,pp3,thetau,thetad,thetas);
  dm4 = mkop4TV(mass1,L,pp4,thetau,thetad,thetas);
  dm5 = mkop4TV(mass1,L,-pp2,thetau,thetad,thetas);
  cout << "k0          "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3
       <<setw(16)<<dm4<<setw(16)<<dm5<<'\n';
  dm1 = metap4TV(mass1,L,pp1,thetau,thetad,thetas);
  dm2 = metap4TV(mass1,L,pp2,thetau,thetad,thetas);
  dm3 = metap4TV(mass1,L,pp3,thetau,thetad,thetas);
  dm4 = metap4TV(mass1,L,pp4,thetau,thetad,thetas);
  dm5 = metap4TV(mass1,L,-pp2,thetau,thetad,thetas);
  cout << "eta         "<< mpiL
       <<setw(16)<<dm1<<setw(16)<<dm2<<setw(16)<<dm3
       <<setw(16)<<dm4<<setw(16)<<dm5<<'\n';
}


