// massdecayvevTV.cc is part of the 
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



#include<cmath>
#include "massdecayvevTV.h"
#include "finitevolumeonelooptwist.h"


//////////////// masses ////////////////////////////////////////////////////
double mpipp4TV(const lomass mass, const double L,
		const fourvector p, const fourvector thetau,
		const fourvector thetad, const fourvector thetas){
  double mpi2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-mpi2/3.;
  double fpi = mass.getf0();

  fourvector theta0;
  double part1 = mpi2/(fpi*fpi)*(-0.5*AbVtwistt(1,0,0,mpi2,L,theta0)
				 +1./6.*AbVtwistt(1,0,0,me2,L,theta0));
  double part2 = 0;
  for(int i=1; i<=3;i++){
    // remember the Minkowski signs!
    part2 += -p.out(i)/(fpi*fpi)*(
            -2.*AbVtwistt(1,2,i,mpi2,L,thetau-thetad)// pip
	    -AbVtwistt(1,2,i,mk2,L,thetau-thetas)    // kp
	    +AbVtwistt(1,2,i,mk2,L,thetad-thetas));   //ko
  }
  return part1+part2;
}


double mpiop4TV(const lomass mass, const double L,
		const fourvector p, const fourvector thetau,
		const fourvector thetad, const fourvector thetas){
  double mpi2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-mpi2/3.;
  double fpi = mass.getf0();

  fourvector theta0;
  double part1 = mpi2/(fpi*fpi)*(-AbVtwistt(1,0,0,mpi2,L,thetau-thetad)
				 +0.5*AbVtwistt(1,0,0,mpi2,L,theta0)
				 +1./6.*AbVtwistt(1,0,0,me2,L,theta0));
  return part1;
}

double mkpp4TV(const lomass mass, const double L,
	       const fourvector p, const fourvector thetau,
	       const fourvector thetad, const fourvector thetas){
  double mpi2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-mpi2/3.;
  double fpi = mass.getf0();

  fourvector theta0;
  double part1 = mk2/(fpi*fpi)*(-1./3.*AbVtwistt(1,0,0,me2,L,theta0));

  double part2 = 0;
  for(int i=1; i<=3;i++){
    // remember the Minkowski signs!
    part2 += -p.out(i)/(fpi*fpi)*(-AbVtwistt(1,2,i,mpi2,L,thetau-thetad)
				-2.*AbVtwistt(1,2,i,mk2,L,thetau-thetas)
				-AbVtwistt(1,2,i,mk2,L,thetad-thetas));
  }
  return part1+part2;
}

double mkop4TV(const lomass mass, const double L,
	       const fourvector p, const fourvector thetau,
	       const fourvector thetad, const fourvector thetas){
  double mpi2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-mpi2/3.;
  double fpi = mass.getf0();

  fourvector theta0;
  double part1 = mk2/(fpi*fpi)*(-1./3.*AbVtwistt(1,0,0,me2,L,theta0));

  double part2 = 0;
  for(int i=1; i<=3;i++){
    // remember the Minkowski signs!
    part2 += -p.out(i)/(fpi*fpi)*(+AbVtwistt(1,2,i,mpi2,L,thetau-thetad)
				-AbVtwistt(1,2,i,mk2,L,thetau-thetas)
				-2.*AbVtwistt(1,2,i,mk2,L,thetad-thetas));
  }
  return part1+part2;
}

double metap4TV(const lomass mass, const double L,
		const fourvector p, const fourvector thetau,
		const fourvector thetad, const fourvector thetas){
  double mpi2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-mpi2/3.;
  double fpi = mass.getf0();

  fourvector theta0;
  double part1 = -mk2/(fpi*fpi)*2./3.*( AbVtwistt(1,0,0,mk2,L,thetau-thetas)
				       +AbVtwistt(1,0,0,mk2,L,thetad-thetas))
         +(2./3.*me2-1./6.*mpi2)/(fpi*fpi)*( AbVtwistt(1,0,0,me2,L,theta0))
         +mpi2/(fpi*fpi)/6.*( 2.*AbVtwistt(1,0,0,mpi2,L,thetau-thetad)
			      +AbVtwistt(1,0,0,mpi2,L,theta0));
  return part1;
}


