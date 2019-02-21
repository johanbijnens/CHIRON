// massdecayvevTV.h is part of the 
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

#ifndef MASSDECAYVEVTV_H
#define MASSDECAYVEVTV_H

#include "inputs.h"
#include "fourvector.h"

double mpipp4TV(const lomass mass, const double L,
		const fourvector p, const fourvector thetau,
		const fourvector thetad, const fourvector thetas);
double mpiop4TV(const lomass mass, const double L,
		const fourvector p, const fourvector thetau,
		const fourvector thetad, const fourvector thetas);
double mkpp4TV(const lomass mass, const double L,
	       const fourvector p, const fourvector thetau,
	       const fourvector thetad, const fourvector thetas);
double mkop4TV(const lomass mass, const double L,
	       const fourvector p, const fourvector thetau,
	       const fourvector thetad, const fourvector thetas);
double metap4TV(const lomass mass, const double L,
		const fourvector p, const fourvector thetau,
		const fourvector thetad, const fourvector thetas);


#endif // MASSDECAYVEVTV_H
