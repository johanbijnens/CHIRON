// finitevolumeonelooptwist.h is part of the 
// CHIRON ChPT program collection
// Copyright (C) 2015-2016 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

//
// finitevolume one-loop integrals
// integral over jacobi theta function
// case with a general twist
// MINKOWSKI
// for the method and definitions see:
// J.~Bijnens and J.~Relefors, JHEP 1405 (2014) 015, [arXiv:1402.1385].


#ifndef FINITEVOLUMEONELOOPTWIST_H
#define FINITEVOLUMEONELOOPTWIST_H

#include "fourvector.h"

////////////// with theta functions //////////////////////////////////////

double AbVtwistt(const int nprop,const int ntype,const int ncomponent,
		 const double msq,const double xl,
		 const fourvector theta);

double BbVtwistt(const int nprop,const int ntype,const int ncomponent,
  	         const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q);


void setprecisionAbVtwistt(const double releps=1e-10);
void setprecisionBbVtwistt(const double releps=1e-10);
double getprecisionAbVtwistt(void);
double getprecisionBbVtwistt(void);
void setprecisionfinitevolumeonelooptwistt(const double Aeps=1e-10,
	  const double Beps=1e-10, const bool printout=true);

// shortcuts
double AbVtwistt(const int nprop, const double msq,const double xl,
		 const fourvector theta);
double A22bVtwistt(const int nprop, const double msq,const double xl,
		   const fourvector theta);
double A23bVtwistt(const int nprop, const int mu1, const int mu2,
		   const double msq,const double xl,
		   const fourvector theta);
double A33bVtwistt(const int nprop, const int mu1, const int mu2,
		   const int mu3, const double msq,const double xl,
		   const fourvector theta);
double BbVtwistt(const int nprop,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q);
double B1bVtwistt(const int nprop,
		  const double m1sq, const double m2sq, const double qsq,
		  const double xl,
		  const fourvector theta, const fourvector q);
double B2bVtwistt(const int nprop,const int mu1,
		  const double m1sq, const double m2sq, const double qsq,
		  const double xl,
		  const fourvector theta, const fourvector q);
double B21bVtwistt(const int nprop,
		   const double m1sq, const double m2sq, const double qsq,
		   const double xl,
		   const fourvector theta, const fourvector q);
double B22bVtwistt(const int nprop,
		   const double m1sq, const double m2sq, const double qsq,
		   const double xl,
		   const fourvector theta, const fourvector q);
double B23bVtwistt(const int nprop,const int mu1, const int mu2,
		   const double m1sq, const double m2sq, const double qsq,
		   const double xl,
		   const fourvector theta, const fourvector q);
double B31bVtwistt(const int nprop,
		   const double m1sq, const double m2sq, const double qsq,
		   const double xl,
		   const fourvector theta, const fourvector q);
double B32bVtwistt(const int nprop,
 		   const double m1sq, const double m2sq, const double qsq,
		   const double xl,
		   const fourvector theta, const fourvector q);
double B33bVtwistt(const int nprop,const int mu1, const int mu2, const int mu3,
		   const double m1sq, const double m2sq, const double qsq,
		   const double xl,
		   const fourvector theta, const fourvector q);

#endif // FINITEVOLUMEONELOOPTWIST_H
