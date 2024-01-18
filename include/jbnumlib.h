// jbnumlib.h is part of the numerical library jbdnumlib
// included with the CHIRON ChPT at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#ifndef JBNUMLIB_H
#define JBNUMLIB_H
#include <cmath>
#include <complex>

typedef std::complex<double> dcomplex; 

//**************************************************************************
//  integration routines
// eps is error but can be relative or absolute 
// real
double jbdgauss(double (*f)(const double),const double a,const double b,
		const double eps);
double jbdgauss(double (*f)(const double, void*),const double a,const double b,
		const double eps, void* ap);
double jbdgauss2(double (*f)(const double),const double a,const double b,
		const double eps);
double jbdgauss2(double (*f)(const double, void*),const double a,const double b,
		 const double eps, void* ap);
double jbdquad15(double (*f)(const double),const double a,const double b,
		const double eps);
double jbdquad15(double (*f)(const double, void*),const double a,const double b,
		 const double eps, void* ap);
double jbdquad21(double (*f)(const double),const double a,const double b,
		const double eps);
double jbdquad21(double (*f)(const double, void*),const double a,const double b,
		 const double eps, void* ap);

// real with singularity
double jbdcauch(double (*f)(const double),const double a,const double b,
		const double s,const double eps);
double jbdcauch(double (*f)(const double,void*),const double a,const double b,
		const double s,const double eps, void* ap);
double jbdcauch2(double (*f)(const double),const double a,const double b,
		const double s,const double eps);
double jbdcauch2(double (*f)(const double,void*),const double a,const double b,
		const double s,const double eps, void* ap);
double jbdsing15(double (*f)(const double),const double a,const double b,
		const double s,const double eps);
double jbdsing15(double (*f)(const double,void*),const double a,const double b,
		const double s,const double eps, void* ap);
double jbdsing21(double (*f)(const double),const double a,const double b,
		const double s,const double eps);
double jbdsing21(double (*f)(const double,void*),const double a,const double b,
		const double s,const double eps, void* ap);

// complex
dcomplex jbwgauss(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b,const double eps);
dcomplex jbwgauss(dcomplex (*f)(const dcomplex,void*),const dcomplex a,
		  const dcomplex b,const double eps, void* ap);
dcomplex jbwgauss2(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b,const double eps);
dcomplex jbwgauss2(dcomplex (*f)(const dcomplex,void*),const dcomplex a,
		   const dcomplex b,const double eps, void* ap);
dcomplex jbwquad15(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b,const double eps);
dcomplex jbwquad15(dcomplex (*f)(const dcomplex, void*),const dcomplex a,
		   const dcomplex b,const double eps, void* ap);
dcomplex jbwquad21(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b,const double eps);
dcomplex jbwquad21(dcomplex (*f)(const dcomplex, void*),const dcomplex a,
		   const dcomplex b,const double eps, void* ap);

// complex with singularity
//dcomplex jbwcauch2(dcomplex (*f)(const dcomplex),const dcomplex a,
//		 const dcomplex b,const dcomplex s,const double eps);

// dimension 2 and 3 real integration
double jbdad2(double (*fcn)(double x[]),double a[],double b[],
	      const double eps, double &relerr,int &ifail);
double jbdad3(double (*fcn)(double x[]),double a[],double b[],
	      const double eps, double &relerr,int &ifail);
double jbdad2(double (*fcn)(double x[], void*),void*, double a[],double b[],
	      const double eps, double &relerr,int &ifail);
double jbdad3(double (*fcn)(double x[], void*),void*, double a[],double b[],
	      const double eps, double &relerr,int &ifail);
double jbdad4(double (*fcn)(double x[], void*),void*, double a[],double b[],
	      const double eps, double &relerr,int &ifail);
// dimension 2 and 3 complex integration
dcomplex jbdad2(dcomplex (*fcn)(dcomplex x[],void*),void* ap,dcomplex a[],
		dcomplex b[],const double releps, double &relerr,int &ifail);
dcomplex jbdad3(dcomplex (*fcn)(dcomplex x[],void*),void* ap, dcomplex a[],
		dcomplex b[],const double releps, double &relerr,int &ifail);
dcomplex jbdad4(dcomplex (*fcn)(dcomplex x[],void*),void* ap,dcomplex a[],
		dcomplex b[],const double releps, double &relerr,int &ifail);
//**************************************************************************
// special functions
dcomplex jbdli2(const dcomplex x);
dcomplex jbdli2p(const dcomplex x);
dcomplex jbdli3(const dcomplex x);
dcomplex jbdli4(const dcomplex x);
double jbdbesi0(const double x); 
double jbdbesi1(const double x); 
double jbdbesk0(const double x); 
double jbdbesk1(const double x); 
// bessel functions implemented via recursion:
double jbdbesk2(const double x); 
double jbdbesk3(const double x); 
double jbdbesk4(const double x); 
double jbdtheta3(const double u,const double q);
double jbderivutheta3(const double u,const double q);
double jbderiv2utheta3(const double u,const double q);
double jbderiv3utheta3(const double u,const double q);
double jbdtheta3(const double u,const double q);
double jbdtheta30(const double q);
double jbdtheta30m1(const double q);
double jbdtheta32(const double q);
double jbdtheta34(const double q);
double jbdtheta2d0(const double alpha, const double beta, const double gamma);
double jbdtheta2d0m1(const double alpha, const double beta, const double gamma);
double jbdtheta2d02(const double alpha, const double beta, const double gamma);

double jbdPlm(const int el, const int em, const double x);
double jbdP00(const double x);
double jbdP1m1(const double x);
double jbdP10(const double x);
double jbdP11(const double x);
double jbdP30(const double x);
double jbdP40(const double x);
double jbdP41(const double x);
double jbdP4m1(const double x);
dcomplex jbdYlm(const int el, const int em, const double theta,
		const double phi);
dcomplex jbdY00(const double theta, const double phi);
dcomplex jbdY10(const double theta, const double phi);
dcomplex jbdY11(const double theta, const double phi);
dcomplex jbdY1m1(const double theta, const double phi);
dcomplex jbdY20(const double theta, const double phi);
dcomplex jbdY21(const double theta, const double phi);
dcomplex jbdY22(const double theta, const double phi);
dcomplex jbdY2m1(const double theta, const double phi);
dcomplex jbdY2m2(const double theta, const double phi);


//***************************************************************************
// roots and zeros
void jbdrteq3(const double r, const double s,const double t,
	      double x[3],double & d );

double jbdzerox(double (*f)(const double),const double a0,const double b0,
		const double eps,const int maxf=10000, const int mode=1);



#endif // JBNUMLIB_H
