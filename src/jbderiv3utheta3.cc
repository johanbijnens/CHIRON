// jbderiv3utheta3.cc is part of the numerical library jbnumlib
// included with the CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include <cmath>
#include <iostream>
#include "jbnumlib.h"

double jbderiv3utheta3(const double u,const double q)
  /*
     computes -16pi^3\sum_{n=1,\infty} n^3 q^(n^2) sin(2n\pi u)

     for q<= 0.5 simply sums the above series, up to n=10 terms

     for q >= 0.5 uses the deivative w.r.t. u of
                theta3(u,q) = sqrt(pi/|lnq|)*exp(-pi^2 u^2/|lnq|)
	*theta3(-iu pi/|lnq|,exp(-pi^2/|lnq|))

     precision works to machine precision (tested by comparing it against
     running with much higher nmaxu and nmaxl and comparing the two
     algorithms against each other) for reasonable u
  */
{
  const int nmaxl = 10, nmaxu = 2; 
  double a[nmaxl+1];// should be the larger of nmaxl/nmaxu
  if( (q < 0.) || (q >=1.)){
    std::cout << "not allowed input value in jbderivtheta3" << std::endl;
    return 0.;}
  if( q < 0.5 ){
    double qsq = q*q;
    double atemp = 1.;
    double p = q;
    for( int i = 0; i<= nmaxl; i++){
      a[i] = atemp;
      atemp *= p;
      p *= qsq;
    }
    double result = 0.;
    // Clenshaw resummation not used u small happens too often
    for (int i=nmaxl; i>0;i--){
      result += a[i]*double(i*i*i)*sin(2.*M_PI*u*double(i));
    }
    return 16.*M_PI*M_PI*M_PI*result;
  }
  // case >= 0.5
  double lam1 = M_PI/fabs(log(q));
  double lam2 = M_PI*lam1;
  // everything put in the exponentials otherwise overflow
  double result = 0.;
  for (int i=nmaxu; i>0;i--){
    double id = double(i);
    result += 
      (  (u+id)*(3.-2.*lam2*pow(u+id,2))*exp(lam2*(-2.*id*u-u*u-id*id))
        +(u-id)*(3.-2.*lam2*pow(u-id,2))*exp(lam2*( 2.*id*u-u*u-id*id)));
  }
  return sqrt(lam1)*lam2*lam2*4.*(u*(3.-2.*u*u*lam2)*exp(-u*u*lam2)
			       +result);// first term is n=0
}
