// jbdYlm.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015-2022 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON and jbnumlib are licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// defines the associated Legendre polynomials and the spherical harmonics
// P_l^m via recursion relations
// Ylm via Plm
#include<iostream>
#include<cstdlib>  // used for abs(m)
#include<cmath>
#include<complex>
#include "jbnumlib.h"

typedef std::complex<double> dcomplex;

// convention (-1)^m (Condon-Shortley) NOT included in Plm
//  as wikipedia mentions for quantum mechanics for Ylm,
// the (-1)^m is included in the Ylm definition instead
// calculation with with recursion relation
double jbdPlm(const int el, const int em, const double x){
  int m = abs(em);
  if (el < m) return 0.;
  // start with Pmm
  // (2m-1)!!
  double y = sqrt(1.-x*x);
  double mfacfac = 1.;
  for(int i = 2*m-1; i>1; i = i-2) mfacfac *= double(i); 
  double Pmm = mfacfac*pow(y,m);
  //if (m%2==1) Pmm = -Pmm;// extra - sign for Condon-Shortley if wanted
  double Pml=0.;
  if (m==el){
    Pml = Pmm;}
  else{
    double Pmm1 = x*double(2*m+1)*Pmm;
    if (m==(el-1)){
      Pml = Pmm1;}
    else{
      // now the recursion
      for(int l=m;  l < (el-1); l++){
	Pml = (double(2*l+3)*x*Pmm1-double(l+m+1)*Pmm)/double(l-m+2);
	Pmm = Pmm1;
	Pmm1 = Pml;
      }
    }
  }
  if(em>=0) return Pml;
  double lmfac = 1.;
  for(int i=el+m; i > (el-m) ; i--){
    lmfac *= double(i);
  }
  Pml /= lmfac;
  if (m%2==1) Pml = -Pml;
  return Pml;
}
// a number of explicit expressions
// (-1)^m removed from the Wikipedia expressions
double jbdP00(const double x){
  return 1.;
}
double jbdP1m1(const double x){
  return -sqrt(1.-x*x)/2.;
}
double jbdP10(const double x){
  return x;
}
double jbdP11(const double x){
  return sqrt(1.-x*x);
}
double jbdP30(const double x){
  return 0.5*(5.*pow(x,3)-3.*x);
}
double jbdP40(const double x){
  return 0.125*(35.*pow(x,4)-30.*pow(x,2)+3.);
}
double jbdP41(const double x){
  return 5./2.*(7.*pow(x,3)-3.*x)*sqrt(1.-x*x);
}
double jbdP4m1(const double x){
  return -1./20.*jbdP41(x);
}

// Spherical harmonics
// Wikipedia quantum mechanics convention
// uses general Plm defined above
dcomplex jbdYlm(const int el, const int em, const double theta, const double phi){
  int m=abs(em);
  double lmfac = 1.;
  for(int i=el+m; i > (el-m) ; i--){
    lmfac *= double(i);
  }
  double Ylm1 = sqrt(double(2*el+1)/(4.*M_PI*lmfac))*
    jbdPlm(el,m,cos(theta));
  if (m%2==1) Ylm1 = -Ylm1;// related to Condon-Shortley phases
  if(em >= 0) return dcomplex(Ylm1*cos(double(m)*phi),Ylm1*sin(double(m)*phi));
  if (m%2==1) Ylm1 = -Ylm1;
  return dcomplex(Ylm1*cos(double(m)*phi),-Ylm1*sin(double(m)*phi));
}

//PDG version which say Condon-Shortley conventions
// some explicit cases
dcomplex jbdY00(const double theta, const double phi){
  return dcomplex(sqrt(1./(4.*M_PI)));
}
dcomplex jbdY10(const double theta, const double phi){
  return dcomplex(sqrt(3./(4.*M_PI))*cos(theta));
}
dcomplex jbdY11(const double theta, const double phi){
  double Y11r = -sqrt(3./(8.*M_PI))*sin(theta);
  return dcomplex(Y11r*cos(phi),Y11r*sin(phi));
}
dcomplex jbdY1m1(const double theta, const double phi){
  double Y11r = sqrt(3./(8.*M_PI))*sin(theta);
  return dcomplex(Y11r*cos(phi),-Y11r*sin(phi));
}
dcomplex jbdY20(const double theta, const double phi){
  return dcomplex(sqrt(5./(4.*M_PI))*(1.5*pow(cos(theta),2)-0.5));
}
dcomplex jbdY21(const double theta, const double phi){
  double Y11r = -sqrt(15./(32.*M_PI))*sin(2.*theta);
  return dcomplex(Y11r*cos(phi),Y11r*sin(phi));
}
dcomplex jbdY22(const double theta, const double phi){
  double Y11r = sqrt(15./(32.*M_PI))*pow(sin(theta),2);
  return dcomplex(Y11r*cos(2.*phi),Y11r*sin(2.*phi));
}
dcomplex jbdY2m1(const double theta, const double phi){
  double Y11r = sqrt(15./(32.*M_PI))*sin(2.*theta);
  return dcomplex(Y11r*cos(phi),-Y11r*sin(phi));
}
dcomplex jbdY2m2(const double theta, const double phi){
  double Y11r = sqrt(15./(32.*M_PI))*pow(sin(theta),2);
  return dcomplex(Y11r*cos(2.*phi),-Y11r*sin(2.*phi));
}

// testing code, explicit expressions versus general case
//int main(void){
//  std::cout <<"some output\n";
//  for(int i=-10; i < 11; i++){
//    double x = 0.1*double(i);
//    std::cout <<"0 0: "<<jbdPlm(0,0,x)<<' '<<jbdP00(x)<<'\n';
//    std::cout <<"1 1: "<<jbdPlm(1,1,x)<<' '<<jbdP11(x)<<'\n';
//    std::cout <<"1 0: "<<jbdPlm(1,0,x)<<' '<<jbdP10(x)<<'\n';
//    std::cout <<"1-1: "<<jbdPlm(1,-1,x)<<' '<<jbdP1m1(x)<<'\n';
//    std::cout <<"3 0: "<<jbdPlm(3,0,x)<<' '<<jbdP30(x)<<'\n';
//    std::cout <<"4 0: "<<jbdPlm(4,0,x)<<' '<<jbdP40(x)<<'\n';
//    std::cout <<"4 1: "<<jbdPlm(4,1,x)<<' '<<jbdP41(x)<<'\n';
//    std::cout <<"4-1: "<<jbdPlm(4,-1,x)<<' '<<jbdP4m1(x)<<'\n';
//  }
//
//  double theta = 0.23;
//  double phi = 1.3;
//  std::cout <<"0 0: "<<jbdYlm(0,0,theta,phi) <<' '<<jbdY00(theta,phi)<<'\n';
//  std::cout <<"1 1: "<<jbdYlm(1,1,theta,phi) <<' '<<jbdY11(theta,phi)<<'\n';
//  std::cout <<"1 0: "<<jbdYlm(1,0,theta,phi) <<' '<<jbdY10(theta,phi)<<'\n';
//  std::cout <<"1-1: "<<jbdYlm(1,-1,theta,phi)<<' '<<jbdY1m1(theta,phi)<<'\n';
//  std::cout <<"2 2: "<<jbdYlm(2,2,theta,phi) <<' '<<jbdY22(theta,phi)<<'\n';
//  std::cout <<"2 1: "<<jbdYlm(2,1,theta,phi) <<' '<<jbdY21(theta,phi)<<'\n';
//  std::cout <<"2 0: "<<jbdYlm(2,0,theta,phi) <<' '<<jbdY20(theta,phi)<<'\n';
//  std::cout <<"2-1: "<<jbdYlm(2,-1,theta,phi)<<' '<<jbdY2m1(theta,phi)<<'\n';
//  std::cout <<"2-2: "<<jbdYlm(2,-2,theta,phi)<<' '<<jbdY2m2(theta,phi)<<'\n';
//  std::cout <<"0 0: "<<jbdYlm(0,0,theta,phi) /jbdY00(theta,phi)<<'\n';
//  std::cout <<"1 1: "<<jbdYlm(1,1,theta,phi) /jbdY11(theta,phi)<<'\n';
//  std::cout <<"1 0: "<<jbdYlm(1,0,theta,phi) /jbdY10(theta,phi)<<'\n';
//  std::cout <<"1-1: "<<jbdYlm(1,-1,theta,phi)/jbdY1m1(theta,phi)<<'\n';
//  std::cout <<"2 2: "<<jbdYlm(2,2,theta,phi) /jbdY22(theta,phi)<<'\n';
//  std::cout <<"2 1: "<<jbdYlm(2,1,theta,phi) /jbdY21(theta,phi)<<'\n';
//  std::cout <<"2 0: "<<jbdYlm(2,0,theta,phi) /jbdY20(theta,phi)<<'\n';
//  std::cout <<"2-1: "<<jbdYlm(2,-1,theta,phi)/jbdY2m1(theta,phi)<<'\n';
//  std::cout <<"2-2: "<<jbdYlm(2,-2,theta,phi)/jbdY2m2(theta,phi)<<'\n';
//
//  return 0;
//}
	       
