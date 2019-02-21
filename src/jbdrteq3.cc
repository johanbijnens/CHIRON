// jbdrteq3.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// C++ version of the CERNLIB routine drteq3
// roots of third order equation
// x^3 + r x^2 + s x + t

#include <cmath>
#include <complex>
#include <iostream>

#include "jbnumlib.h"

typedef std::complex<double> dcomplex;

// this is the fortran sign function
inline double signfortran(const double a, const double b){
  if (b >= 0.) return fabs(a);
  return -fabs(a);
}

void jbdrteq3(const double r, const double s,const double t,
		       double x[3],double & d ){
  dcomplex z[3];
  double y[3];
  const double eps = 1e-6, delta = 1e-15;
  //const dcomplex i = dcomplex(0,1);
  const double zd1 = 1., zq1 = 1.,
    r1 = 2.*zd1/27., r2 = zd1/2., r3 = zd1/3.,
    w3 = 1.732050807568877294, r4 = w3/2.,
    q1 = 2.*zq1/27., q2 = zq1/2., q3 = zq1/3.;
  if( (s == 0.) && (t == 0.)){
    x[0] = -r;
    x[1]=0.;
    x[2]=0.;
    d=0.;
    return;
  }
  double p=s-r3*r*r;
  double q=(r1*r*r-r3*s)*r+t;
  d=pow(r2*q,2)+pow(r3*p,3);
  if(fabs(d) <= eps){
    double pp=s-q3*r*r;
    double qq=(q1*r*r-q3*s)*r+t;
    d=pow(q2*qq,2)+pow(q3*pp,3);
    p=pp;
    q=qq;
  }
  double h=r3*r;
  double h1=r2*q;
  if(d >=  delta){
    double h2=sqrt(d);
    double u0=-h1+h2;
    double v0=-h1-h2;
    double u=signfortran(pow(fabs(u0),r3),u0);
    double v=signfortran(pow(fabs(v0),r3),v0);
    x[0]=u+v-h;
    x[1]=-r2*(u+v)-h;
    x[2]=r4*fabs(u-v);
    // possible Newton improvement
    if((fabs(u0) <= eps) || (fabs(v0) <= eps)){
      y[0]=x[0];
      for(int k=0; k <2; k++){
	y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/((3.*y[k]+2.*r)*y[k]+s);
      }
      x[0]=y[2];
      z[0]= dcomplex(x[1],x[2]);
      for(int k = 0; k< 2; k++){
	z[k+1]=z[k]-(((z[k]+r)*z[k]+s)*z[k]+t)/((3.*z[k]+2.*r)*z[k]+s);
      }
      x[1]= real(z[2]);
      x[2]= imag(z[2]);
    }
    return;
  }
  if(fabs(d) <= delta){
    d=0.;
    double u=signfortran(pow(fabs(h1),r3),-h1);
    x[0]=u+u-h;
    x[1]=-u-h;
    x[2]=x[1];
    if(fabs(h1) <=  eps){
      y[0]=x[1];
      for(int k = 0; k < 2; k++){
        h1=(3*y[k]+2.*r)*y[k]+s;
        if(fabs(h1) > delta){
	  y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/h1;
	}
        else{
	  x[0]= -r3*r;
	  x[1]=x[0];
	  x[2]=x[0];
	  return;
	}
        x[0]=y[2];
        x[1]=-r2*(r+x[0]);
	x[2]=x[1];
      }
    }
    return;
  }
  double h3=sqrt(pow(fabs(r3*p),3));
  double h2=r3*acos(-h1/h3);
  h1=pow(h3,r3);
  double u=h1*cos(h2);
  double v=w3*h1*sin(h2);
  x[0]=u+u-h;
  x[1]=-u-v-h;
  x[2]=-u+v-h;
  if( (h3 <= eps) || (x[0] <= eps) || (x[1] <= eps) ||    
      (x[2] <= eps)){
    for(int j = 0; j< 3; j++){
      y[0]=x[j];
      for(int k = 0; k < 2; k++){
	y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/((3*y[k]+2*r)*y[k]+s);
      }
      x[j]=y[2];
    }
  }
}
