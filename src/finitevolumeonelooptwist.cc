// finitevolumeonelooptwist.cc is part of the 
// CHIRON ChPT program collection
// Copyright (C) 2015-2016 Johan Bijnens, v1.1
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


// definitions of the components are done in integer
// but beware leading zeroes, since they give octal interpretration
// only produces problems from the three index case on,
// here leading zeroes suppressed

// the accuracies can be set in the various namespaces
#include <cmath>
#include <iostream>
#include <vector>
#include "jbnumlib.h"

#include "finitevolumeonelooptwist.h"

#ifndef DINTEGRAL
#define DINTEGRAL jbdgauss2
#endif


const double pi = M_PI;
const double pi16 = 1./(16.*pi*pi);
const double pi2 = M_PI*M_PI;


///////// tadpoles with twist ////////////////////////////////////////////
namespace Abvtwisttspace{
  double releps=1e-10;
  double msq,xl;
  std::vector<double> theta(4);
  int nprop,ntype,ncomponent;
}

double Abvtwisttinternal(double x){
  using namespace Abvtwisttspace;
  const double alpha = 0.5;
  double lt = pow(x/(1.-x),alpha);
  double ovfac = alpha*pow(lt,-1.-1./alpha)/pow(1.-x,2);
  // ovfac is dlt/lt^2 for the dimensionless lambda
  switch(nprop){
  case 1:
    ovfac *= -4./(xl*xl);
    break;
  case 2:
    ovfac *= lt;
    break;
  case 3:
    ovfac *= -lt*(lt*xl*xl/4.)/2.;
    break;
  case 4:
    ovfac *= lt*pow(lt*xl*xl/4.,2)/6.;
    break;
  default:
    std::cout << "Bad nprop in Abvtwistt, nprop = "<<nprop<<std::endl;
    break;
  }
  ovfac *= exp(-lt*msq*xl*xl/4.);

  double q = exp(-1./lt);
  double tx=0.,ty=0.,tz=0.,part1=0.;

  switch(ntype){
  case 0:
    tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
    if (theta[2] == theta[1]) ty = tx;
    else ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
    if (theta[3] == theta[1]) tz = tx;
    else{
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
    }
    return ovfac*(tx*ty*tz-1.);
    break;
  case 2:
    switch(ncomponent){
    case 1:
      tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return ovfac/(M_PI*lt*xl)*tx*ty*tz;
    break;
    case 2:
      ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[3] == theta[1]) tz = tx;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return ovfac/(M_PI*lt*xl)*tx*ty*tz;
    case 3:
      tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[2] == theta[1]) ty = tx;
      else ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      return ovfac/(M_PI*lt*xl)*tx*ty*tz;
    case 0:
      return 0.;
      break;
    default:
    std::cout << "Bad ncomponent in Abvtwistt, nprop = "<<nprop
	      <<" ntype = "<<ntype<<" ncomponent = "<<ncomponent<<std::endl;
    break;
    }
    break;
  case 22:
    tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
    if (theta[2] == theta[1]) ty = tx;
    else ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
    if (theta[3] == theta[1]) tz = tx;
    else{
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
    }
    return -ovfac*(2./(lt*xl*xl))*(tx*ty*tz-1.);
    break;
  case 23:
    switch(ncomponent){
    case 11:
      tx = jbderiv2utheta3(-theta[1]/(2.*M_PI),q);
      ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return ovfac/pow(M_PI*xl*lt,2)*tx*ty*tz;
      break;
    case 22:
      ty = jbderiv2utheta3(-theta[2]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[3] == theta[1]) tz = tx;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return ovfac/pow(M_PI*xl*lt,2)*tx*ty*tz;
      break;
    case 33:
      tz = jbderiv2utheta3(-theta[3]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[2] == theta[1]) ty = tx;
      else ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      return ovfac/pow(M_PI*xl*lt,2)*tx*ty*tz;
      break;
    case 12:
    case 21:
      tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      if (theta[2] == theta[1]) ty = tx;
      else ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return ovfac/pow(M_PI*xl*lt,2)*tx*ty*tz;
      break;
    case 13:
    case 31:
      tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      if (theta[3] == theta[1]) tz = tx;
      else tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      return ovfac/pow(M_PI*xl*lt,2)*tx*ty*tz;
      break;
    case 23:
    case 32:
      ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      return ovfac/pow(M_PI*xl*lt,2)*tx*ty*tz;
      break;
    case  0:
    case  1:
    case 10:
    case  2:
    case 20:
    case  3:
    case 30:
      // note that since there is no 0,1,2,3 this works out as expected
      return 0.;
      break;
    default:
      std::cout << "Bad ncomponent in Abvtwistt, nprop = "<<nprop<<std::endl;
      break;
    }
    break;
  case 33:
    switch(ncomponent){
    case   0:
    case  11:
    case  12:
    case  13:
    case  21:
    case  22:
    case  23:
    case  31:
    case  32:
    case  33:
    case 101:
    case 102:
    case 103:
    case 201:
    case 202:
    case 203:
    case 301:
    case 302:
    case 303:
    case 110:
    case 120:
    case 130:
    case 210:
    case 220:
    case 230:
    case 310:
    case 320:
    case 330:
      return 0.;
      break;
    case   1:
    case  10:
    case 100:
      tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return -2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      break;
    case   2:
    case  20:
    case 200:
      ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[3] == theta[1]) tz = tx;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return -2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      break;
    case   3:
    case  30:
    case 300:
      tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[2] == theta[1]) ty = tx;
      else ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      return -2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      break;
    case 111:
      tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      part1 = 6.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      tx = jbderiv3utheta3(-theta[1]/(2.*M_PI),q);
      //ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      //tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 222:
      ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[3] == theta[1]) tz = tx;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      part1 = 6.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      ty = jbderiv3utheta3(-theta[2]/(2.*M_PI),q);
      //tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      //tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 333:
      tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[2] == theta[1]) ty = tx;
      else ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      part1 = 6.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      tz = jbderiv3utheta3(-theta[3]/(2.*M_PI),q);
      //tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      //tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 112:
    case 121:
    case 211:
      ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[3] == theta[1]) tz = tx;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      part1 = 2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      //ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      tx = jbderiv2utheta3(-theta[1]/(2.*M_PI),q);
      //tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 113:
    case 131:
    case 311:
      tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[2] == theta[1]) ty = tx;
      else ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      part1 = 2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      //tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      tx = jbderiv2utheta3(-theta[1]/(2.*M_PI),q);
      //ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 221:
    case 212:
    case 122:
      tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      part1 = 2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      //tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      ty = jbderiv2utheta3(-theta[2]/(2.*M_PI),q);
      //tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 223:
    case 232:
    case 322:
      tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[2] == theta[1]) ty = tx;
      else ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      part1 = 2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      //tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      ty = jbderiv2utheta3(-theta[2]/(2.*M_PI),q);
      //tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 331:
    case 313:
    case 133:
      tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      if (theta[3] == theta[2]) tz = ty;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      part1 = 2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      //tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      tz = jbderiv2utheta3(-theta[3]/(2.*M_PI),q);
      //ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 332:
    case 323:
    case 233:
      ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      tx = jbdtheta3(-theta[1]/(2.*M_PI),q);
      if (theta[3] == theta[1]) tz = tx;
      else tz = jbdtheta3(-theta[3]/(2.*M_PI),q);
      part1 = 2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tz;
      //tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      tz = jbderiv2utheta3(-theta[3]/(2.*M_PI),q);
      //ty = jbdtheta3(-theta[2]/(2.*M_PI),q);
      return part1 +ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    case 123:
    case 132:
    case 213:
    case 231:
    case 312:
    case 321:
      tx = jbderivutheta3(-theta[1]/(2.*M_PI),q);
      if (theta[2] == theta[1]) ty = tx;
      else ty = jbderivutheta3(-theta[2]/(2.*M_PI),q);
      if (theta[3] == theta[1]) tz = tx;
      else{
	if (theta[3] == theta[2]) tz = ty;
	tz = jbderivutheta3(-theta[3]/(2.*M_PI),q);
      }
      return ovfac/(pow(M_PI*lt*xl,3))*tx*ty*tz;
      break;
    default:
      std::cout << "Bad ncomponent in Abvtwistt, nprop = "<<nprop<<std::endl;
      break;
    }
    break;
    // YYYYYYYYYYYYYYYYYYYYYY
  default:
    std::cout << "Bad ntype in Abvtwistt, nprop = "<<nprop
	      << " ntype = "<<ntype<<std::endl;
    break;
  }
  return 0.;
}

double AbVtwistt(const int nprop,const int ntype,const int ncomponent,
		 const double msq,const double xl,
		 const fourvector theta){
  Abvtwisttspace::nprop = nprop;
  Abvtwisttspace::ntype = ntype;
  Abvtwisttspace::ncomponent = ncomponent;
  Abvtwisttspace::msq = msq;
  Abvtwisttspace::xl = xl;
  Abvtwisttspace::theta = theta;

  double result = DINTEGRAL(Abvtwisttinternal,0.,1.,Abvtwisttspace::releps);
  return pi16*result;
}

// some other abbreviations
double AbVtwistt(const int nprop, const double msq,const double xl,
		 const fourvector theta){
  return AbVtwistt(nprop, 0,0,msq,xl,theta);
}
double A22bVtwistt(const int nprop, const double msq,const double xl,
		 const fourvector theta){
  return AbVtwistt(nprop,22,0,msq,xl,theta);
}
double A23bVtwistt(const int nprop, const int mu1, const int mu2,
		   const double msq,const double xl,
		 const fourvector theta){
  int munu = 10*mu1+mu2;
  return AbVtwistt(nprop,23,munu,msq,xl,theta);
}
double A33bVtwistt(const int nprop, const int mu1, const int mu2,
		   const int mu3, const double msq,const double xl,
		 const fourvector theta){
  int munurho = 100*mu1+10*mu2+mu3;
  return AbVtwistt(nprop,33,munurho,msq,xl,theta);
}

////////////////// Bubbles with twist ///////////////////////////////////////
namespace Bbvtwisttspace{
  double m1sq,m2sq,qsq,xl,releps=1e-10;
  int nprop,ntype,ncomponent;
  std::vector<double> theta(4),qext(4);
}

double Bbvtwisttinternal(double xx[2]){
  using namespace Bbvtwisttspace;
  double x = xx[0];
  double y = 1.-x;
  double z = xx[1];
  const double alpha = 0.5;
  double lt = pow(z/(1.-z),alpha);
  double ovfac = alpha*pow(lt,-1./alpha)/pow(1.-z,2);// includes a 1/lambda

  //double xl2 = xl*xl;

  switch(nprop){
  case 1:
    ovfac *= pi16;
    break;
  case 2:
    ovfac *= -pi16*xl*xl/4.*x*lt;
    break;
  case 3:
    ovfac *= -pi16*xl*xl/4.*y*lt;
    break;
  case 4:
    ovfac *= pi16*pow(xl*xl/4.,2)*x*y*lt*lt;
    break;
  default:
    std::cout << "Bad nprop in Bbvtwistt, nprop = "<<nprop<<std::endl;
    break;
  }
  std::vector<double> thetap = {0.,theta[1]-y*qext[1]*xl,
				theta[2]-y*qext[2]*xl,theta[3]-y*qext[3]*xl};
  double tx=0.,ty=0.,tz=0.,txp=0.,typ=0.,tzp=0.,txpp=0.,typp=0.,tzpp=0.,
    txppp=0.,typpp=0.,tzppp=0.,part1=0.;
  double q = exp(-1./lt);
  double mhat2 = x*m1sq+y*m2sq-x*y*qsq;
  ovfac *= exp(-lt*mhat2*xl*xl/4.);
  switch(ntype){
  case 0:
    tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
    if (thetap[2] == thetap[1]) ty = tx;
    else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
    if (thetap[3] == thetap[1]) tz = tx;
    else{
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
    }
    return ovfac*(tx*ty*tz-1.);
    break;
  case 1:
    tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
    if (thetap[2] == thetap[1]) ty = tx;
    else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
    if (thetap[3] == thetap[1]) tz = tx;
    else{
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
    }
    return ovfac*y*(tx*ty*tz-1.);
    break;
  case 2:
    switch(ncomponent){
    case 1:
      tx = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);     
      return ovfac/(M_PI*lt*xl)*tx*ty*tz;
      break;
    case 2:
      ty = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);     
      return ovfac/(M_PI*lt*xl)*tx*ty*tz;
      break;
    case 3:
      tz = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);     
      return ovfac/(M_PI*lt*xl)*tx*ty*tz;
      break;
    case 0:
      return 0.;
      break;
    default:
      std::cout <<"Wrong component in BbVtwistt\n";
      return 0.;
      break;
    }
    break;
  case 21:
    tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
    if (thetap[2] == thetap[1]) ty = tx;
    else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
    if (thetap[3] == thetap[1]) tz = tx;
    else{
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
    }
    return ovfac*y*y*(tx*ty*tz-1.);
    break;
  case 22:
    tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
    if (thetap[2] == thetap[1]) ty = tx;
    else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
    if (thetap[3] == thetap[1]) tz = tx;
    else{
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
    }
    return -ovfac*2./(lt*xl*xl)*(tx*ty*tz-1.);
    break;
  case 23:
    switch(ncomponent){
    case 11:
      tx = jbderiv2utheta3(-thetap[1]/(2.*M_PI),q);
      ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      return ovfac*
	(2.*y/(M_PI*lt*xl)*qext[1]*txp*ty*tz+tx*ty*tz/pow(M_PI*lt*xl,2));
      break;
    case 22:
      ty = jbderiv2utheta3(-thetap[2]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      return ovfac*
	(2.*y/(M_PI*lt*xl)*qext[2]*tx*typ*tz+tx*ty*tz/pow(M_PI*lt*xl,2));
      break;
    case 33:
      tz = jbderiv2utheta3(-thetap[3]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      return ovfac*
	(2.*y/(M_PI*lt*xl)*qext[3]*tx*ty*tzp+tx*ty*tz/pow(M_PI*lt*xl,2));
      break;
    case 12:
    case 21:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
	if (thetap[3] == thetap[2]) tz = ty;
	else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp =  jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) typ = txp;
      else typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      return ovfac*
	(y/(M_PI*lt*xl)*(qext[1]*tx*typ+qext[2]*txp*ty)*tz
	 +txp*typ*tz/pow(M_PI*lt*xl,2));
      break;
    case 13:
    case 31:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
	if (thetap[3] == thetap[2]) tz = ty;
	else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp =  jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tzp = txp;
      else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      return ovfac*
	(y/(M_PI*lt*xl)*(qext[1]*tx*tzp+qext[3]*txp*tz)*ty
	 +txp*ty*tzp/pow(M_PI*lt*xl,2));
      break;
    case 23:
    case 32:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
	if (thetap[3] == thetap[2]) tz = ty;
	else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      typ =  jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tzp = typ;
      else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      return ovfac*
	(y/(M_PI*lt*xl)*(qext[2]*ty*tzp+qext[3]*typ*tz)*tx
	 +tx*typ*tzp/pow(M_PI*lt*xl,2));
      break;
    case  1:
    case 10:
      tx = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      return ovfac*y/(M_PI*lt*xl)*qext[0]*tx*ty*tz;
      break;
    case  2:
    case 20:
      ty = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      return ovfac*y/(M_PI*lt*xl)*qext[0]*tx*ty*tz;
      break;
    case  3:
    case 30:
      tz = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      return ovfac*y/(M_PI*lt*xl)*qext[0]*tx*ty*tz;
      break;
    case  0:
      return 0.;
      break;
    default:
      std::cout << "Bad ncomponent in Bbvtwistt, ntype = "<<ntype<<std::endl;
      break;
    }
  case 31:
    tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
    if (thetap[2] == thetap[1]) ty = tx;
    else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
    if (thetap[3] == thetap[1]) tz = tx;
    else{
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
    }
    return ovfac*y*y*y*(tx*ty*tz-1.);
    break;
  case 32:
    tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
    if (thetap[2] == thetap[1]) ty = tx;
    else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
    if (thetap[3] == thetap[1]) tz = tx;
    else{
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
    }
    return -ovfac*2.*y/(lt*xl*xl)*(tx*ty*tz-1.);
    break;
  case 33:
    switch(ncomponent){
    case 000:
      return 0.;
      break;
    case   1:
    case  10:
    case 100:
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      part1 =  -2.*ovfac/(M_PI*lt*lt*pow(xl,3))*txp*ty*tz;
      return part1+ovfac*y*y/(M_PI*lt*xl)*txp*ty*tz*pow(qext[0],2);
      break;
    case   2:
    case  20:
    case 200:
      typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      part1 = -2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*typ*tz;
      return part1+ovfac*y*y/(M_PI*lt*xl)*tx*typ*tz*pow(qext[0],2);
      break;
    case   3:
    case  30:
    case 300:
      tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      part1 = -2.*ovfac/(M_PI*lt*lt*pow(xl,3))*tx*ty*tzp;
      return part1+ovfac*y*y/(M_PI*lt*xl)*tx*ty*tzp*pow(qext[0],2);
      break;
    case  11:
    case 101:
    case 110:
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tz = ty;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      txpp = jbderiv2utheta3(-thetap[1]/(2.*M_PI),q);
      part1 = ovfac*y*y/(M_PI*lt*xl)*txp*ty*tz*2.*qext[0]*qext[1];
      return part1+ovfac*y/(pow(M_PI*lt*xl,2))*txpp*ty*tz*qext[0];
      break;
    case  22:
    case 202:
    case 220:
      typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      typp = jbderiv2utheta3(-thetap[2]/(2.*M_PI),q);
      part1 = ovfac*y*y/(M_PI*lt*xl)*tx*typ*tz*2.*qext[0]*qext[2];
      return part1+ovfac*y/(pow(M_PI*lt*xl,2))*tx*typp*tz*qext[0];
      break;
    case  33:
    case 303:
    case 330:
      tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      tzpp = jbderiv2utheta3(-thetap[3]/(2.*M_PI),q);
      part1 = ovfac*y*y/(M_PI*lt*xl)*tx*ty*tzp*2.*qext[0]*qext[3];
      return part1+ovfac*y/(pow(M_PI*lt*xl,2))*tx*ty*tzpp*qext[0];
      break;
    case  12:
    case  21:
    case 102:
    case 120:
    case 201:
    case 210:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) typ = txp;
      else typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      part1 = qext[0]*ovfac*(
	     y*y/(M_PI*lt*xl)*(txp*ty*tz*qext[2]+tx*typ*tz*qext[1])
	     +y/(pow(M_PI*lt*xl,2))*txp*typ*tz
			     );
      return part1;
      break;
    case  13:
    case  31:
    case 103:
    case 130:
    case 301:
    case 310:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tzp = txp;
      else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      part1 = qext[0]*ovfac*(
	     y*y/(M_PI*lt*xl)*(txp*ty*tz*qext[3]+tx*ty*tzp*qext[1])
	     +y/(pow(M_PI*lt*xl,2))*txp*ty*tzp
			     );
      return part1;
      break;
    case  23:
    case  32:
    case 203:
    case 230:
    case 302:
    case 320:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tzp = typ;
      else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      part1 = qext[0]*ovfac*(
	     y*y/(M_PI*lt*xl)*(tx*typ*tz*qext[3]+tx*ty*tzp*qext[2])
	     +y/(pow(M_PI*lt*xl,2))*tx*typ*tzp
			     );
      return part1;
      break;
    case 112:
    case 121:
    case 211:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) typ = txp;
      else typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      txpp = jbderiv2utheta3(-thetap[1]/(2.*M_PI),q);
      part1 = 	2./(M_PI*lt*lt*pow(xl,3))*tx*typ*tz
        +y*y/(M_PI*lt*xl)*(2.*txp*ty*tz*qext[1]*qext[2]+tx*typ*tz*qext[1]*qext[1])
	+y/(pow(M_PI*lt*xl,2))*(txpp*ty*tz*qext[2]+2.*txp*typ*tz*qext[1])
        +1/pow(M_PI*lt*xl,3)*txpp*typ*tz;
      return ovfac*part1;
      break;
    case 113:
    case 131:
    case 311:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tzp = txp;
      else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      txpp = jbderiv2utheta3(-thetap[1]/(2.*M_PI),q);
      part1 = 2./(M_PI*lt*lt*pow(xl,3))*tx*ty*tzp
       +y*y/(M_PI*lt*xl)*(2.*txp*ty*tz*qext[1]*qext[3]+tx*ty*tzp*qext[1]*qext[1])
	+y/(pow(M_PI*lt*xl,2))*(txpp*ty*tz*qext[3]+2.*txp*ty*tzp*qext[1])
        +1/pow(M_PI*lt*xl,3)*txpp*ty*tzp;
      return ovfac*part1;
      break;
    case 221:
    case 212:
    case 122:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) typ = txp;
      else typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      typp = jbderiv2utheta3(-thetap[2]/(2.*M_PI),q);
      part1 = 2./(M_PI*lt*lt*pow(xl,3))*txp*ty*tz
        +y*y/(M_PI*lt*xl)*(2.*tx*typ*tz*qext[1]*qext[2]+txp*ty*tz*qext[2]*qext[2])
	+y/(pow(M_PI*lt*xl,2))*(tx*typp*tz*qext[1]+2.*txp*typ*tz*qext[2])
        +1/pow(M_PI*lt*xl,3)*txp*typp*tz;
      return ovfac*part1;
      break;
    case 223:
    case 232:
    case 322:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tzp = typ;
      else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      typp = jbderiv2utheta3(-thetap[2]/(2.*M_PI),q);
      part1 =  2./(M_PI*lt*lt*pow(xl,3))*tx*ty*tzp
        +y*y/(M_PI*lt*xl)*(2.*tx*typ*tz*qext[2]*qext[3]+tx*ty*tzp*qext[2]*qext[2])
	+y/(pow(M_PI*lt*xl,2))*(tx*typp*tz*qext[3]+2.*tx*typ*tzp*qext[2])
        +1/pow(M_PI*lt*xl,3)*tx*typp*tzp;
      return ovfac*part1;
      break;
    case 331:
    case 313:
    case 133:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tzp = txp;
      else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      tzpp = jbderiv2utheta3(-thetap[3]/(2.*M_PI),q);
      part1 = 2./(M_PI*lt*lt*pow(xl,3))*txp*ty*tz
        +y*y/(M_PI*lt*xl)*(2.*tx*ty*tzp*qext[1]*qext[3]+txp*ty*tz*qext[3]*qext[3])
	+y/(pow(M_PI*lt*xl,2))*(tx*ty*tzpp*qext[1]+2.*txp*ty*tzp*qext[3])
        +1/pow(M_PI*lt*xl,3)*txp*ty*tzpp;
      return ovfac*part1;
      break;
    case 332:
    case 323:
    case 233:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[2]) tzp = typ;
      else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      tzpp = jbderiv2utheta3(-thetap[3]/(2.*M_PI),q);
      part1 = 2./(M_PI*lt*lt*pow(xl,3))*tx*typ*tz
        +y*y/(M_PI*lt*xl)*(2.*tx*ty*tzp*qext[2]*qext[3]+tx*typ*tz*qext[3]*qext[3])
	+y/(pow(M_PI*lt*xl,2))*(tx*ty*tzpp*qext[2]+2.*tx*typ*tzp*qext[3])
        +1/pow(M_PI*lt*xl,3)*tx*typ*tzpp;
      return ovfac*part1;
      break;
    case 123:
    case 132:
    case 213:
    case 231:
    case 312:
    case 321:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) typ = txp;
      else typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tzp = txp;
      else {
	if (thetap[3] == thetap[2]) tzp = typ;
	else tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      }
      part1 = 
        +y*y/(M_PI*lt*xl)*(txp*ty*tz*qext[2]*qext[3]+tx*typ*tz*qext[1]*qext[3]
                           +tx*ty*tzp*qext[1]*qext[2])
	+y/(pow(M_PI*lt*xl,2))*(tx*typ*tzp*qext[1]+txp*ty*tzp*qext[2]
                                +txp*typ*tz*qext[3])
        +1/pow(M_PI*lt*xl,3)*txp*typ*tzp;
      return ovfac*part1;
      break;
    case 111:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      txp = jbderivutheta3(-thetap[1]/(2.*M_PI),q);
      txpp = jbderiv2utheta3(-thetap[1]/(2.*M_PI),q);
      txppp = jbderiv3utheta3(-thetap[1]/(2.*M_PI),q);
      part1 =
 	6./(M_PI*lt*lt*pow(xl,3))*txp*ty*tz
        +y*y/(M_PI*lt*xl)*(3.*txp*ty*tz*qext[1]*qext[1])
	+y/(pow(M_PI*lt*xl,2))*(3.*txpp*ty*tz*qext[1])
        +1/pow(M_PI*lt*xl,3)*txppp*ty*tz;
      return ovfac*part1;
      break;
    case 222:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      typ = jbderivutheta3(-thetap[2]/(2.*M_PI),q);
      typp = jbderiv2utheta3(-thetap[2]/(2.*M_PI),q);
      typpp = jbderiv3utheta3(-thetap[2]/(2.*M_PI),q);
      part1 = 	
 	6./(M_PI*lt*lt*pow(xl,3))*tx*typ*tz
        +y*y/(M_PI*lt*xl)*(3.*tx*typ*tz*qext[2]*qext[2])
	+y/(pow(M_PI*lt*xl,2))*(3.*tx*typp*tz*qext[2])
        +1/pow(M_PI*lt*xl,3)*tx*typpp*tz;
      return ovfac*part1;
      break;
    case 333:
      tx = jbdtheta3(-thetap[1]/(2.*M_PI),q);
      if (thetap[2] == thetap[1]) ty = tx;
      else ty = jbdtheta3(-thetap[2]/(2.*M_PI),q);
      if (thetap[3] == thetap[1]) tz = tx;
      else{
      if (thetap[3] == thetap[2]) tz = ty;
      else  tz = jbdtheta3(-thetap[3]/(2.*M_PI),q);
      }
      tzp = jbderivutheta3(-thetap[3]/(2.*M_PI),q);
      tzpp = jbderiv2utheta3(-thetap[3]/(2.*M_PI),q);
      tzppp = jbderiv3utheta3(-thetap[3]/(2.*M_PI),q);
      part1 = 
 	6./(M_PI*lt*lt*pow(xl,3))*tx*ty*tzp
        +y*y/(M_PI*lt*xl)*(3.*tx*ty*tzp*qext[3]*qext[3])
	+y/(pow(M_PI*lt*xl,2))*(3.*tx*ty*tzpp*qext[3])
        +1/pow(M_PI*lt*xl,3)*tx*ty*tzppp;
      return ovfac*part1;
      break;

    default:
      std::cout << "Bad ncomponent in Bbvtwistt, ntype = "<<ntype<<std::endl;
      break;
    }
    break;
  default:
    std::cout << "Bad ntype in Bbvtwistt, ntype = "<<ntype<<std::endl;
    break;
  }
  return 0.;
}

double BbVtwistt(const int nprop,const int ntype,const int ncomponent,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  Bbvtwisttspace::m1sq = m1sq;
  Bbvtwisttspace::m2sq = m2sq;
  Bbvtwisttspace::qsq = qsq;
  Bbvtwisttspace::xl = xl;
  Bbvtwisttspace::nprop = nprop;
  Bbvtwisttspace::ntype = ntype;
  Bbvtwisttspace::ncomponent = ncomponent;
  Bbvtwisttspace::theta = theta;
  Bbvtwisttspace::qext = q;

  double a[2] = {0.,0.};
  double b[2] = {1.,1.};
  double releps = Bbvtwisttspace::releps;
  double relerr;
  int ifail;
  double result =  jbdad2(Bbvtwisttinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"precision not reached in Bbvtwistt, nprop, ntype, ncomponent ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ncomponent<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}

// shortcuts
double BbVtwistt(const int nprop,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  return BbVtwistt(nprop,0,0,m1sq,m2sq,qsq,xl,theta,q);
}
double B1bVtwistt(const int nprop,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  return BbVtwistt(nprop,1,0,m1sq,m2sq,qsq,xl,theta,q);
}
double B2bVtwistt(const int nprop,const int mu1,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  return BbVtwistt(nprop,2,mu1,m1sq,m2sq,qsq,xl,theta,q);
}
double B21bVtwistt(const int nprop,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  return BbVtwistt(nprop,21,0,m1sq,m2sq,qsq,xl,theta,q);
}
double B22bVtwistt(const int nprop,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  return BbVtwistt(nprop,22,0,m1sq,m2sq,qsq,xl,theta,q);
}
double B23bVtwistt(const int nprop,const int mu1, const int mu2,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  int munu = 10*mu1+mu2;
  return BbVtwistt(nprop,23,munu,m1sq,m2sq,qsq,xl,theta,q);
}
double B31bVtwistt(const int nprop,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  return BbVtwistt(nprop,31,0,m1sq,m2sq,qsq,xl,theta,q);
}
double B32bVtwistt(const int nprop,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  return BbVtwistt(nprop,32,0,m1sq,m2sq,qsq,xl,theta,q);
}
double B33bVtwistt(const int nprop,const int mu1, const int mu2, const int mu3,
		 const double m1sq, const double m2sq, const double qsq,
		 const double xl,
		 const fourvector theta, const fourvector q){
  int munurho = 100*mu1+10*mu2+mu3;
  return BbVtwistt(nprop,33,munurho,m1sq,m2sq,qsq,xl,theta,q);
}



void setprecisionAbVtwistt(const double releps){
  Abvtwisttspace::releps = releps;
}

void setprecisionBbVtwistt(const double releps){
  Bbvtwisttspace::releps = releps;
}

double getprecisionAbVtwistt(void){
  return Abvtwisttspace::releps;
}

double getprecisionBbVtwistt(void){
  return Bbvtwisttspace::releps;
}

void setprecisionfinitevolumeonelooptwistt(const double Aeps, const double Beps, const bool printout){
  Abvtwisttspace::releps = Aeps;
  Bbvtwisttspace::releps = Beps;
  if(printout){
    std::cout << "#accuracies AbVtwistt BbVtwistt : "
	      << Abvtwisttspace::releps <<' '<< Bbvtwisttspace::releps <<'\n';
  }
} 

void getprecisionfinitevolumeonelooptwistt(double & Aeps, double & Beps){
  Aeps = Abvtwisttspace::releps;
  Beps = Bbvtwisttspace::releps;
}
