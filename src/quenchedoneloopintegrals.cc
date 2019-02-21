// quenchedoneloopintegrals.cc is part of the CHIRON ChPT at two loops program
// collection
// Copyright (C) 2016 Johan Bijnens, v1.00
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of the one-loop integrals

#include<cmath>
#include<iostream>
#include "quenchedoneloopintegrals.h"
#include "oneloopintegrals.h"

const double pi = M_PI;
const double pi16 = 1./(16.*pow(M_PI,2));

// bare bubble integral
dcomplex Bb(const int nprop, const double m1sq,const double m2sq,
	    const double qsq, const double mu2){
  double d;
  dcomplex b0intana;
  double a,b,c,dsq,y0,y1;

  switch(nprop){
  case 1:
      return Bb(m1sq,m2sq,qsq,mu2);
  case 2:
    if(qsq == 0.){
      if(m1sq == m2sq) return -pi16/(2.*m1sq);
      return -pi16/(m1sq-m2sq)*(1.-m2sq/(m1sq-m2sq)*log(m1sq/m2sq));
    }
    a=qsq;
    b=m1sq-m2sq-qsq;
    c=m2sq;
    dsq = (b*b-4.*c*a)/(4.*a*a);
    y0= b/(2.*a);
    y1= 1.+y0;
    b0intana = 0.5*log(fabs((y1*y1-dsq)/(y0*y0-dsq)));
    if (dsq >= 0){
      d=sqrt(dsq);
      b0intana += (0.5*y0/d)*(log(fabs((y1+d)/(y1-d)))
    			    -log(fabs((y0+d)/(y0-d))));
    }
    else{
      d=sqrt(-dsq);
      b0intana += (-y0/d)*(atan(y1/d)-atan(y0/d));
    }
    b0intana *= (-pi16/qsq);
    if (qsq >= pow(sqrt(m1sq)+sqrt(m2sq),2)){
      b0intana = b0intana
        +dcomplex(0.0,pi16*pi/qsq/sqrt(pow((qsq-m1sq-m2sq),2)-4.0*m1sq*m2sq)
    		*(m1sq-m2sq-qsq));
    }
    return b0intana;
  case 3:
    if(qsq == 0.){
      if(m1sq == m2sq) return -pi16/(2.*m2sq);
      return -pi16/(m2sq-m1sq)*(1.-m1sq/(m2sq-m1sq)*log(m2sq/m1sq));
    }
    a=qsq;
    b=m2sq-m1sq-qsq;
    c=m1sq;
    dsq = (b*b-4.*c*a)/(4.*a*a);
    y0= b/(2.*a);
    y1= 1.+y0;
    b0intana = 0.5*log(fabs((y1*y1-dsq)/(y0*y0-dsq)));
    if (dsq >= 0){
      d=sqrt(dsq);
      b0intana += (0.5*y0/d)*(log(fabs((y1+d)/(y1-d)))
    			    -log(fabs((y0+d)/(y0-d))));
    }
    else{
      d=sqrt(-dsq);
      b0intana += (-y0/d)*(atan(y1/d)-atan(y0/d));
    }
    b0intana *= (-pi16/qsq);
    if (qsq >= pow(sqrt(m1sq)+sqrt(m2sq),2)){
      b0intana = b0intana
        +dcomplex(0.0,pi16*pi/qsq/sqrt(pow((qsq-m1sq-m2sq),2)-4.0*m1sq*m2sq)
    		*(m2sq-m1sq-qsq));
    }
    return b0intana;
  }
  std::cout << "case nprop = "<<nprop<<" not implemented in Bb";
  return nan("");
}

dcomplex B1b(const int nprop, const double m1sq,const double m2sq,
	       const double qsq, const double mu2){
  switch(nprop){
  case 1:
    return B1b(m1sq,m2sq,qsq,mu2);
  case 2:
    if (qsq == 0.){
      if (m1sq == m2sq) return -pi16/(6.*m1sq);
      return pi16*((-pow(m1sq,2) + pow(m2sq,2) + 2*m1sq*m2sq*log(m1sq) - 
		    2*m1sq*m2sq*log(m2sq))/(2.*pow(m1sq - m2sq,3)));
    }
    return -0.5/qsq*(Bb(m1sq,mu2)-Bb(m1sq,m2sq,qsq,mu2)+
       (m2sq-m1sq-qsq)*Bb(2,m1sq,m2sq,qsq,mu2));
  case 3:
    if (qsq == 0.){
      if (m1sq == m2sq) return -pi16/(3.*m1sq);
      return pi16*((3*pow(m1sq,2) - 4*m1sq*m2sq + pow(m2sq,2) - 
      2*pow(m1sq,2)*log(m1sq) + 2*pow(m1sq,2)*log(m2sq))/
    (2.*pow(m1sq - m2sq,3)));
    }
    return -0.5/qsq*(-Bb(m2sq,mu2)+Bb(m1sq,m2sq,qsq,mu2)+
       (m2sq-m1sq-qsq)*Bb(3,m1sq,m2sq,qsq,mu2));
   }
  std::cout << "case nprop = "<<nprop<<" not implemented in B1b";
  return nan("");
}

dcomplex B21b(const int nprop, const double m1sq,const double m2sq,
		const double qsq, const double mu2){
  switch(nprop){
  case 1:
    return B21b(m1sq,m2sq,qsq,mu2);
  case 2:
    if(qsq == 0.){
      if (m1sq == m2sq) return -pi16/(12.*m1sq);
      return(-(2*pow(m1sq,3) + 3*pow(m1sq,2)*m2sq - 
       6*m1sq*pow(m2sq,2) + pow(m2sq,3) - 
       6*pow(m1sq,2)*m2sq*log(m1sq) + 6*pow(m1sq,2)*m2sq*log(m2sq))/
    (6.*pow(m1sq - m2sq,4)))*pi16;
    }
    return 1./qsq*(Bb(m1sq,m2sq,qsq,mu2)+m1sq*Bb(2,m1sq,m2sq,qsq,mu2)
        -4.*B22b(2,m1sq,m2sq,qsq,mu2)+pi16/2.);
  case 3:
    if(qsq == 0.){
      if (m1sq == m2sq) return -pi16/(4.*m1sq);
      return ((11*pow(m1sq,3) - 18*pow(m1sq,2)*m2sq + 9*m1sq*pow(m2sq,2) - 
      2*pow(m2sq,3) - 6*pow(m1sq,3)*log(m1sq) + 
      6*pow(m1sq,3)*log(m2sq))/(6.*pow(m1sq - m2sq,4)))*pi16;
    }
    return 1./qsq*(Bb(m2sq,mu2)+m1sq*Bb(3,m1sq,m2sq,qsq,mu2)
        -4.*B22b(3,m1sq,m2sq,qsq,mu2)+pi16/2.);
  }
  std::cout << "case nprop = "<<nprop<<" not implemented in B21b";
  return nan("");
}

dcomplex B22b(const int nprop, const double m1sq,const double m2sq,
		const double qsq, const double mu2){
  switch(nprop){
  case 1:
    return B22b(m1sq,m2sq,qsq,mu2);
  case 2:
    return 0.5*B1b(m2sq,m1sq,qsq,mu2);
  case 3:
    return 0.5*B1b(m1sq,m2sq,qsq,mu2);
  }
  std::cout << "case nprop = "<<nprop<<" not implemented in Bb";
  return nan("");
}

dcomplex B31b(const int nprop, const double m1sq,const double m2sq,
		const double qsq, const double mu2){
  switch(nprop){
  case 1:
    return B31b(m1sq,m2sq,qsq,mu2);
  case 2:
    if(qsq == 0.){
      if (m1sq == m2sq) return -pi16/(20.*m1sq);
      return ((-3*pow(m1sq,4) - 10*pow(m1sq,3)*m2sq + 
      18*pow(m1sq,2)*pow(m2sq,2) - 6*m1sq*pow(m2sq,3) + pow(m2sq,4) + 
      12*pow(m1sq,3)*m2sq*log(m1sq) - 12*pow(m1sq,3)*m2sq*log(m2sq))/
    (12.*pow(m1sq - m2sq,5)))*pi16;
    }
    return 1./qsq*(-0.5*B1b(m1sq,m2sq,qsq,mu2)
                  -0.5*m1sq*B1b(2,m1sq,m2sq,qsq,mu2)
                  +3./4.*B21b(m1sq,m2sq,qsq,mu2)
		   -3./4.*(m2sq-m1sq-qsq)*B21b(2,m1sq,m2sq,qsq,mu2)
                  -pi16/12.);
  case 3:
    if(qsq == 0.){
      if (m1sq == m2sq) return -pi16/(5.*m1sq);
      return ((25*pow(m1sq,4) - 48*pow(m1sq,3)*m2sq + 
      36*pow(m1sq,2)*pow(m2sq,2) - 16*m1sq*pow(m2sq,3) + 
      3*pow(m2sq,4) - 12*pow(m1sq,4)*log(m1sq) + 
      12*pow(m1sq,4)*log(m2sq))/(12.*pow(m1sq - m2sq,5)))*pi16;
    }
    return 1./qsq*(0.25*Bb(m2sq,mu2)-0.5*m1sq*B1b(3,m1sq,m2sq,qsq,mu2)
                   -3./4.*B21b(m1sq,m2sq,qsq,mu2)
		   -3./4.*(m2sq-m1sq-qsq)*B21b(3,m1sq,m2sq,qsq,mu2)
                  -pi16/6.);
  }
  std::cout << "case nprop = "<<nprop<<" not implemented in B31b";
  return nan("");
}

dcomplex B32b(const int nprop, const double m1sq,const double m2sq,
		const double qsq, const double mu2){
  switch(nprop){
  case 1:
    return B32b(m1sq,m2sq,qsq,mu2);
  case 2:
    return 0.25*B1b(m1sq,m2sq,qsq,mu2)
           +0.25*m1sq*B1b(2,m1sq,m2sq,qsq,mu2)
           -1./8.*B21b(m1sq,m2sq,qsq,mu2)
	   +1./8.*(m2sq-m1sq-qsq)*B21b(2,m1sq,m2sq,qsq,mu2)
           +pi16/24.;
  case 3:
    return +1./8.*Bb(m2sq,mu2)+0.25*m1sq*B1b(3,m1sq,m2sq,qsq,mu2)
      +1./8.*B21b(m1sq,m2sq,qsq,mu2)
	+1./8.*(m2sq-m1sq-qsq)*B21b(3,m1sq,m2sq,qsq,mu2)
        +pi16/12.;
  }
  std::cout << "case nprop = "<<nprop<<" not implemented in Bb";
  return nan("");
}
