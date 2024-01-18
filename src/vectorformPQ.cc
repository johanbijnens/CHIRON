// vectorformPQ.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Partially quenched vector form-factors

#include<cmath>
#include<iostream>
#include<complex>

#include "vectorformPQ.h"
#include "quenchedoneloopintegrals.h"

typedef std::complex<double> dcomplex;

//*************************************************************************
// auxiliary functions
// R^[m,n]_j({m_i},{\mu_k}) 
//   = \frac{\Pi_{l=1,n}(m_j-mu_l)/\frac{\Pi_{l=1,m; l \ne j }(m_j-m_l)
// my definition is with masses squared
// BUT COUNTING STARTS FROM ZERO

// the mass sets are done as vector here, so I can use size
// eliminates mistyping values for m,n
double Rmnj(const int j,const std::vector<double> Msq,
	    const std::vector<double> Musq){
  int m = int(Msq.size());
  int n = int(Musq.size());
  double mj = Msq[j];
  double rmnu = 1.,rmnd=1.;
  for(int k = 0; k < n; k++) rmnu *= (mj-Musq[k]);
  for(int k = 0; k < m; k++) if(j != k) rmnd *= (mj-Msq[k]);
  return rmnu/rmnd;
}

// derivative w.r.t. m^2_i of Rmnj
double Rimnj(const int i,const int j,const std::vector<double> Msq,
	     const std::vector<double> Musq){
  int m = int(Msq.size());
  int n = int(Musq.size());
  double mj = Msq[j];
  double rmnu = 1.,rmnd=1.;
  for(int k = 0; k < n; k++) rmnu *= (mj-Musq[k]);
  for(int k = 0; k < m; k++) if(j != k) rmnd *= (mj-Msq[k]);

  if( i != j) return rmnu/(rmnd*(mj-Msq[i]));

  double deriv = 0.;
  for(int k = 0; k < n; k++) deriv += 1./(mj-Musq[k]);
  for(int k = 0; k < m; k++) if(j != k) deriv += -1./(mj-Msq[k]);
  return rmnu/rmnd*deriv;
}

//***************************************************************************
// auxiliary functions for the mass ratios (old version)
double rn(const double m11,const double m22){
  return 1./(m11-m22);
}

double rs(const double m11,const double m22,const double m33){
  return (m11-m22)/(m11-m33);
}

double rs(const double m11,const double m22,const double m33,
	   const double m44,const double m55){
  return (m11-m22)*(m11-m33)/((m11-m44)*(m11-m55));
}

double rs(const double m11,const double m22,const double m33,
	  const double m44,const double m55, const double m66,
	  const double m77){
  return (m11-m22)*(m11-m33)*(m11-m44)/((m11-m55)*(m11-m66)*(m11-m77));
}

double rd(const double m11,const double m22,const double m33,
	   const double m44){
  return (m11-m22)*(m11-m33)/(m11-m44);
}

double rd(const double m11,const double m22,const double m33,
	   const double m44,const double m55, const double m66){
  return (m11-m22)*(m11-m33)*(m11-m44)/((m11-m55)*(m11-m66));
}

double fvpnf3p4L(const double qsq, const quarkmassnf mass,const Linf Liin){
  return 2.*qsq*Liin.out(9)/pow(mass.getf0(),2);
}


//**************************************************************************
dcomplex fvpv2s2nf3p4R(const double qsq, const quarkmassnf mass){
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 4) || (B0mq.size() != 4)){
    std::cout << "ERROR: fvpv2s2nf3p4R probably called with wrong number"
	 << "of quark masses\n";}
  // 2 valence and 2 sea quark masses
  double m11 = 2.*B0mq[0];
  double m33 = 2.*B0mq[1];
  double m44 = 2.*B0mq[2];
  double m66 = 2.*B0mq[3];

  double mu2 = mu0*mu0;
  // derived masses
  double m13,m14,m16,m34,m36,me2;
  me2 = (m44+2.*m66)/3.;
  m13 = (m11+m33)/2.;
  m14 = (m11+m44)/2.;
  m16 = (m11+m66)/2.;
  m34 = (m33+m44)/2.;
  m36 = (m33+m66)/2.;
  // version follwowong Bernard, Bijnens, Gamiz 2013 notation
  std::vector<double> Mu2  = {m44,m66};
  std::vector<double> M2xx = {m33,me2};
  std::vector<double> M2yy = {m11,me2};
  std::vector<double> M3xy = {m33,m11,me2};
  std::vector<double> Mxf  = {m34,m36};
  std::vector<double> Myf  = {m14,m16};
  std::vector<double> nsea = {2.,1.};// first sea mass twice

  dcomplex fp432 = 0.;

  // Note x = strange quark valence = 3
  //      y   up                    = 1
  // The spectator is down quark valence = 2

  // sum over sea quarks
  for(int k=0; k<int(Mxf.size()); k++)
    fp432 += 
      nsea[k]*(Ab(Mxf[k],mu2)+Ab(Myf[k],mu2)-4.*B22b(Mxf[k],Myf[k],qsq,mu2));

  // Dxx and Dyy terms, derivative on R^[m,n]_k
  for(int k=0; k<int(M2xx.size()); k++){
    fp432 += (-1./3.)*
      ( Rimnj(0,k,M2xx,Mu2)*(Ab(M2xx[k],mu2)-4.*B22b(m13,M2xx[k],qsq,mu2))
	+Rimnj(0,k,M2yy,Mu2)*(Ab(M2yy[k],mu2)-4.*B22b(m13,M2yy[k],qsq,mu2))
	);
  }

  // Dxx and Dyy terms
  // derivative on the integral, needs no sum, only first term contributes
  fp432 += (-1./3.)*
    (
     +Rmnj(0,M2xx,Mu2)*(Bb(M2xx[0],mu2)-4.*B22b(3,m13,M2xx[0],qsq,mu2))
     +Rmnj(0,M2yy,Mu2)*(Bb(M2yy[0],mu2)-4.*B22b(3,m13,M2yy[0],qsq,mu2))
     );
  // Dxy terms
  for(int k=0; k<int(M3xy.size()); k++)
    fp432 += (2./3.)*
      Rmnj(k,M3xy,Mu2)*(Ab(M3xy[k],mu2)-4.*B22b(m13,M3xy[k],qsq,mu2));

  return fp432/(4.*f0*f0);
}

//**************************************************************************
// partially quenched 23 case
dcomplex fvpv2s3nf3p4R(const double qsq, const quarkmassnf mass){
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 5) || (B0mq.size() != 5)){
    std::cout << "ERROR: fvpv2s3nf3p4R probably called with wrong number"
	 << "of quark masses\n";}
  // 2 valence and 3 sea quark masses
  double m11 = 2.*B0mq[0];
  double m33 = 2.*B0mq[1];
  double m44 = 2.*B0mq[2];
  double m55 = 2.*B0mq[3];
  double m66 = 2.*B0mq[4];

  double mu2 = mu0*mu0;

  //derived masses
  double m13,m14,m15,m16,m34,m35,m36,mp2,me2;
  // see notes for mp2,me2
  double a = (m44+m55)/2.;
  double b = (m44+m55+4.*m66)/6.;
  double c = (m44-m55)/sqrt(12.);
  mp2 = ((a+b)-sqrt(pow(a-b,2)+4.*c*c))/2.;
  me2 = ((a+b)+sqrt(pow(a-b,2)+4.*c*c))/2.;
  m13 = (m11+m33)/2.;
  m14 = (m11+m44)/2.;
  m15 = (m11+m55)/2.;
  m16 = (m11+m66)/2.;
  m34 = (m33+m44)/2.;
  m35 = (m33+m55)/2.;
  m36 = (m33+m66)/2.;

  std::vector<double> Mu3  = {m44,m55,m66};
  std::vector<double> M3xx = {m33,mp2,me2};
  std::vector<double> M3yy = {m11,mp2,me2};
  std::vector<double> M4xy = {m33,m11,mp2,me2};
  std::vector<double> Mxf  = {m34,m35,m36};
  std::vector<double> Myf  = {m14,m15,m16};

  dcomplex fp433 = 0.;

  // Note x = strange quark valence = 3
  //      y   up                    = 1
  // The spectator is down quark valence = 2

  // sum over sea quarks
  for(int k=0; k<int(Mxf.size()); k++)
    fp433 += Ab(Mxf[k],mu2)+Ab(Myf[k],mu2)-4.*B22b(Mxf[k],Myf[k],qsq,mu2);

  // Dxx and Dyy terms, derivative on R^[m,n]_k
  for(int k=0; k<int(M3xx.size()); k++){
    fp433 += (-1./3.)*
      ( Rimnj(0,k,M3xx,Mu3)*(Ab(M3xx[k],mu2)-4.*B22b(m13,M3xx[k],qsq,mu2))
	+Rimnj(0,k,M3yy,Mu3)*(Ab(M3yy[k],mu2)-4.*B22b(m13,M3yy[k],qsq,mu2))
	);
  }
  // Dxx and Dxy terms
  // derivative on the integral, needs no sum, only first term contributes
  fp433 += (-1./3.)*
    (
     +Rmnj(0,M3xx,Mu3)*(Bb(M3xx[0],mu2)-4.*B22b(3,m13,M3xx[0],qsq,mu2))
     +Rmnj(0,M3yy,Mu3)*(Bb(M3yy[0],mu2)-4.*B22b(3,m13,M3yy[0],qsq,mu2))
     );

  // Dxy terms
  for(int k=0; k<int(M4xy.size()); k++)
    fp433 += (2./3.)*
      Rmnj(k,M4xy,Mu3)*(Ab(M4xy[k],mu2)-4.*B22b(m13,M4xy[k],qsq,mu2));

  return fp433/(4.*f0*f0);

}

// f_+
// partially quenched 32 case
dcomplex fvpv3s2nf3p4R(const double qsq, const quarkmassnf mass){
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 5) || (B0mq.size() != 5)){
    std::cout << "ERROR: fvpv3s2nf3p4R probably called with wrong number"
	 << "of quark masses\n";}
  quarkmassnf mass2({B0mq[0],B0mq[2],B0mq[3],B0mq[4]},f0,mu0);
  return fvpv2s2nf3p4R(qsq,mass2);
}

// partially quenched 33 case
dcomplex fvpv3s3nf3p4R(const double qsq, const quarkmassnf mass){
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 6) || (B0mq.size() != 6)){
    std::cout << "ERROR: fvpv3s3nf3p4R probably called with wrong number"
	 << "of quark masses\n";}
  quarkmassnf mass2({B0mq[0],B0mq[2],B0mq[3],B0mq[4],B0mq[5]},f0,mu0);
  return fvpv2s3nf3p4R(qsq,mass2);
}


//************************************************************************
// alternative versions using the RS and RD functions more directly
// from PQ calculation using Bijnens, Danielsson,Lähde methods
dcomplex fvpv2s2nf3p4Rp(const double qsq, const quarkmassnf mass){
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 4) || (int(B0mq.size()) != 4)){
    std::cout << "ERROR: fvpv2s2nf3p4R probably called with wrong number"
	 << "of quark masses\n";}
  // 2 valence and 2 sea quark masses
  double m11 = 2.*B0mq[0];
  double m33 = 2.*B0mq[1];
  double m44 = 2.*B0mq[2];
  double m66 = 2.*B0mq[3];

  double mu2 = mu0*mu0;
  // derived masses
  double m13,m14,m16,m34,m36,me2;
  me2 = (m44+2.*m66)/3.;
  m13 = (m11+m33)/2.;
  m14 = (m11+m44)/2.;
  m16 = (m11+m66)/2.;
  m34 = (m33+m44)/2.;
  m36 = (m33+m66)/2.;

  dcomplex fpp4 =
    + Ab(me2,mu2) * (  - 1./12.*rs(me2,m44,m66,m11,m11)
  		       + 1./6.*rs(me2,m44,m66,m11,m33)
  		       - 1./12.*rs(me2,m44,m66,m33,m33) )
    + Ab(m11,mu2) * (  - 1./12.*rs(m11,m44,me2) 
  		       + 1./12.*rs(m11,m44,m66,me2,me2)
  		       + 1./6*rs(m11,m44,m66,m33,me2)
  		       - 1./12.*rs(m11,m66,me2) )
    + Ab(m14,mu2) * ( 1./2. )
    + Ab(m16,mu2) * ( 1./4. )
    + Ab(m33,mu2) * (  - 1./12.*rs(m33,m44,me2)
  		       + 1./12.*rs(m33,m44,m66,me2,me2)
  		       + 1./6.* rs(m33,m44,m66,m11,me2)
  		       - 1./12.*rs(m33,m66,me2) )
    + Ab(m34,mu2) * ( 1./2. )
    + Ab(m36,mu2) * ( 1./4. )
    + Bb(m11,mu2) * (  - 1./12.*rd(m11,m44,m66,me2) )
    + Bb(m33,mu2) * (  - 1./12.*rd(m33,m44,m66,me2) )
    + (B22b(m13,me2,qsq,mu2)) * ( 1./3.*rs(me2,m44,m66,m11,m11)
  				      - 2./3.*rs(me2,m44,m66,m11,m33)
  				      + 1./3.*rs(me2,m44,m66,m33,m33) )
    + (B22b(m13,m11,qsq,mu2)) * (  1./3.*rs(m11,m44,me2)
  				       - 1./3.*rs(m11,m44,m66,me2,me2)
  				       -2./3.*rs(m11,m44,m66,m33,me2) 
  				       + 1./3.*rs(m11,m66,me2) )
    + (B22b(m13,m33,qsq,mu2)) * (  1./3.*rs(m33,m44,me2)
  				       -1./3.*rs(m33,m44,m66,me2,me2)
  				       - 2/3.*rs(m33,m44,m66,m11,me2)
  				       + 1./3.*rs(m33,m66,me2) )             
    + (B22b(m14,m34,qsq,mu2)) * (  - 2. )
    + (B22b(m16,m36,qsq,mu2)) * (  - 1. )
    + (B22b(3,m13,m11,qsq,mu2)) * ( 1./3.*rd(m11,m44,m66,me2) )
    + (B22b(3,m13,m33,qsq,mu2)) * ( 1./3.*rd(m33,m44,m66,me2) );
    return fpp4/(f0*f0);
}

// 2 valence 3 sea masses
dcomplex fvpv2s3nf3p4Rp(const double qsq, const quarkmassnf mass){
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 5) || (B0mq.size() != 5)){
    std::cout << "ERROR: fvpv2s3nf3p4R probably called with wrong number"
	 << "of quark masses\n";}
  // 2 valence and 3 sea quark masses
  double m11 = 2.*B0mq[0];
  double m33 = 2.*B0mq[1];
  double m44 = 2.*B0mq[2];
  double m55 = 2.*B0mq[3];
  double m66 = 2.*B0mq[4];

  double mu2 = mu0*mu0;

  //derived masses
  double m13,m14,m15,m16,m34,m35,m36,mp2,me2;
  // see notes for mp2,me2
  double a = (m44+m55)/2.;
  double b = (m44+m55+4.*m66)/6.;
  double c = (m44-m55)/sqrt(12.);
  mp2 = ((a+b)-sqrt(pow(a-b,2)+4.*c*c))/2.;
  me2 = ((a+b)+sqrt(pow(a-b,2)+4.*c*c))/2.;
  m13 = (m11+m33)/2.;
  m14 = (m11+m44)/2.;
  m15 = (m11+m55)/2.;
  m16 = (m11+m66)/2.;
  m34 = (m33+m44)/2.;
  m35 = (m33+m55)/2.;
  m36 = (m33+m66)/2.;
  dcomplex FP4133221x23 =
       + Ab(mp2,mu2) * (  - 1./12.*rs(mp2,m44,m55,m66,m11,m11,me2) + 1./
         6.*rs(mp2,m44,m55,m66,m11,m33,me2) - 1./12.*rs(mp2,m44,m55,m66
         ,m33,m33,me2) );

      FP4133221x23 +=  + Ab(me2,mu2) * (  - 1./12.*rs(me2,m44,m55,m66,
         m11,m11,mp2) + 1./6.*rs(me2,m44,m55,m66,m11,m33,mp2) - 1./12.*
         rs(me2,m44,m55,m66,m33,m33,mp2) );

      FP4133221x23 +=  + Ab(m11,mu2) * (  - 1./12.*rs(m11,m44,m55,mp2,
         me2) + 1./12.*rs(m11,m44,m55,m66,mp2,mp2,me2) + 1./12.*rs(m11,
         m44,m55,m66,mp2,me2,me2) + 1./6.*rs(m11,m44,m55,m66,m33,mp2,
         me2) - 1./12.*rs(m11,m44,m66,mp2,me2) - 1./12.*rs(m11,m55,m66,
         mp2,me2) );

      FP4133221x23 +=  + Ab(m14,mu2) * ( 1./4. );

      FP4133221x23 +=  + Ab(m15,mu2) * ( 1./4. );

      FP4133221x23 +=  + Ab(m16,mu2) * ( 1./4. );

      FP4133221x23 +=  + Ab(m33,mu2) * (  - 1./12.*rs(m33,m44,m55,mp2,
         me2) + 1./12.*rs(m33,m44,m55,m66,mp2,mp2,me2) + 1./12.*rs(m33,
         m44,m55,m66,mp2,me2,me2) + 1./6.*rs(m33,m44,m55,m66,m11,mp2,
         me2) - 1./12.*rs(m33,m44,m66,mp2,me2) - 1./12.*rs(m33,m55,m66,
         mp2,me2) );

      FP4133221x23 +=  + Ab(m34,mu2) * ( 1./4. );

      FP4133221x23 +=  + Ab(m35,mu2) * ( 1./4. );

      FP4133221x23 +=  + Ab(m36,mu2) * ( 1./4. );

      FP4133221x23 +=  + B22b(m13,mp2,qsq,mu2) * ( 1./3.*rs(mp2,m44,m55
         ,m66,m11,m11,me2) - 2./3.*rs(mp2,m44,m55,m66,m11,m33,me2) + 1./
         3.*rs(mp2,m44,m55,m66,m33,m33,me2) );

      FP4133221x23 +=  + B22b(m13,me2,qsq,mu2) * ( 1./3.*rs(me2,m44,m55
         ,m66,m11,m11,mp2) - 2./3.*rs(me2,m44,m55,m66,m11,m33,mp2) + 1./
         3.*rs(me2,m44,m55,m66,m33,m33,mp2) );

      FP4133221x23 +=  + B22b(m13,m11,qsq,mu2) * ( 1./3.*rs(m11,m44,m55
         ,mp2,me2) - 1./3.*rs(m11,m44,m55,m66,mp2,mp2,me2) - 1./3.*rs(
         m11,m44,m55,m66,mp2,me2,me2) - 2./3.*rs(m11,m44,m55,m66,m33,
         mp2,me2) + 1./3.*rs(m11,m44,m66,mp2,me2) + 1./3.*rs(m11,m55,
         m66,mp2,me2) );

      FP4133221x23 +=  + B22b(m13,m33,qsq,mu2) * (  - 1./2. + 1./3.*rs(
         m33,m44,m55,mp2,me2) - 1./3.*rs(m33,m44,m55,m66,mp2,mp2,me2)
          - 1./3.*rs(m33,m44,m55,m66,mp2,me2,me2) - 2./3.*rs(m33,m44,
         m55,m66,m11,mp2,me2) + 1./3.*rs(m33,m44,m66,mp2,me2) + 1./3.*
         rs(m33,m55,m66,mp2,me2) );

      FP4133221x23 +=  + B22b(m14,m34,qsq,mu2) * (  - 1. );

      FP4133221x23 +=  + B22b(m15,m35,qsq,mu2) * (  - 1. );

      FP4133221x23 +=  + B22b(m16,m36,qsq,mu2) * (  - 1. );

      FP4133221x23 +=  + B22b(m33,m13,qsq,mu2) * ( 1./2. );

      FP4133221x23 +=  + Bb(m11,mu2) * (  - 1./12.*rd(m11,m44,m55,m66,
         mp2,me2) );

      FP4133221x23 +=  + Bb(m33,mu2) * (  - 1./12.*rd(m33,m44,m55,m66,
         mp2,me2) );

      FP4133221x23 +=  + B22b(3,m13,m11,qsq,mu2) * ( 1./3.*rd(m11,
         m44,m55,m66,mp2,me2) );

      FP4133221x23 +=  + B22b(3,m13,m33,qsq,mu2) * ( 1./3.*rd(m33,
         m44,m55,m66,mp2,me2) );

       return FP4133221x23/(f0*f0);
}
