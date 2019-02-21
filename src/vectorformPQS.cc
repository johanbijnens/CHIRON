// vectorformPQS.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Partially quenched vector form-factors staggered case

#include<cmath>
#include<iostream>
#include<complex>

#include "vectorformPQS.h"
#include "quenchedoneloopintegrals.h"
#include "jbnumlib.h"

typedef std::complex<double> dcomplex;

//*************************************************************************
// auxiliary functions
// R^[m,n]_j({m_i},{\mu_k}) 
//   = \frac{\Pi_{l=1,n}(m_j-mu_l)/\frac{\Pi_{l=1,m; l \ne j }(m_j-m_l)
// my definition is with masses squared
// BUT COUNTING STARTS FROM ZERO
// extra s added (to avoid conflicts with version in vectorformPQ.cc)

// these are defined in vectorformPQ.cc
// the mass sets are done as vector here, so I can use size
// eliminates mistyping values for m,n
double Rmnj(const int j,const std::vector<double> Msq,
	    const std::vector<double> Musq);

// derivative w.r.t. m^2_i of Rmnj
double Rimnj(const int i,const int j,const std::vector<double> Msq,
	     const std::vector<double> Musq);

//************************************************************************
// f_+

// staggered partially quenched 22 case
dcomplex fvpv2s2nf3p4RS(const double qsq, const quarkmassnf mass,
		     const double DS, const double DV, const double DA,
		     const double DT, const double dV, const double dA,
		     const bool printmasses){
  // Di are the mass shifts for the different traces (\Delta_i)
  // di are the vertices (\delta_i)
  const double root = 1./4.;// 1/4 if rooted, else 1
  const double DP = 0;
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 4) || (B0mq.size() != 4)){
    std::cout << "ERROR: fvpv2s2nf3p4RS probably called with wrong number"
	 << "of quark masses\n";}
  // 2 valence and 2 sea quark masses
  double m11 = 2.*B0mq[0];
  double m33 = 2.*B0mq[1];
  double m44 = 2.*B0mq[2];
  double m66 = 2.*B0mq[3];

  double mu2 = mu0*mu0;

  double m13,m14,m16,m34,m36,me2,mep2;
  // see notes for mp2,me2
  me2 = (m44+2.*m66)/3.;
  m13 = (m11+m33)/2.;
  m14 = (m11+m44)/2.;
  m16 = (m11+m66)/2.;
  m34 = (m33+m44)/2.;
  m36 = (m33+m66)/2.;
  std::vector<double> DXi = {DP,DS,DV,DA,DT};
  std::vector<double> aXi = {1.,1.,4.,4.,6.};
  // needed for all 5 tastes
  std::vector<double> Mu3[5],Mxf[5],Myf[5];
  for(int taste = 0; taste < 5; taste++){
    Mxf[taste] = {m34+DXi[taste],m36+DXi[taste]};
    Myf[taste] = {m14+DXi[taste],m16+DXi[taste]};
  }
  std::vector<double> nsea = {2.,1.};// first sea mass twice
  // needed for S,V,A
  std::vector<double> DXiD = {DS,DV,DA};
  std::vector<double> aXiD = {1.,4.,4.};
  std::vector<double> dXiD = {-1./(3.*root),dV,dA};
  std::vector<double> M3xx[3],M3yy[3],M4xy[3];
  for(int taste = 0; taste < 3; taste++)
    Mu3[taste] = {m44+DXiD[taste],m66+DXiD[taste]};
  // taste singlet masses
  M3xx[0] = {m33+DXiD[0],me2+DXiD[0]};
  M3yy[0] = {m11+DXiD[0],me2+DXiD[0]};
  M4xy[0] = {m33+DXiD[0],m11+DXiD[0],me2+DXiD[0]};
  // calculate eta, etap mass for V
  double a = 3.*dV*root+m44+m66;
  double b = sqrt(pow(m44-m66,2)+2.*(m44-m66)*dV*root+pow(3.*dV*root,2)); 
  me2 = (a-b)/2.;
  mep2 = (a+b)/2.;
  M3xx[1] = {m33+DXiD[1],me2+DXiD[1],mep2+DXiD[1]};
  M3yy[1] = {m11+DXiD[1],me2+DXiD[1],mep2+DXiD[1]};
  M4xy[1] = {m33+DXiD[1],m11+DXiD[1],me2+DXiD[1],mep2+DXiD[1]};
  // calculate eta, etap mass for A
  a = 3.*dA*root+m44+m66;
  b = sqrt(pow(m44-m66,2)+2.*(m44-m66)*dA*root+pow(3.*dA*root,2)); 
  me2 = (a-b)/2.;
  mep2 = (a+b)/2.;
  M3xx[2] = {m33+DXiD[2],me2+DXiD[2],mep2+DXiD[2]};
  M3yy[2] = {m11+DXiD[2],me2+DXiD[2],mep2+DXiD[2]};
  M4xy[2] = {m33+DXiD[2],m11+DXiD[2],me2+DXiD[2],mep2+DXiD[2]};

  if(printmasses){
    char taste5[5]={'P','S','V','A','T'};
    char taste3[3]={'S','V','A'};
    std::cout << "Mxf"<<std::endl;
    for(int k=0;k<5;k++){
      std::cout <<taste5[k];
      for (int i=0; i < int(Mxf[k].size()); i++)
	std::cout<<' '<< Mxf[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "Myf"<<std::endl;
    for(int k=0;k<5;k++){
      std::cout <<taste5[k];
      for (int i=0; i < int(Myf[k].size()); i++)
	std::cout<<' '<< Myf[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "M3xx"<<std::endl;
    for(int k=0;k<3;k++){
      std::cout <<taste3[k];
      for (int i=0; i < int(M3xx[k].size()); i++)
	std::cout<<' '<< M3xx[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "M3yy"<<std::endl;
    for(int k=0;k<3;k++){
      std::cout <<taste3[k];
      for (int i=0; i < int(M3yy[k].size()); i++)
	std::cout<<' '<< M3yy[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "M4xy"<<std::endl;
    for(int k=0;k<3;k++){
      std::cout <<taste3[k];
      for (int i=0; i < int(M4xy[k].size()); i++)
	std::cout<<' '<< M4xy[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "Mu3"<<std::endl;
    for(int k=0;k<3;k++){
      std::cout <<taste3[k];
      for (int i=0; i < int(Mu3[k].size()); i++)
	std::cout<<' '<< Mu3[k][i];
      std::cout<<std::endl;;
    }
  }


  dcomplex fp432 = 0.;
  dcomplex fp1,fp2,fp3;
  // Note x = strange quark valence = 3
  //      y   up                    = 1
  // The spectator is down quark valence = 2

  // sum over sea quarks
  for(int taste = 0; taste < 5; taste++)
    for(int k=0; k<int(Mxf[taste].size()); k++)
      fp432 += aXi[taste]*root*nsea[k]*
	(Ab(Mxf[taste][k],mu2)+Ab(Myf[taste][k],mu2)
	 -4.*B22b(Mxf[taste][k],Myf[taste][k],qsq,mu2));

  if (printmasses){
    fp1 = fp432/(16.*f0*f0);
    std::cout << "from sea     "<<fp1<<std::endl;}

  for(int taste=0; taste < 3; taste++){
    // Dxx and Dyy terms, derivative on R^[m,n]_k
    for(int k=0; k<int(M3xx[taste].size()); k++){
      fp432 += aXiD[taste]*dXiD[taste]*
	( Rimnj(0,k,M3xx[taste],Mu3[taste])*
	  (Ab(M3xx[taste][k],mu2)
	   -4.*B22b(m13+DXiD[taste],M3xx[taste][k],qsq,mu2))
	  +Rimnj(0,k,M3yy[taste],Mu3[taste])*
	  (Ab(M3yy[taste][k],mu2)
	   -4.*B22b(m13+DXiD[taste],M3yy[taste][k],qsq,mu2))
	  );
    }
    // Dxx and Dyy terms
    // derivative on the integral, needs no sum, only first term contributes
    fp432 += aXiD[taste]*dXiD[taste]*
      (
       +Rmnj(0,M3xx[taste],Mu3[taste])*
       (Bb(M3xx[taste][0],mu2)
	-4.*B22b(3,m13+DXiD[taste],M3xx[taste][0],qsq,mu2))
       +Rmnj(0,M3yy[taste],Mu3[taste])*
       (Bb(M3yy[taste][0],mu2)
	-4.*B22b(3,m13+DXiD[taste],M3yy[taste][0],qsq,mu2))
       );
  }

  if (printmasses){
    fp2 = fp432/(16.*f0*f0)-fp1;
    std::cout << "from Dxx Dyy "<<fp2<<std::endl;}

  // Dxy terms
  for(int taste=0; taste < 3; taste++){
    for(int k=0; k<int(M4xy[taste].size()); k++){
      fp432 += (-2.)*aXiD[taste]*dXiD[taste]*
	Rmnj(k,M4xy[taste],Mu3[taste])*
	(Ab(M4xy[taste][k],mu2)
	 -4.*B22b(m13+DXiD[taste],M4xy[taste][k],qsq,mu2));
    }
  }
  if (printmasses){
    fp3 = fp432/(16.*f0*f0)-fp2-fp1;
    std::cout << "from Dxy     "<<fp3<<std::endl;}
  return fp432/(16.*f0*f0);
}

// staggered partially quenched 23 case
dcomplex fvpv2s3nf3p4RS(const double qsq, const quarkmassnf mass,
		     const double DS, const double DV, const double DA,
		     const double DT, const double dV, const double dA,
		     const bool printmasses){
  // Di are the mass shifts for the different traces (\Delta_i)
  // di are the vertices (\delta_i)
  const double root = 1./4.;// 1/4 if rooted, else 1
  const double DP = 0;

  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 5) || (B0mq.size() != 5)){
    std::cout << "ERROR: fvpv2s3nf3p4RS probably called with wrong number"
	 << "of quark masses\n";}
  // 2 valence and 3 sea quark masses
  double m11 = 2.*B0mq[0];
  double m33 = 2.*B0mq[1];
  double m44 = 2.*B0mq[2];
  double m55 = 2.*B0mq[3];
  double m66 = 2.*B0mq[4];

  double mu2 = mu0*mu0;

  double m13,m14,m15,m16,m34,m35,m36,mp2,me2,mep2;
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
  std::vector<double> DXi = {DP,DS,DV,DA,DT};
  std::vector<double> aXi = {1.,1.,4.,4.,6.};
  // needed for all 5 tastes
  std::vector<double> Mu3[5],Mxf[5],Myf[5];
  for(int taste = 0; taste < 5; taste++){
    Mxf[taste] = {m34+DXi[taste],m35+DXi[taste],m36+DXi[taste]};
    Myf[taste] = {m14+DXi[taste],m15+DXi[taste],m16+DXi[taste]};
  }
  // needed for S,V,A
  std::vector<double> DXiD = {DS,DV,DA};
  std::vector<double> aXiD = {1.,4.,4.};
  std::vector<double> dXiD = {-1./(3.*root),dV,dA};
  std::vector<double> M3xx[3],M3yy[3],M4xy[3];
  for(int taste = 0; taste < 3; taste++)
    Mu3[taste] = {m44+DXiD[taste],m55+DXiD[taste],m66+DXiD[taste]};
  // taste singlet masses
  M3xx[0] = {m33+DXiD[0],mp2+DXiD[0],me2+DXiD[0]};
  M3yy[0] = {m11+DXiD[0],mp2+DXiD[0],me2+DXiD[0]};
  M4xy[0] = {m33+DXiD[0],m11+DXiD[0],mp2+DXiD[0],me2+DXiD[0]};
  // calculate pi, eta, etap mass for V
  double x[3];
  double r,s,t,d;
  r = -3.*dV*root-m44-m55-m66;
  s = 2.*dV*root*(m44+m55+m66)+(m44*m55+m55*m66+m66*m44);
  t = -m44*m55*m66-dV*root*(m44*m55+m55*m66+m66*m44);
  jbdrteq3(r,s,t,x,d);
  if(d > 0.) std::cout << "something wrong with masses in fp4staggered\n"; 
  mp2 = x[0];
  me2 = x[1];
  mep2 = x[2];
  M3xx[1] = {m33+DXiD[1],mp2+DXiD[1],me2+DXiD[1],mep2+DXiD[1]};
  M3yy[1] = {m11+DXiD[1],mp2+DXiD[1],me2+DXiD[1],mep2+DXiD[1]};
  M4xy[1] = {m33+DXiD[1],m11+DXiD[1],mp2+DXiD[1],me2+DXiD[1],mep2+DXiD[1]};
  // calculate pi, eta, etap mass for A
  r = -3.*dA*root-m44-m55-m66;
  s = 2.*dA*root*(m44+m55+m66)+(m44*m55+m55*m66+m66*m44);
  t = -m44*m55*m66-dA*root*(m44*m55+m55*m66+m66*m44);
  jbdrteq3(r,s,t,x,d);
  if(d > 0.) std::cout << "something wrong with masses in fp4staggered\n"; 
  mp2 = x[0];
  me2 = x[1];
  mep2 = x[2];
  M3xx[2] = {m33+DXiD[2],mp2+DXiD[2],me2+DXiD[2],mep2+DXiD[2]};
  M3yy[2] = {m11+DXiD[2],mp2+DXiD[2],me2+DXiD[2],mep2+DXiD[2]};
  M4xy[2] = {m33+DXiD[2],m11+DXiD[2],mp2+DXiD[2],me2+DXiD[2],mep2+DXiD[2]};


  if(printmasses){
    char taste5[5]={'P','S','V','A','T'};
    char taste3[3]={'S','V','A'};
    std::cout << "Mxf"<<std::endl;
    for(int k=0;k<5;k++){
      std::cout <<taste5[k];
      for (int i=0; i < int(Mxf[k].size()); i++)
	std::cout<<' '<< Mxf[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "Myf"<<std::endl;
    for(int k=0;k<5;k++){
      std::cout <<taste5[k];
      for (int i=0; i < int(Myf[k].size()); i++)
	std::cout<<' '<< Myf[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "M3xx"<<std::endl;
    for(int k=0;k<3;k++){
      std::cout <<taste3[k];
      for (int i=0; i < int(M3xx[k].size()); i++)
	std::cout<<' '<< M3xx[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "M3yy"<<std::endl;
    for(int k=0;k<3;k++){
      std::cout <<taste3[k];
      for (int i=0; i < int(M3yy[k].size()); i++)
	std::cout<<' '<< M3yy[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "M4xy"<<std::endl;
    for(int k=0;k<3;k++){
      std::cout <<taste3[k];
      for (int i=0; i < int(M4xy[k].size()); i++)
	std::cout<<' '<< M4xy[k][i];
      std::cout<<std::endl;;
    }
    std::cout << "Mu3"<<std::endl;
    for(int k=0;k<3;k++){
      std::cout <<taste3[k];
      for (int i=0; i < int(Mu3[k].size()); i++)
	std::cout<<' '<< Mu3[k][i];
      std::cout<<std::endl;;
    }
  }


  dcomplex fp433 = 0.;
  dcomplex fp1,fp2,fp3;

  // Note x = strange quark valence = 3
  //      y   up                    = 1
  // The spectator is down quark valence = 2

  // sum over sea quarks
  for(int taste = 0; taste < 5; taste++)
    for(int k=0; k<int(Mxf[taste].size()); k++)
      fp433 += aXi[taste]*root*(Ab(Mxf[taste][k],mu2)+Ab(Myf[taste][k],mu2)
				-4.*B22b(Mxf[taste][k],Myf[taste][k],qsq,mu2));

  if (printmasses){
    fp1 = fp433/(16.*f0*f0);
    std::cout << "from sea     "<<fp1<<std::endl;}


  for(int taste=0; taste < 3; taste++){
    // Dxx and Dyy terms, derivative on R^[m,n]_k
    for(int k=0; k<int(M3xx[taste].size()); k++){
      fp433 += aXiD[taste]*dXiD[taste]*
	( Rimnj(0,k,M3xx[taste],Mu3[taste])*
	  (Ab(M3xx[taste][k],mu2)
	   -4.*B22b(m13+DXiD[taste],M3xx[taste][k],qsq,mu2))
	  +Rimnj(0,k,M3yy[taste],Mu3[taste])*
	  (Ab(M3yy[taste][k],mu2)
	   -4.*B22b(m13+DXiD[taste],M3yy[taste][k],qsq,mu2))
	  );
    }
    // Dxx and Dyy terms
    // derivative on the integral, needs no sum, only first term contributes
    fp433 += aXiD[taste]*dXiD[taste]*
      (
       +Rmnj(0,M3xx[taste],Mu3[taste])*
       (Bb(M3xx[taste][0],mu2)
	-4.*B22b(3,m13+DXiD[taste],M3xx[taste][0],qsq,mu2))
       +Rmnj(0,M3yy[taste],Mu3[taste])*
       (Bb(M3yy[taste][0],mu2)
	-4.*B22b(3,m13+DXiD[taste],M3yy[taste][0],qsq,mu2))
       );
  }
  if (printmasses){
    fp2 = fp433/(16.*f0*f0)-fp1;
    std::cout << "from Dxx Dyy "<<fp2<<std::endl;}

  // Dxy terms
  for(int taste=0; taste < 3; taste++){
    for(int k=0; k<int(M4xy[taste].size()); k++)
      fp433 += (-2.)*aXiD[taste]*dXiD[taste]*
	Rmnj(k,M4xy[taste],Mu3[taste])*
	(Ab(M4xy[taste][k],mu2)
	 -4.*B22b(m13+DXiD[taste],M4xy[taste][k],qsq,mu2));
  }

  if (printmasses){
    fp3 = fp433/(16.*f0*f0)-fp2-fp1;
    std::cout << "from Dxy     "<<fp3<<std::endl;}
  return fp433/(16.*f0*f0);
}

// versions with three valence masses
// call the versions with two masses since result does not depend
// on the spectator mass

dcomplex fvpv3s2nf3p4RS(const double qsq, const quarkmassnf mass,
			const double DS, const double DV, const double DA,
			const double DT, const double dV, const double dA,
			const bool printmasses){
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 5) || (B0mq.size() != 5)){
    std::cout << "ERROR: fvpv3s2nf3p4RS probably called with wrong number"
	 << "of quark masses\n";}
  quarkmassnf mass2({B0mq[0],B0mq[2],B0mq[3],B0mq[4]},f0,mu0);
  return fvpv2s2nf3p4RS(qsq, mass2, DS, DV, DA, DT, dV, dA, printmasses);
}



dcomplex fvpv3s3nf3p4RS(const double qsq, const quarkmassnf mass,
			const double DS, const double DV, const double DA,
			const double DT, const double dV, const double dA,
			const bool printmasses){
  vector<double> B0mq;
  double f0,mu0;
  int nq;
  mass.out(B0mq,f0,mu0,nq);
  if ( (nq != 6) || (B0mq.size() != 6)){
    std::cout << "ERROR: fvpv3s2nf3p4RS probably called with wrong number"
	 << "of quark masses\n";}
  quarkmassnf mass2({B0mq[0],B0mq[2],B0mq[3],B0mq[4],B0mq[5]},f0,mu0);
  return fvpv2s3nf3p4RS(qsq, mass2, DS, DV, DA, DT, dV, dA, printmasses);
}
