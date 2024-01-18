// li.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2017 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// implementation of the li and libar classes

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<cmath>

#include "li.h"

using namespace std;

const double pi16 = 1./(16.*M_PI*M_PI);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

li::li(const double li1r,const double li2r,const double li3r,
       const double li4r,const double li5r,const double li6r,
       const double li7r,const double hi1r,const double hi2r,
       const double hi3r,
       const double muin,
       const std::string Name){
    l1r = li1r;
    l2r = li2r; 
    l3r = li3r;
    l4r = li4r;
    l5r = li5r;
    l6r = li6r;
    l7r = li7r;
    h1r = hi1r;
    h2r = hi2r; 
    h3r = hi3r;
    mu = muin;
    name = Name;
}

li::li(const libar libarin, const double muin, const double mpi){
  const double pim32 = 1./(32.*M_PI*M_PI);
  const double g1 = 1./3.;
  const double g2 = 2./3.;
  const double g3 = -1./2.;
  const double g4 = 2.;
  const double g5 = -1./6.;
  const double g6 = -1./3.;
  const double g7 = 0.;
  const double g8 = 2.;
  const double g9 = 1./12.;
  const double g10= 0.;
  double l1s,l2s,l3s,l4s,l5s,l6s,l7s,h1s,h2s,h3s;
  string nameout;
  libarin.out(l1s,l2s,l3s,l4s,l5s,l6s,l7s,h1s,h2s,h3s,nameout);
  double lp = 2.*log(mpi/muin);

  l1r = g1*pim32*(l1s+lp);
  l2r = g2*pim32*(l2s+lp);
  l3r = g3*pim32*(l3s+lp);
  l4r = g4*pim32*(l4s+lp);
  l5r = g5*pim32*(l5s+lp);
  l6r = g6*pim32*(l6s+lp);
  l7r = l7s;
  h1r = g8*pim32*(h1s+lp);
  h2r = g9*pim32*(h2s+lp);
  h3r = h3s;
  mu = muin;
  name = nameout;
}

li::~li(void){}

std::ostream & operator<<(std::ostream & os, const li & bb){
  os   << std::setprecision(15)
       << "# li set : " << bb.name << "\n"
       << "#mu : "<<bb.mu <<"\n"
       << "#l1r  l2r  l3r :" <<bb.l1r <<' '<<bb.l2r<<' '<<bb.l3r<<"\n"
       << "#l4r  l5r  l6r :" <<bb.l4r <<' '<<bb.l5r<<' '<<bb.l6r<<"\n"
       << "#l7r  l8r  l9r :" <<bb.l7r <<' '<<bb.h1r<<' '<<bb.h2r<<"\n"
       << "#h3r           :" <<bb.h3r<<"\n";
  os.unsetf(std::ios_base::floatfield);
  os << std::setprecision(6);
  return os;
}

std::istream & operator>>(std::istream & is, li & liout){
  std::string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  std::string temps2;
  double l1r,l2r,l3r,l4r,l5r,l6r,l7r,h1r,h2r,h3r,mu;
  std::string name;
  std::getline(is,temps2,':');// reads in characters up to ':'
  std::getline(is,temps); // to remove end of line
  name = temps;
  if (temps2 != "# li set ") std::cout << "trouble reading in li\n";
  std::getline(is,temps2,':');   is >> mu;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l1r >> l2r >> l3r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l4r >> l5r >> l6r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l7r >> h1r >> h2r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> h3r;  getline(is,temps);
  liout = li(l1r,l2r,l3r,l4r,l5r,l6r,l7r,h1r,h2r,h3r,mu,name);
  return is;
}

void li::setli(const int n, const double lin){
  switch(n){
  case 1: l1r = lin; break;
  case 2: l2r = lin; break;  
  case 3: l3r = lin; break;  
  case 4: l4r = lin; break;  
  case 5: l5r = lin; break;  
  case 6: l6r = lin; break;  
  case 7: l7r = lin; break;  
  case 8: h1r = lin; break;  
  case 9: h2r = lin; break;  
  case 10: h3r = lin; break;
  default: std::cout << "Trying to set li number "<<n<<" in set "<<name<<'\n';
  }
}

void li::setli(const double lin,const int n){
  setli(n,lin);
}

void li::setmu(const double muin){
  mu = muin;
}
void li::setname(const std::string inputname){
  name = inputname;
}
// summing two sets of li
li li::operator+(const li & bb) const{
  double l1s,l2s,l3s,l4s,l5s,l6s,l7s,h1s,h2s,h3s,muout;
  l1s = l1r+bb.l1r;
  l2s = l2r+bb.l2r;
  l3s = l3r+bb.l3r;
  l4s = l4r+bb.l4r;
  l5s = l5r+bb.l5r;
  l6s = l6r+bb.l6r;
  l7s = l7r+bb.l7r;
  h1s = h1r+bb.h1r;
  h2s = h2r+bb.h2r;
  h3s = h3r+bb.h3r;
  muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
    std::cout << "WARNING: summing li \""<< name << "\" and \"" << bb.name 
	   << "\" with different mu " << "\n";}
  return li(l1s,l2s,l3s,l4s,l5s,l6s,l7s,h1s,h2s,h3s,muout);
}
// difference of two sets of li
li li::operator-(const li & bb) const{
  double l1s,l2s,l3s,l4s,l5s,l6s,l7s,h1s,h2s,h3s,muout;
  l1s = l1r-bb.l1r;
  l2s = l2r-bb.l2r;
  l3s = l3r-bb.l3r;
  l4s = l4r-bb.l4r;
  l5s = l5r-bb.l5r;
  l6s = l6r-bb.l6r;
  l7s = l7r-bb.l7r;
  h1s = h1r-bb.h1r;
  h2s = h2r-bb.h2r;
  h3s = h3r-bb.h3r;
  muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
    std::cout << "WARNING: summing li \""<< name << "\" and \"" << bb.name 
	   << "\" with different mu " << "\n";}
  return li(l1s,l2s,l3s,l4s,l5s,l6s,l7s,h1s,h2s,h3s,muout);
}

// multiplying an li with a number
li li::operator*(const double & aa) const{
  return li(aa*l1r,aa*l2r,aa*l3r,aa*l4r,aa*l5r,aa*l6r,aa*l7r,aa*h1r,
	    aa*h2r,aa*h3r,mu);
}

li operator*(const double & aa,const li & bb){
  return li(aa*bb.l1r,aa*bb.l2r,aa*bb.l3r,aa*bb.l4r,aa*bb.l5r,aa*bb.l6r,
	    aa*bb.l7r,aa*bb.h1r,aa*bb.h2r,aa*bb.h3r,bb.mu);
}

// changing scale
void li::changescale(const double newmu){
  const double g1 = 1./3.;
  const double g2 = 2./3.;
  const double g3 = -1./2.;
  const double g4 = 2.;
  const double g5 = -1./6.;
  const double g6 = -1./3.;
  const double g7 = 0.;
  const double g8 = 2.;
  const double g9 = 1./12.;
  const double g10= 0.;

  double logmm = log(mu/newmu);

  l1r =  l1r+pi16*g1*logmm;
  l2r =  l2r+pi16*g2*logmm;
  l3r =  l3r+pi16*g3*logmm;
  l4r =  l4r+pi16*g4*logmm;
  l5r =  l5r+pi16*g5*logmm;
  l6r =  l6r+pi16*g6*logmm;
  l7r =  l7r+pi16*g7*logmm;
  h1r =  h1r+pi16*g8*logmm;
  h2r =  h2r+pi16*g9*logmm;
  h3r =  h3r+pi16*g10*logmm;
  mu = newmu;
}

// output functions
void li::out(double & l1t,double & l2t,double & l3t,
	     double & l4t,double & l5t,double & l6t,
	     double & l7t,
	     double & h1t,double & h2t,double & h3t,
	     double & mut, std::string nameout) const{
  l1t = l1r;
  l2t = l2r;
  l3t = l3r;
  l4t = l4r;
  l5t = l5r;
  l6t = l6r;
  l7t = l7r;
  h1t = h1r;
  h2t = h2r;
  h3t = h3r;
  mut =mu;
  nameout = name;
}

void li::out(double & l1t,double & l2t,double & l3t,
	     double & l4t,double & l5t,double & l6t,
	     double & l7t,
	     double & h1t,double & h2t,double & h3t,
	     double & mut) const{
  l1t = l1r;
  l2t = l2r;
  l3t = l3r;
  l4t = l4r;
  l5t = l5r;
  l6t = l6r;
  l7t = l7r;
  h1t = h1r;
  h2t = h2r;
  h3t = h3r;
  mut =mu;
}
void li::out(double & l1t,double & l2t,double & l3t,
	     double & l4t,double & l5t,double & l6t,
	     double & l7t,
	     double & h1t,double & h2t,double & h3t) const{
  l1t = l1r;
  l2t = l2r;
  l3t = l3r;
  l4t = l4r;
  l5t = l5r;
  l6t = l6r;
  l7t = l7r;
  h1t = h1r;
  h2t = h2r;
  h3t = h3r;
}

void li::out(double & l1t,double & l2t,double & l3t,
	     double & l4t,double & l5t,double & l6t,
	     double & l7t,
	     double & mut) const{
  l1t = l1r;
  l2t = l2r;
  l3t = l3r;
  l4t = l4r;
  l5t = l5r;
  l6t = l6r;
  l7t = l7r;
  mut =mu;
}
void li::out(double & l1t,double & l2t,double & l3t,
	     double & l4t,double & l5t,double & l6t,
	     double & l7t) const{
  l1t = l1r;
  l2t = l2r;
  l3t = l3r;
  l4t = l4r;
  l5t = l5r;
  l6t = l6r;
  l7t = l7r;
}

double li::out(const int n) const{
  switch(n){
  case 1:  return l1r ; break;
  case 2:  return l2r ; break;
  case 3:  return l3r ; break;
  case 4:  return l4r ; break;
  case 5:  return l5r ; break;
  case 6:  return l6r ; break;
  case 7:  return l7r ; break;
  case 8:  return h1r ; break;
  case 9:  return h2r ; break;
  case 10: return h3r ; break;
  default:
      std::cout << "Trying to output li number "<<n<<" in set "<<name<<'\n';
  }
  return 0.;
}

double li::getmu(void) const{
  return mu;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

libar::libar(const double l1barin,const double l2barin,const double l3barin,
	     const double l4barin,const double l5barin,const double l6barin,
	     const double l7rin,const double h1barin,const double h2barin,
	     const double h3rin,
	     const std::string Namein){
  l1bar = l1barin;
  l2bar = l2barin;
  l3bar = l3barin;
  l4bar = l4barin;
  l5bar = l5barin;
  l6bar = l6barin;
  l7r   = l7rin; 
  h1bar = h1barin;
  h2bar = h2barin;
  h3r   = h3rin;
  name = Namein;
}

libar::libar(const li liin, const double mpi){
  const double pi32 = 32.*M_PI*M_PI;
  const double g1 = 1./3.;
  const double g2 = 2./3.;
  const double g3 = -1./2.;
  const double g4 = 2.;
  const double g5 = -1./6.;
  const double g6 = -1./3.;
  const double g7 = 0.;
  const double g8 = 2.;
  const double g9 = 1./12.;
  const double g10= 0.;
  double l1s,l2s,l3s,l4s,l5s,l6s,l7s,h1s,h2s,h3s,mu;
  liin.out(l1s,l2s,l3s,l4s,l5s,l6s,l7s,h1s,h2s,h3s,mu);
  double lp = 2.*log(mpi/mu);

  l1bar = -lp+pi32/g1*l1s;
  l2bar = -lp+pi32/g2*l2s;
  l3bar = -lp+pi32/g3*l3s;
  l4bar = -lp+pi32/g4*l4s;
  l5bar = -lp+pi32/g5*l5s;
  l6bar = -lp+pi32/g6*l6s;
  l7r = l7s;
  h1bar = -lp+pi32/g8*h1s;
  h2bar = -lp+pi32/g9*h2s;
  h3r = h3s;
}

libar::~libar(void){}

ostream & operator<<(ostream & os, const libar & bb){
  os  << setprecision(15)
      << "# libar set : " << bb.name << "\n"
      << "#l1bar l2bar l3bar :"<<bb.l1bar <<' '<<bb.l2bar<<' '<<bb.l3bar<<"\n"
      << "#l4bar l5bar l6bar :"<<bb.l4bar <<' '<<bb.l5bar<<' '<<bb.l6bar<<"\n"
      << "#l7r   h1bar h2bar :" <<bb.l7r <<' '<<bb.h1bar<<' '<<bb.h2bar<<"\n"
      << "#h3r               :" <<bb.h3r<<"\n";
  return os;
}

istream & operator>>(istream & is, libar & Liout){
  string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  string temps2;
  double l1bar,l2bar,l3bar,l4bar,l5bar,l6bar,l7r;
  double h1bar,h2bar,h3r;
  string name;
  getline(is,temps2,':');// reads in characters up to ':'
  getline(is,temps); // to remove end of line
  name = temps;
  if (temps2 != "# libar set ") cout << "trouble reading in libar\n";
  getline(is,temps2,':');   is >> l1bar >> l2bar >> l3bar;  getline(is,temps);
  getline(is,temps2,':');   is >> l4bar >> l5bar >> l6bar;  getline(is,temps);
  getline(is,temps2,':');   is >> l7r >> h1bar >> h2bar;  getline(is,temps);
  getline(is,temps2,':');   is >> h3r;  getline(is,temps);
  Liout = libar(l1bar,l2bar,l3bar,l4bar,l5bar,l6bar,l7r,h1bar,h2bar,h3r,name);
  return is;
}

void libar::setlibar(const int n, const double lbarin){
  switch(n){
  case 1: l1bar = lbarin; break;
  case 2: l2bar = lbarin; break;  
  case 3: l3bar = lbarin; break;  
  case 4: l4bar = lbarin; break;  
  case 5: l5bar = lbarin; break;  
  case 6: l6bar = lbarin; break;  
  case 7: l7r = lbarin; break;  
  case 8: h1bar = lbarin; break;  
  case 9: h2bar = lbarin; break;  
  case 10: h3r = lbarin; break;
  default: cout << "Trying to set Li number "<<n<<" in set "<<name<<endl;
  }
}

void libar::setlibar(const double lin,const int n){
  setlibar(n,lin);
}

void libar::setname(const string inputname){
  name = inputname;
}

void libar::out(double & l1t,double & l2t,double & l3t,
	     double & l4t,double & l5t,double & l6t,
	     double & l7t,double & h1t,double & h2t,
	     double & h3t, string nameout) const{
  l1t = l1bar;
  l2t = l2bar;
  l3t = l3bar;
  l4t = l4bar;
  l5t = l5bar;
  l6t = l6bar;
  l7t = l7r;
  h1t = h1bar;
  h2t = h2bar;
  h3t = h3r;
  nameout = name;
}

void libar::out(double & l1t,double & l2t,double & l3t,
		     double & l4t,double & l5t,double & l6t,
		     double & l7t) const{
  l1t = l1bar;
  l2t = l2bar;
  l3t = l3bar;
  l4t = l4bar;
  l5t = l5bar;
  l6t = l6bar;
  l7t = l7r;
}

double libar::out(const int n) const{
  switch(n){
  case 1: return l1bar; break;
  case 2: return l2bar; break;  
  case 3: return l3bar; break;  
  case 4: return l4bar; break;  
  case 5: return l5bar; break;  
  case 6: return l6bar; break;  
  case 7: return l7r  ; break;  
  case 8: return h1bar; break;  
  case 9: return h2bar; break;  
  case 10: return h3r ; break;
  default: cout << "Trying to output libar number "<<n<<" in set "<<name<<endl;
    return 0.;
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// from three to two flavours
//
//const double pi16 = 1./(16.*M_PI*M_PI);
//
//libar libarfromLip4(const physmass massin, const Li Liin){
//  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
//  double l1r,l2r,l3r,l4r,l5r,l6r;//,l7r,h1r,h2r,h3r;
//  Liin.Liout(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
//  double mk2 = pow(massin.mk,2);
//  double me2 = pow(massin.meta,2);
//  double mp2 = pow(massin.mpi,2);
//  double mu2 = pow(mu,2);
//  double mk2bar = mk2-1./2.*mp2;
//  double me2bar = me2-1./3.*mp2;
//  const double g1 = 1./3.;
//  const double g2 = 2./3.;
//  const double g3 = -1./2.;
//  const double g4 = 2.;
//  const double g5 = -1./6.;
//  const double g6 = -1./3.;
//  //const double g7 = 0.;
//  //const double g8 = 2.;
//  //const double g9 = 1./12.;
//  //const double g10= 0.;
//  double le = log(me2bar/mu2);
//  double lk = log(mk2bar/mu2);
//  double nue = pi16/2.*(le+1.);
//  double nuk = pi16/2.*(lk+1);
//
//  l1r = -nuk/24.+4.*L1r+2.*L3r;
//  l2r = -nuk/12.+4.*L2r;
//  l3r = -nue/18.-8.*L4r-4.*L5r+16.*L6r+8.*L8r;
//  l4r = -nuk/2.+8.*L4r+4.*L5r;
//  l5r = nuk/12.+L10r;
//  l6r = nuk/6.-2.*L9r;
//  //l7r = 0.;
//  //l7r, h1r,h2r,h3r NOT implemented
//
//  // turning them in barred versions
//  double lp = log(mp2/mu2);
//  const double pi32 = 32.*M_PI*M_PI;
//
//  l1r = -lp+pi32/g1*l1r;
//  l2r = -lp+pi32/g2*l2r;
//  l3r = -lp+pi32/g3*l3r;
//  l4r = -lp+pi32/g4*l4r;
//  l5r = -lp+pi32/g5*l5r;
//  l6r = -lp+pi32/g6*l6r;
//  return libar(l1r,l2r,l3r,l4r,l5r,l6r);
//}
//
//double FoverF0p4(const physmass massin, const Li Liin){
//  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
//  Liin.Liout(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
//  double mk2 = pow(massin.mk,2);
//  double me2 = pow(massin.meta,2);
//  double mp2 = pow(massin.mpi,2);
//  double mu2 = pow(mu,2);
//
//  //double F0 = massin.fpi/(1.+fpi4(massin,Liin));// only needed to NLO
//          
//  //double le = log(me2bar2/mu2);
//  double lk = log(mk2/mu2);
//  //double nue = pi16/2.*(le+1.);
//  //double nuk = pi16/2.*(lk+1.);
//
//
//
//  double aF = -1./2.*lk+8.*L4r/pi16;
//  double x = pi16*mk2/pow(massin.fpi,2);
//  double FoverF0 = 1.+aF*x;
//
//  return FoverF0;
//}
//
//double FoverF0p6(const physmass massin, const Li Liin, const Ci Ciin){
//  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
//  Liin.Liout(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
//  double Cir[95],mu6;
//  Ciin.Ciout(Cir,mu6);
//  double mk2 = pow(massin.mk,2);
//  double me2 = pow(massin.meta,2);
//  double mp2 = pow(massin.mpi,2);
//  double mu2 = pow(mu,2);
//  double mk2bar1 = mk2-1./2.*mp2;
//  double me2bar1 = me2-1./3.*mp2;
//  physmass tempmass(1e-6,sqrt(mk2bar1),sqrt(me2bar1),massin.fpi,massin.mu);
//  double mp2temp = mp2-mpi4(massin,Liin);
//  double mk2bar2 = mk2+mk4(tempmass,Liin)-mk4(massin,Liin)-1./2.*mp2temp;
//  //double me2bar2 = me2+meta4(tempmass,Liin)-meta4(massin,Liin)-1./3.*mp2temp;
//
//  double F0 = massin.fpi/(1.+fpi4(massin,Liin));// only needed to NLO
//          
//  //double le = log(me2bar2/mu2);
//  double lk = log(mk2bar2/mu2);
//  //double nue = pi16/2.*(le+1.);
//  //double nuk = pi16/2.*(lk+1.);
//
//
//
//  double aF = -1./2.*lk+8.*L4r/pi16;
//
//  // F and Sigma still fully present below in libarfromLip6
//
//  // notation pij = pi from the paper for lbarj
//  const double rho1 = 1.41602;// some constant from Clausius Cl_2 function
//  //const double rho2 = 1.75793;//  same
//  const double ln43 = log(4./3.);
//  const double NN = 16.*M_PI*M_PI;
//  double p0F =
//    - 73./32
//    + 2./3.*rho1
//    + 1./3.*ln43 
//    + NN * (
//	    - 52./9.*L2r
//	    - 43./27.*L3r
//	    )
//    + pow(NN,2) * (
//		   96.*pow(L4r,2)
//		   + 64.*L4r*L5r
//		   - 256.*L4r*L6r
//		   - 128.*L4r*L8r
//		   + 32.*Cir[16]
//		   )
//    + ln43*NN * (
//		 128./9.*L1r
//		 + 32./9.*L2r
//		 + 32./9.*L3r
//		 - 32./3.*L4r
//		 );
//  double p1F =
//    7./3.
//    - 1./3.*ln43 
//    + NN * (
//            416./9.*L1r
//	    + 104./9.*L2r
//	    + 122./9.*L3r
//	    - 68./3.*L4r
//	    )
//    ;
//
//  const double p2F = - 11./12.;
//  double x = pi16*mk2bar2/pow(F0,2);
//  double FoverF0 = 1.+aF*x+x*x*(p0F+p1F*lk+p2F*lk*lk);
//
//  return FoverF0;
//}
//
//
//libar libarfromLip6(const physmass massin, const Li Liin, const Ci Ciin){
//  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
//  Liin.Liout(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
//  //F, Sigma, h1r,h2r,h3r NOT implemented
//  // l7r does not contain the F_0^2/(8 B0 m_s) term
//  // the constants needed for F and Sigma are present but not used
//  double Cir[95],mu6;
//  Ciin.Ciout(Cir,mu6);
//  if (mu6 != mu){
//    cout << "different mu's in libarfromLip6\n";
//  }
//  double mk2 = pow(massin.mk,2);
//  double me2 = pow(massin.meta,2);
//  double mp2 = pow(massin.mpi,2);
//  double mu2 = pow(mu,2);
//  double mk2bar1 = mk2-1./2.*mp2;
//  double me2bar1 = me2-1./3.*mp2;
//  physmass tempmass(1e-6,sqrt(mk2bar1),sqrt(me2bar1),massin.fpi,massin.mu);
//  double mp2temp = mp2-mpi4(massin,Liin);
//  double mk2bar2 = mk2+mk4(tempmass,Liin)-mk4(massin,Liin)-1./2.*mp2temp;
//  double me2bar2 = me2+meta4(tempmass,Liin)-meta4(massin,Liin)-1./3.*mp2temp;
//          
//  const double g1 = 1./3.;
//  const double g2 = 2./3.;
//  const double g3 = -1./2.;
//  const double g4 = 2.;
//  const double g5 = -1./6.;
//  const double g6 = -1./3.;
//  //const double g7 = 0.;
//  //const double g8 = 2.;
//  //const double g9 = 1./12.;
//  //const double g10= 0.;
//  double le = log(me2bar2/mu2);
//  double lk = log(mk2bar2/mu2);
//  double nue = pi16/2.*(le+1.);
//  double nuk = pi16/2.*(lk+1.);
//
//  //double aF = -1./2.*lk+8.*L4r/pi16;
//  //double asigma = -lk-2./9.*le+32.*L6r/pi16;
//  double a1r = -nuk/24.+4.*L1r+2.*L3r;
//  double a2r = -nuk/12.+4.*L2r;
//  double a3r = -nue/18.-8.*L4r-4.*L5r+16.*L6r+8.*L8r;
//  double a4r = -nuk/2.+8.*L4r+4.*L5r;
//  double a5r = nuk/12.+L10r;
//  double a6r = nuk/6.-2.*L9r;
//  double a7r = -5./18.*pi16+nuk/2.+5./9.*nue+4.*L4r-4.*L6r-36.*L7r-12.*L8r;
//
//  // notation pij = pi from the paper for lbarj
//  const double rho1 = 1.41602;// some constant from Clausius Cl_2 function
//  const double rho2 = 1.75793;//  same
//  const double ln43 = log(4./3.);
//  const double NN = 16.*M_PI*M_PI;
//  //double p0F =
//  //  - 73./32
//  //  + 2./3.*rho1
//  //  + 1./3.*ln43 
//  //  + NN * (
//  //	    - 52./9.*L2r
//  //	    - 43./27.*L3r
//  //	    )
//  //  + pow(NN,2) * (
//  //		   96.*pow(L4r,2)
//  //		   + 64.*L4r*L5r
//  //		   - 256.*L4r*L6r
//  //		   - 128.*L4r*L8r
//  //		   + 32.*Cir[16]
//  //		   )
//  //  + ln43*NN * (
//  //		 128./9.*L1r
//  //		 + 32./9.*L2r
//  //		 + 32./9.*L3r
//  //		 - 32./3.*L4r
//  //		 );
//  //double p0sigma =
//  //  26./81.*pow(ln43,2)
//  //  + pow(NN,2) * (
//  //		   512.*L4r*L6r
//  //		   + 256.*L5r*L6r
//  //		   - 1024.*pow(L6r,2)
//  //		   - 512.*L6r*L8r
//  //		   + 64.*Cir[20]
//  //		   + 192.*Cir[21]
//  //		   )
//  //  + ln43*NN * (
//  //		 160./9.*L4r
//  //		 + 16./3.*L5r
//  //		 - 448./9.*L6r
//  //		 + 64./3.*L7r
//  //		 - 32./9.*L8r
//  //		 );
//  double p01 =
//    pi16 * (
//	    - 73./288.
//	    + 1./24.*ln43
//	    + 1./16.*rho1
//	    )
//    + 8.*L1r
//    + 2.*L3r
//    - 4.*L4r
//    + NN * (
//            8.*Cir[6]
//	    - 8.*Cir[11]
//	    + 32.*Cir[13]
//	    );
//  double p02 =
//    pi16 * (
//            433./288.
//	    - 1./24.*ln43
//	    + 1./16.*rho1
//	    )
//    + NN * (
//            16.*Cir[11]
//	    - 32.*Cir[13]
//	    );
//  double p03 =
//    pi16 * (
//            1075./2592.
//	    - 79./288.*rho1
//	    )
//    - 383./1296.*ln43*pi16
//    + 5./81.*pow(ln43,2)*pi16 
//    + NN * (
//	    - 128.*pow(L4r,2)
//	    - 64.*L4r*L5r
//	    + 768.*L4r*L6r
//	    + 256.*L4r*L8r
//	    + 128.*L5r*L6r
//	    - 1024.*pow(L6r,2)
//	    - 512.*L6r*L8r
//	    - 32.*Cir[13]
//	    - 16.*Cir[15]
//	    + 32.*Cir[20]
//	    + 192.*Cir[21]
//	    + 32.*Cir[32]
//	    )
//    + ln43 * (
//	      - 64./9.*L1r
//	      - 16./9.*L2r
//	      - 16./9.*L3r
//	      + 80./9.*L4r
//	      + 16./9.*L5r
//	      - 64./9.*L6r
//	      - 32./9.*L8r
//	      )
//    - 176./9.*L1r
//    - 124./27.*L3r
//    + 64./3.*L4r
//    + 140./27.*L5r
//    - 208./9.*L6r
//    + 64./9.*L7r
//    - 8.*L8r
//    ;
//  double p04 =
//    pi16 * (
//            67./144.
//	    + 11./24.*rho1
//	    )
//    - 5./36.*ln43*pi16 
//    + NN * (
//            128.*pow(L4r,2)
//	    + 64.*L4r*L5r
//	    - 256.*L4r*L6r
//	    - 128.*L5r*L6r
//	    + 16.*Cir[15]
//	    )
//    + ln43 * (
//	      64./9.*L1r
//	      + 16./9.*L2r
//	      + 16./9.*L3r
//	      - 16./9.*L4r
//	      )
//    + 176./9.*L1r
//    + 124./27.*L3r
//    - 16./9.*L4r
//    + 4.*L5r
//    - 16.*L6r
//    - 8.*L8r
//    ;
//  double p05 =
//    pi16 * (
//	    - 67./576.
//	    + 7./64.*rho1
//	    )
//    + 5./96.*ln43*pi16 
//    + NN * (
//	    - 8.*Cir[13]
//	    + 8.*Cir[62]
//	    - 8.*Cir[81]
//	    )
//    ;
//  double p06 =
//    pi16 * (
//	    - 163./288.
//	    - 1./16.*rho1
//	    )
//    + 1./24.*ln43*pi16 
//    + NN * (
//            32.*Cir[13]
//	    + 8.*Cir[64]
//	    )
//    ;
//  double p07 =
//    pi16 * (
//	    - 1937./576.
//	    + 5./9.*rho1
//	    + 2./27.*rho2
//	    )
//    - 1./288.
//    + 25./18.*ln43*pi16
//    - 22./81.*pow(ln43,2)*pi16
//    + NN * (
//            1152.*pow(L7r,2)
//	    - 1152.*L7r*L4r
//	    + 2304.*L7r*L6r
//	    + 768.*L7r*L8r
//	    + 32.*pow(L4r,2)
//	    - 128.*L4r*L6r
//	    - 384.*L4r*L8r
//	    + 128.*pow(L6r,2)
//	    + 768.*L6r*L8r
//	    + 128.*pow(L8r,2)
//	    + 16.*Cir[16]
//	    - 24.*Cir[19]
//	    - 56.*Cir[20]
//	    - 24.*Cir[21]
//	    - 16.*Cir[31]
//	    - 48.*Cir[32]
//	    - 48.*Cir[33]
//	    )
//    + ln43 * (
//	      64./9.*L1r
//	      + 16./9.*L2r
//	      + 16./9.*L3r
//	      - 16./3.*L4r
//	      - 80./9.*L5r
//	      + 32./9.*L6r
//	      - 32./3.*L7r
//	      + 128./9.*L8r
//	      )
//    - 26./9.*L2r
//    - 43./54.*L3r
//    - 8.*L5r
//    + 16.*L8r
//    ;
//  //double p1F =
//  //  7./3.
//  //  - 1./3.*ln43 
//  //  + NN * (
//  //          416./9.*L1r
//  //	    + 104./9.*L2r
//  //	    + 122./9.*L3r
//  //	    - 68./3.*L4r
//  //	    )
//  //  ;
//  //double p1sigma =
//  //  - 2./81.*ln43 
//  //  + NN * (
//  //          592./9.*L4r
//  //	    + 64./3.*L5r
//  //	    - 1312./9.*L6r
//  //	    + 64./3.*L7r
//  //	    - 320./9.*L8r
//  //	    )
//  //  ;
//  double p11 = 
//    5./6.*pi16 
//    + 8.*L1r
//    + 4.*L2r
//    + 3.*L3r
//    - 4.*L4r
//    ;
//  double p12 =
//    13./24.*pi16 
//    - 8.*L2r
//    - 2.*L3r
//    ;
//  double p13 =
//    - 1501./648.*pi16 
//    + 11./162.*ln43*pi16 
//    - 352./9.*L1r
//    - 88./9.*L2r
//    - 106./9.*L3r
//    + 440./9.*L4r
//    + 124./9.*L5r
//    - 496./9.*L6r
//    - 248./9.*L8r
//    ;
//  double p14 =
//    37./18.*pi16 
//    - 7./18.*ln43*pi16 
//    + 352./9.*L1r
//    + 88./9.*L2r
//    + 106./9.*L3r
//    - 88./9.*L4r
//    - 16.*L6r
//    - 8.*L8r
//    ;
//  double p15 =
//    17./48.*pi16 
//    - L9r
//    - 2.*L10r
//    ;
//  double p16 =
//    - 7./24.*pi16
//    + 2.*L3r
//    + 2.*L9r
//    ;
//  double p17 =
//    73./18.*pi16
//    + 101./162.*ln43*pi16
//    + 208./9.*L1r
//    + 52./9.*L2r
//    + 61./9.*L3r
//    - 52./3.*L4r
//    - 224./9.*L5r
//    + 104./9.*L6r
//    + 184./3.*L7r
//    + 632./9.*L8r
//    ;
//  //const double p2F = - 11./12.;
//  //const double p2sigma = - 28./81.;
//  const double p21 = pi16*(- 1./4.);
//  const double p22 = pi16*(3./8.);
//  const double p23 = pi16*(211./648.);
//  const double p24 = pi16*(- 5./9.);
//  const double p25 = pi16*(- 1./16.);
//  const double p26 = pi16*( - 1./8.);
//  const double p27 = pi16*17./1296.;
//
//  double l1r,l2r,l3r,l4r,l5r,l6r,l7r;//,h1r,h2r,h3r;
//  double x = pi16*mk2bar2/pow(massin.fpi,2);
//
//  l1r = a1r + x*(p01+p11*lk+p21*pow(lk,2));
//  l2r = a2r + x*(p02+p12*lk+p22*pow(lk,2));
//  l3r = a3r + x*(p03+p13*lk+p23*pow(lk,2));
//  l4r = a4r + x*(p04+p14*lk+p24*pow(lk,2));
//  l5r = a5r + x*(p05+p15*lk+p25*pow(lk,2));
//  l6r = a6r + x*(p06+p16*lk+p26*pow(lk,2));
//  l7r = a7r + x*(p07+p17*lk+p27*pow(lk,2));
//
//  // turning them in barred versions
//  double lp = log(mp2/mu2);
//  const double pi32 = 32.*M_PI*M_PI;
//
//  l1r = -lp+pi32/g1*l1r;
//  l2r = -lp+pi32/g2*l2r;
//  l3r = -lp+pi32/g3*l3r;
//  l4r = -lp+pi32/g4*l4r;
//  l5r = -lp+pi32/g5*l5r;
//  l6r = -lp+pi32/g6*l6r;
//  return libar(l1r,l2r,l3r,l4r,l5r,l6r);
//}
