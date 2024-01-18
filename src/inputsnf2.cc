// inputsnf2.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2019 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.


// implementation of input classes
// and setting them to particular values

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "inputsnf2.h"

//++++++++ physmassnf2 ++++++++++++++++++++++++++++++++++++++++++++++++++++++

physmassnf2::~physmassnf2(void){}

physmassnf2::physmassnf2(const double mpiin, const double fpiin,
		   const double muin){
  mpi = mpiin;
  fpi = fpiin;
  mu = muin;
}

void physmassnf2::setmpi(const double mpiin){
  mpi = mpiin;
}
void physmassnf2::setfpi(const double fpiin){
  fpi = fpiin;
}
void physmassnf2::setmu(const double muin){
  mu = muin;
}

void physmassnf2::out(double &mpiout, double &fpiout, double &muout) const{
  mpiout = mpi;
  fpiout = fpi;
  muout = mu;
}

double physmassnf2::getmpi(void ) const{
  return mpi;
}
double physmassnf2::getfpi(void) const{
  return fpi;
}
double physmassnf2::getmu(void) const{
  return mu;
}

std::ostream & operator<<(std::ostream & os, const physmassnf2 & bb){
  os   << setprecision(15) << fixed
       << "# mpi :  "<< bb.mpi << "\n"
       << "# fpi :  "<<bb.fpi <<"\n"
       << "# mu  :  "<<bb.mu <<"\n";
  os.unsetf(std::ios_base::floatfield);
  return os;
}

std::istream & operator>>(std::istream & is, physmassnf2 & mass){
  string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in massnf2 output
  // up to the :
  string temps2;
  double mpi,fpi,mu;
  std::getline(is,temps2,':');// reads in characters up to ':'
  is >> mpi;
  std::getline(is,temps); // to remove end of line
  if (temps2 != "# mpi ") std::cout << "trouble reading in massnf2 mpi "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> fpi;
  std::getline(is,temps);
  if (temps2 != "# fpi ") std::cout <<"trouble reading in massnf2 fpi "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mu;
  std::getline(is,temps);
  if (temps2 != "# mu  ") std::cout << "trouble reading in massnf2 mu "<< temps2<<'\n';
  mass = physmassnf2(mpi,fpi,mu);
  return is;
}

bool operator==(const physmassnf2 mass1,const physmassnf2 mass2){
  const double massprecision = 1e-7;
  if ( (abs((mass1.mpi-mass2.mpi)/mass1.mpi) < massprecision) &&
        (abs((mass1.fpi-mass2.fpi)/mass1.fpi) < massprecision) &&
       (abs((mass1.mu-mass2.mu)/mass1.mu) < massprecision) ) return true;
  return false;
}

//+++++++++++lomassnf2+++++++++++++++++++++++++++++++++++++++++++++++++++++++

lomassnf2::~lomassnf2(void){}

lomassnf2::lomassnf2(const double mp0in, const double f0in,
		     const double muin){
  mp = mp0in;
  f = f0in;
  mu = muin;
}

lomassnf2::lomassnf2(const quarkmassnf2 mass){
  mp = sqrt(2.*mass.getBmhat());
  f = mass.getf();
  mu = mass.getmu();
}
void lomassnf2::setmp(const double mp0in){
  mp = mp0in;
}
void lomassnf2::setf(const double f0in){
  f = f0in;
}
void lomassnf2::setmu(const double muin){
  mu = muin;
}
void lomassnf2::out(double &mp0out, double &f0out, double &muout) const{
  mp0out = mp;
  f0out = f;
  muout = mu;
}

double lomassnf2::getmp(void ) const{
  return mp;
}
double lomassnf2::getf(void) const{
  return f;
}
double lomassnf2::getmu(void) const{
  return mu;
}

std::ostream & operator<<(std::ostream & os, const lomassnf2 & bb){
  os   << setprecision(15) << fixed
       << "# mp  :  "<< bb.mp << "\n"
       << "# f   :  "<<bb.f  <<"\n"
       << "# mu  :  "<<bb.mu <<"\n";
  os.unsetf(std::ios_base::floatfield);
  return os;
}

std::istream & operator>>(std::istream & is, lomassnf2 & mass){
  string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  string temps2;
  double mp,f,mu;
  std::getline(is,temps2,':');// reads in characters up to ':'
  is >> mp;
  std::getline(is,temps); // to remove end of line
  if (temps2 != "# mp  ") std::cout << "trouble reading in lomassnf2 mp "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> f;
  std::getline(is,temps);
  if (temps2 != "# f   ") std::cout <<"trouble reading in lomassnf2 f "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mu;
  std::getline(is,temps);
  if (temps2 != "# mu  ") std::cout << "trouble reading in lomassnf2 mu "<< temps2<<'\n';
  mass = lomassnf2(mp,f,mu);
  return is;
}

bool operator==(const lomassnf2 mass1,const lomassnf2 mass2){
  const double massprecision = 1e-7;
  if ( (abs((mass1.mp-mass2.mp)/mass1.mp) < massprecision) &&
       (abs((mass1.f-mass2.f)/mass1.f) < massprecision) &&
       (abs((mass1.mu-mass2.mu)/mass1.mu) < massprecision) ) return true;
  return false;
}

//+++++++++++quarkmassnf2+++++++++++++++++++++++++++++++++++++++++++++++++++++++

quarkmassnf2::~quarkmassnf2(void){}

quarkmassnf2::quarkmassnf2(const double B0mhatin,
	       const double f0in, const double muin){
  Bmhat = B0mhatin;
  f = f0in;
  mu = muin;
}
quarkmassnf2::quarkmassnf2(const lomassnf2 mass){
  Bmhat = pow(mass.getmp(),2)/2.;
  f = mass.getf();
  mu = mass.getmu();
}
void quarkmassnf2::setBmhat(const double B0mhatin){
  Bmhat = B0mhatin;
}
void quarkmassnf2::setf(const double f0in){
  f = f0in;
}
void quarkmassnf2::setmu(const double muin){
  mu = muin;
}
void quarkmassnf2::out(double &B0mhatout,
		   double &f0out, double &muout) const{
  B0mhatout = Bmhat;
  f0out = f;
  muout = mu;
}

double quarkmassnf2::getBmhat(void ) const{
  return Bmhat;
}
double quarkmassnf2::getf(void) const{
  return f;
}
double quarkmassnf2::getmu(void) const{
  return mu;
}

std::ostream & operator<<(std::ostream & os, const quarkmassnf2 & bb){
  os   << setprecision(15) << fixed
       << "# Bmhat  :  "<< bb.Bmhat << "\n"
       << "# f      :  "<<bb.f <<"\n"
       << "# mu     :  "<<bb.mu <<"\n";
  os.unsetf(std::ios_base::floatfield);
  return os;
}

std::istream & operator>>(std::istream & is, quarkmassnf2 & mass){
  string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  string temps2;
  double B0mhat,f0,mu;
  std::getline(is,temps2,':');// reads in characters up to ':'
  is >> B0mhat;
  std::getline(is,temps); // to remove end of line
  if (temps2 != "# Bmhat  ") std::cout << "trouble reading in quarkmassnf2 Bmhat "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> f0;
  std::getline(is,temps);
  if (temps2 != "# f      ") std::cout <<"trouble reading in quarkmassnf2 f "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mu;
  std::getline(is,temps);
  if (temps2 != "# mu     ") std::cout << "trouble reading in quarkmassnf2 mu "<< temps2<<'\n';
  mass = quarkmassnf2(B0mhat,f0,mu);
  return is;
}

bool operator==(const quarkmassnf2 mass1,const quarkmassnf2 mass2){
  const double massprecision = 1e-7;
  if ( (abs((mass1.Bmhat-mass2.Bmhat)/mass1.Bmhat) < massprecision) &&
       (abs((mass1.f-mass2.f)/mass1.f) < massprecision) &&
       (abs((mass1.mu-mass2.mu)/mass1.mu) < massprecision) ) return true;
  return false;
}

