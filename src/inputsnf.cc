// inputsnf.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.


// implementation of input classes
// and setting them to particular values

#include <iostream>
#include <iomanip>
#include <string>
#include<vector>
#include <cmath>
#include "inputsnf.h"

//+++++++++++quarkmassnf+++++++++++++++++++++++++++++++++++++++++++++++++++++++

quarkmassnf::~quarkmassnf(void){}

quarkmassnf::quarkmassnf(const double f0in, const double muin, const int nqin){
  nq = nqin;
  f0 = f0in;
  mu = muin;
  B0mq.clear();
  for(int i=0;i<nq;i++) B0mq.push_back(0.);
}
quarkmassnf::quarkmassnf(const vector<double> B0mqin, const double f0in,
			 const double muin){
  nq = B0mqin.size();
  f0 = f0in;
  mu = muin;
  B0mq = B0mqin;
}
void quarkmassnf::setB0mq(const vector<double> B0mqin){
  nq = B0mqin.size();
  B0mq = B0mqin;
}
void quarkmassnf::setB0mq(const double B0mqin, const int i){
  if (i > nq || i < 1){
    std::cout << "attempt to set nonexisting quarkmass in quarkmassnf\n";
  }
  else{
    B0mq[i-1] = B0mqin;
  }
}
void quarkmassnf::setB0mq(const int i, const double B0mqin){
  if (i > nq || i < 1){
    std::cout << "attempt to set nonexisting quarkmass in quarkmassnf\n";
  }
  else{
    B0mq[i-1] = B0mqin;
  }
}

void quarkmassnf::setf0(const double f0in){
  f0 = f0in;
}
void quarkmassnf::setmu(const double muin){
  mu = muin;
}

void quarkmassnf::out(vector<double> &B0mqout) const{
  B0mqout = B0mq;
}

void quarkmassnf::out(vector<double> &B0mqout, double &f0out,
		      double &muout) const{
  B0mqout = B0mq;
  f0out = f0;
  muout = mu;
}

void quarkmassnf::out(vector<double> &B0mqout, double &f0out,
		      double &muout, int &nqout) const{
  B0mqout = B0mq;
  f0out = f0;
  muout = mu;
  nqout = nq;
}

int quarkmassnf::getnq(void) const{
  return nq;
}

vector<double> quarkmassnf::getB0mq(void) const{
  return B0mq;
}
double quarkmassnf::getB0mq(const int i) const{
  return B0mq[i-1];
}

double quarkmassnf::getf0(void) const{
  return f0;
}
double quarkmassnf::getmu(void) const{
  return mu;
}

std::ostream & operator<<(std::ostream & os, const quarkmassnf & bb){
  os   << setprecision(15) << fixed
       << "# nq     :  "<< bb.nq <<'\n';
  for(int i=0; i<bb.nq; i++){
    os  << "# B0mq : "<<(i+1)<<' '<< bb.B0mq[i] << "\n";
  }
  os   << "# f0     :  "<<bb.f0 <<"\n"
       << "# mu     :  "<<bb.mu <<"\n";
  os.unsetf(std::ios_base::floatfield);
  return os;
}

std::istream & operator>>(std::istream & is, quarkmassnf & mass){
  string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  string temps2;
  double B0mi,f0,mu;
  vector<double> B0mq;
  int nqin,iqin;
  std::getline(is,temps2,':');// reads in characters up to ':'
  is >> nqin;
  std::getline(is,temps); // to remove end of line
  if (temps2 != "# nq     ") std::cout << "trouble reading in quarkmassnf nq "<< temps2<<'\n';
  for(int i=0;i<nqin;i++){
    std::getline(is,temps2,':');
    is >> iqin >> B0mi;
    B0mq.push_back(B0mi);
    std::getline(is,temps);
    if (temps2 != "# B0mq ") std::cout << "trouble reading in quarkmassnf B0mq "
				       <<iqin<<' '<< temps2<<'\n';
  }
  std::getline(is,temps2,':');
  is >> f0;
  std::getline(is,temps);
  if (temps2 != "# f0     ") std::cout <<"trouble reading in quarkmassnf f0 "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mu;
  std::getline(is,temps);
  if (temps2 != "# mu     ") std::cout << "trouble reading in quarkmassnf mu "<< temps2<<'\n';
  mass = quarkmassnf(B0mq,f0,mu);
  return is;
}

bool operator==(const quarkmassnf mass1,const quarkmassnf mass2){
  const double massprecision = 1e-7;
  int nq = mass1.getnq();
  if( nq != mass2.getnq()) return false;
  vector<double> B01,B02;
  mass1.out(B01);
  mass2.out(B02);
  for(int i=0; i < nq;i++){
    if (abs((B01[i]-B02[i])/B01[i]) > massprecision) return false;
  }
  if(abs((mass1.getf0()-mass2.getf0())/mass1.getf0()) > massprecision) return false;
  if(abs((mass1.getmu()-mass2.getmu())/mass1.getmu()) > massprecision) return false;
  return true;
}
