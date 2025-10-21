// fourvector.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.


// A small fourvector class /////////////////////////////////////////////////
// allows for adding, subtracting, scalar product and some implicit conversions
// as well as multiplying and dividing by a number
// NO checks done for consistency !!!

#include <vector>
#include <iostream>
#include <iomanip>

#include "fourvector.h"

fourvector::fourvector(const double a0,const double a1,const double a2,
		       const double a3){
  v[0]=a0;
  v[1]=a1;
  v[2]=a2;
  v[3]=a3;
}

fourvector::fourvector(const std::vector<double> vv){
  if(vv.size() != 4)
    std::cout << "WARNING attempt to set fourvector with wrong number of elements\n";
  v = vv;
}

fourvector fourvector::operator+(const fourvector vv) const{
  return fourvector(v[0]+vv.v[0],v[1]+vv.v[1],v[2]+vv.v[2],v[3]+vv.v[3]);
}

fourvector& fourvector::operator+=(const fourvector vv){
  *this=*this+vv;
  return *this;
}

fourvector fourvector::operator-(const fourvector vv) const{
  return fourvector(v[0]-vv.v[0],v[1]-vv.v[1],v[2]-vv.v[2],v[3]-vv.v[3]);
}

fourvector& fourvector::operator-=(const fourvector vv){
  *this=*this-vv;
  return *this;
}

fourvector fourvector::operator-(void) const{
  return fourvector(-v[0],-v[1],-v[2],-v[3]);
}

fourvector fourvector::operator*(const double aa) const{
  return fourvector(aa*v[0],aa*v[1],aa*v[2],aa*v[3]);
}

double fourvector::operator*(const fourvector vv) const{
  return vv.v[0]*v[0]-vv.v[1]*v[1]-vv.v[2]*v[2]-vv.v[3]*v[3];
}

fourvector operator*(const double aa, const fourvector vv){
  return fourvector(aa*vv.v[0],aa*vv.v[1],aa*vv.v[2],aa*vv.v[3]);
}

std::ostream & operator<<(std::ostream & os,const fourvector & bb){ 
  // does the output
  os << bb.v[0]<<' '<<bb.v[1]<<' '<<bb.v[2]<<' '<<bb.v[3];
  return os;
}

std::istream & operator>>(std::istream & is, fourvector & bb){ 
  // does the input from four numbers
  is >> bb.v[0]>>bb.v[1]>>bb.v[2]>>bb.v[3];
  return is;
}

fourvector::operator std::vector<double>() const{
  return {v[0],v[1],v[2],v[3]};
}

double fourvector::out(const int n) const{
  return v[n];
}

fourvector fourvector::operator/(const double aa) const{
  return fourvector(v[0]/aa,v[1]/aa,v[2]/aa,v[3]/aa);
}

void fourvector::set(const int n, const double xn){
  v[n] = xn;
}

void fourvector::set(const double xn, const int n){
  v[n] = xn;
}

void fourvector::lower(void){
  for(int i=1; i<4;i++){
    v[i] = -v[i];
  }
}
