// testinputs.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// test the masses etc in and output routines in inputs.cc

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "inputsnf.h"

int main(void){
  // a number of different input values
  double  mpi = 0.15;
  double  mk = 0.5;
  double  meta = 0.55;
  double fpi = 0.95;
  double  mu = 0.9;
  double mp0 = 0.12;
  double mk0 = 0.41;
  double f0 = 0.88;
  double B0mhat = 0.011;
  double B0ms = 0.28;
  // creating a set with standard input values and puttin it out;
  quarkmassnf mass1;
  cout << "First mass case\n" << mass1;
  vector<double> B0mi={0.1,0.2,0.3,0.4,0.5};
  quarkmassnf mass2(B0mi,0.85,0.65);
  cout << "Second mass case\n" << mass2;
  quarkmassnf mass3 = mass2;
  mass3.setB0mq(1,0.15);
  vector<double> B0m;
  mass3.out(B0m);
  cout << mass3;
  for(int i=0; i<B0m.size();i++) cout << B0m[i]<<'\n';
  // writing them out to a file and then reading back in
  ofstream ouf("test.dat");
  ouf << mass1 << mass2;
  ouf.close();
  quarkmassnf mass4,mass5;
  ifstream inf("test.dat");
  inf >> mass4 >> mass5;
  inf.close();
  cout << "written to and read from file\n" << mass4 << mass5;
  // checking equality of two masses (note only to about 7 digits
  // to avoid calculated masses to become unequal
  if( mass1 == mass4){
    cout << "mass4 equal to mass1\n";}
  else{
    cout << "mass4 not equal to mass1\n";}
  if( mass1 == mass5){
    cout << "mass5 equal to mass1\n";}
  else{
    cout << "mass5 not equal to mass1\n";}


  cout <<"Note the numbering 1,..,nq in the setting and output\n"
       << "but 0,...,nq-1 in the vector storing them\n";
  B0m = mass3.getB0mq();
  int nq = mass3.getnq();
  for(int i=0; i< nq;i++){
    cout << mass3.getB0mq(i+1)<<' '<<B0m[i]<<'\n';
  }
  return 0;
}
