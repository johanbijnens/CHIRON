// testinputsnf2.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2019 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// test the masses etc in and output routines in inputs.cc

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "inputsnf2.h"

int main(void){
  // a number of different input values
  double  mpi = 0.15;
  double fpi = 0.95;
  double  mu = 0.9;
  double mp = 0.12;
  double f = 0.88;
  double Bmhat = 0.011;
  // creating a set with standard input values and puttin it out;
  physmassnf2 mass1;
  lomassnf2 mass1l;
  quarkmassnf2 mass1q;
  cout << "First mass case\n" << mass1 << mass1l << mass1q;
  // creating a set with the other inputs above and then changin one
  physmassnf2 mass2(mpi,fpi,mu);
  lomassnf2 mass2l(mp,f,mu);
  quarkmassnf2 mass2q(Bmhat,f,mu);
  cout << "Second mass case\n" << mass2 << mass2l << mass2q;
  physmassnf2 mass3 = mass2;
  lomassnf2 mass3l = mass2l;
  quarkmassnf2 mass3q = mass2q;
  // checking the conversion
  lomassnf2 mass3lp = mass2q;
  quarkmassnf2 mass3qp = mass2l;
  cout << "Third mass case\n" << mass3 << mass3l << mass3q
       << mass3lp << mass3qp;
  // writing them out to a file and then reading back in
  ofstream ouf("test.dat");
  ouf << mass2 << mass2l << mass2q;
  ouf.close();
  physmassnf2 mass4; lomassnf2 mass4l; quarkmassnf2 mass4q;
  ifstream inf("test.dat");
  inf >> mass4 >> mass4l >> mass4q;
  inf.close();
  cout << "written to and read from file\n" << mass4 << mass4l << mass4q;
  // reading out the data
  mass4.out(mpi,fpi,mu);
  cout <<"mpi,fpi,mu\n"
       <<mpi<<' '<<fpi<<' '<<mu<<endl;
  mass4l.out(mp,f,mu);
  cout << "mp, f, mu\n"
       << mp<<' '<<f<<' '<<mu<<'\n';
  mass4q.out(Bmhat,f,mu);
  cout << "Bmhat, f, mu\n"
       << Bmhat<<' '<<f<<' '<<mu<<'\n';
  // same but now read out with the separate functions, resetting fpi, f
  mass4.setfpi(0.12);
  mass4l.setf(0.095);
  mass4l.setmu(0.81);
  mass4q.setf(0.084);
  mpi = mass4.getmpi();
  fpi = mass4.getfpi();
  mu  = mass4.getmu();
  cout <<"mpi,fpi,mu\n"
       <<mpi<<' '<<fpi<<' '<<mu<<endl;
  mp = mass4l.getmp();
  f = mass4l.getf();
  mu = mass4l.getmu();
  cout << "mp, f, mu\n"
       << mp<<' '<<f<<' '<<mu<<'\n';
  Bmhat = mass4q.getBmhat();
  f = mass4q.getf();
  mu = mass4q.getmu();
  cout << "Bmhat, f, mu\n"
       << Bmhat<<' '<<f<<' '<<mu<<'\n';
  // checking equality of two masses (note only to about 7 digits
  // to avoid calculated masses to become unequal
  if( mass1 == mass4){
    cout << "mass4 equal to mass1\n";}
  else{
    cout << "mass4 not equal to mass1\n";}
  if( mass1l == mass4l){
    cout << "mass4l equal to mass1l\n";}
  else{
    cout << "mass4l not equal to mass1l\n";}
  if( mass1q == mass4q){
    cout << "mass4q equal to mass1q\n";}
  else{
    cout << "mass4q not equal to mass1q\n";}
  physmassnf2 mass5; lomassnf2 mass5l; quarkmassnf2 mass5q;
  if( mass1 == mass5){
    cout << "mass5 equal to mass1\n";}
  else{
    cout << "mass5 not equal to mass1\n";}
  if( mass1l == mass5l){
    cout << "mass5l equal to mass1l\n";}
  else{
    cout << "mass5l not equal to mass1l\n";}
  if( mass1q == mass5q){
    cout << "mass5q equal to mass1q\n";}
  else{
    cout << "mass5q not equal to mass1q\n";}

  return 0;
}
