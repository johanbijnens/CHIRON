// testli.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2017 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// test the li and libar classes

#include<fstream>
#include<iostream>
#include "linf2.h"


int main(void){
  using namespace std;
  // creating one with zero entries
  li liset1;
  cout << "A zero set\n";
  cout << liset1;
  // input with everything specified
  cout << "one with everything nonzero\n";
  li liset2(1e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3,10e-3,0.78,"Set 2");
  cout << liset2;
  // resetting liset1 to some values
  liset1.setli(1,1e-3);
  liset1.setli(2,2e-3);
  liset1.setli(3,3e-3);
  liset1.setli(4,4e-3);
  liset1.setli(5,5e-3);
  liset1.setli(6,6e-3);
  liset1.setli(7,7e-3);
  liset1.setli(8,8e-3);
  liset1.setli(9,9e-3);
  liset1.setli(10,10e-3);
  liset2.setmu(0.9);
  liset2.setname("Set 1 redefined");
  cout << "First set fully redefined\n"<< liset1;
  // sum of two
  li liset3 = liset1+liset2;
  cout << liset3;
  // difference  of two
  li liset4 = liset1+liset2;
  cout << liset4;
  // multiplying with a number front and back
  li liset5 = 2.*liset1;
  liset5.setname("Multiplied from front");
  li liset6 = liset1*3.;
  liset6.setname("Multiplied from back");
  cout << "multiplying with a number front and back\n";
  cout << liset5;
  cout << liset6;
  // writing to a stream and reading it back in, here with a filestream
  ofstream ouf("temp.dat");
  ouf << liset2;
  ouf.close();
  ifstream inf("temp.dat");
  li liset9;
  inf >> liset9;
  inf.close();
  cout << "Written to and read back from file:\n";
  cout << liset9;
  // changing scale
  liset9.changescale(0.5);
  liset9.setname("Scale changed");
  cout << liset9;
  return 0;
}
