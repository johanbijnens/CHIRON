// testfourvector.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include<iostream>
#include<fstream>
#include "fourvector.h"
#include<vector>

int main(void){
  using namespace std;

  fourvector p1(1.,0.,0.,0.), p2={2.,0.5,0.1,0.2};
  cout <<"outputting two fourvectors\n";
  cout << "p1 : " <<p1<<'\n'
       << "p2 : " <<p2<<'\n';
  vector<double> pp3={1.,2.,3.};
  cout <<"This should give a warning\n";
  fourvector p3 = pp3;
  cout <<"The sum p1+p2        : "<<p1+p2<<'\n';
  cout <<"The difference p1-p2 : "<<p1-p2<<'\n';
  cout <<" -p2                 : "<<(-p2)<<'\n';
  cout <<" p2*0.1              : "<< p2*0.1<<'\n';
  cout <<" 0.1*p2              : "<< 0.1*p2<<'\n';

  ofstream outfile;
  outfile.open("tempfourvector.dat");
  outfile << p2 <<'\n';
  outfile.close();
  fourvector p4;
  ifstream infile;
  infile.open("tempfourvector.dat");
  infile >> p4;
  infile.close();
  cout << "Written to file and read back in : "<<p4<<'\n';
  vector<double> pvect = p4;
  cout << "As a vector of doubles and looking at one component\n";
  cout << pvect[1] <<' '<<p4.out(1)<<'\n';
  fourvector p5(0.1,0.2,0.3,0.4);
  cout << "p5 and p5/5\n";
  cout << p5 << '\n'<< p5/5 <<'\n';
  cout << "scalar product: p5.p5 and p1.p2\n";
  cout << p5*p5 <<' '<<p1*p2<<'\n';

  cout <<"lowering indices\n";
  cout << p5 <<'\n';
  p5.lower();
  cout << p5 <<'\n';

  return 0;
}
