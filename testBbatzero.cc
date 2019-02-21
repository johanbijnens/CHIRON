// g++ testBbatzero.cc -lchiron -ljbnumlib
// testing values at zero for infinute volume bubble integrals
#include<iostream>
#include "oneloopintegrals.h"
#include "quenchedoneloopintegrals.h"

int main(void){
  using namespace std;
  double m1sq = 0.135*0.135;
  double m2sq = 0.495*0.495;
  double mu2 = 0.77*0.77;
  double qsq0 = 0.;
  double qsq1 = 1e-3;
  double qsq2 = -qsq1;

  // normal case
  cout << Bb(m1sq,m2sq,qsq0,mu2)<<' '
       <<(Bb(m1sq,m2sq,qsq1,mu2)+Bb(m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B1b(m1sq,m2sq,qsq0,mu2)<<' '
       <<(B1b(m1sq,m2sq,qsq1,mu2)+B1b(m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B21b(m1sq,m2sq,qsq0,mu2)<<' '
       <<(B21b(m1sq,m2sq,qsq1,mu2)+B21b(m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B22b(m1sq,m2sq,qsq0,mu2)<<' '
       <<(B22b(m1sq,m2sq,qsq1,mu2)+B22b(m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B31b(m1sq,m2sq,qsq0,mu2)<<' '
       <<(B31b(m1sq,m2sq,qsq1,mu2)+B31b(m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B32b(m1sq,m2sq,qsq0,mu2)<<' '
       <<(B32b(m1sq,m2sq,qsq1,mu2)+B32b(m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << Bb(m1sq,m1sq,qsq0,mu2)<<' '
       <<(Bb(m1sq,m1sq,qsq1,mu2)+Bb(m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B1b(m1sq,m1sq,qsq0,mu2)<<' '
       <<(B1b(m1sq,m1sq,qsq1,mu2)+B1b(m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B21b(m1sq,m1sq,qsq0,mu2)<<' '
       <<(B21b(m1sq,m1sq,qsq1,mu2)+B21b(m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B22b(m1sq,m1sq,qsq0,mu2)<<' '
       <<(B22b(m1sq,m1sq,qsq1,mu2)+B22b(m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B31b(m1sq,m1sq,qsq0,mu2)<<' '
       <<(B31b(m1sq,m1sq,qsq1,mu2)+B31b(m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B32b(m1sq,m1sq,qsq0,mu2)<<' '
       <<(B32b(m1sq,m1sq,qsq1,mu2)+B32b(m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << '\n';
  cout << Bb(2,m1sq,m2sq,qsq0,mu2)<<' '
       <<(Bb(2,m1sq,m2sq,qsq1,mu2)+Bb(2,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << Bb(3,m1sq,m2sq,qsq0,mu2)<<' '
       <<(Bb(3,m1sq,m2sq,qsq1,mu2)+Bb(3,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B1b(2,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B1b(2,m1sq,m2sq,qsq1,mu2)+B1b(2,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B1b(3,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B1b(3,m1sq,m2sq,qsq1,mu2)+B1b(3,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B21b(2,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B21b(2,m1sq,m2sq,qsq1,mu2)+B21b(2,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B21b(3,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B21b(3,m1sq,m2sq,qsq1,mu2)+B21b(3,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B22b(2,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B22b(2,m1sq,m2sq,qsq1,mu2)+B22b(2,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B22b(3,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B22b(3,m1sq,m2sq,qsq1,mu2)+B22b(3,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B31b(2,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B31b(2,m1sq,m2sq,qsq1,mu2)+B31b(2,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B31b(3,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B31b(3,m1sq,m2sq,qsq1,mu2)+B31b(3,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B32b(2,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B32b(2,m1sq,m2sq,qsq1,mu2)+B32b(2,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << B32b(3,m1sq,m2sq,qsq0,mu2)<<' '
       <<(B32b(3,m1sq,m2sq,qsq1,mu2)+B32b(3,m1sq,m2sq,qsq2,mu2))/2.
       <<'\n';
  cout << '\n';
  cout << Bb(2,m1sq,m1sq,qsq0,mu2)<<' '
       <<(Bb(2,m1sq,m1sq,qsq1,mu2)+Bb(2,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << Bb(3,m1sq,m1sq,qsq0,mu2)<<' '
       <<(Bb(3,m1sq,m1sq,qsq1,mu2)+Bb(3,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B1b(2,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B1b(2,m1sq,m1sq,qsq1,mu2)+B1b(2,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B1b(3,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B1b(3,m1sq,m1sq,qsq1,mu2)+B1b(3,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B21b(2,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B21b(2,m1sq,m1sq,qsq1,mu2)+B21b(2,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B21b(3,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B21b(3,m1sq,m1sq,qsq1,mu2)+B21b(3,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B22b(2,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B22b(2,m1sq,m1sq,qsq1,mu2)+B22b(2,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B22b(3,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B22b(3,m1sq,m1sq,qsq1,mu2)+B22b(3,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B31b(2,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B31b(2,m1sq,m1sq,qsq1,mu2)+B31b(2,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B31b(3,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B31b(3,m1sq,m1sq,qsq1,mu2)+B31b(3,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B32b(2,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B32b(2,m1sq,m1sq,qsq1,mu2)+B32b(2,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';
  cout << B32b(3,m1sq,m1sq,qsq0,mu2)<<' '
       <<(B32b(3,m1sq,m1sq,qsq1,mu2)+B32b(3,m1sq,m1sq,qsq2,mu2))/2.
       <<'\n';

  return 0;

}
