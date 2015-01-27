// testfinitevolumeoneloopintegrals.cc is part of 
// the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include<iostream>
#include<iomanip>
#include<cmath>
#include "finitevolumeoneloopintegrals.h"

const double pi16 = 1./pow(4.*M_PI,2);

int main(void){
  using namespace std;
  double mpi = 0.135;
  double xl = 2./mpi;
  double msq = mpi*mpi;
  cout <<"A  number of results for mpi L=2, theta and Bessel\n";
  cout << "AbVacc maxsum AbVt AbVb A22bVt A22bVb A23bVt A23bVb BbVt BbVb\n"; 
  double epst[8] = {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8};
  int nb[8] = {1,2,5,10,20,40,60,100};
  for(int i=0; i<8; i++){
    setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
    setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
    cout << epst[i]<<' '<<nb[i]<<' '
	 <<AbVt(msq,xl) <<' '<<AbVb(msq,xl)<<' '
	 <<A22bVt(msq,xl) <<' '<<A22bVb(msq,xl)<<' '
	 <<A23bVt(msq,xl) <<' '<<A23bVb(msq,xl)<<' '
	 <<BbVt(msq,xl) <<' '<<BbVb(msq,xl)<<'\n';
  }
  cout <<"checking the relation 4*A22bV+3*A23bV=msq*AbV\n";
  cout << 4.*A22bVt(msq,xl)+3.*A23bVt(msq,xl) << ' '<<msq*AbVt(msq,xl)<<'\n';
  cout << 4.*A22bVb(msq,xl)+3.*A23bVb(msq,xl) << ' '<<msq*AbVb(msq,xl)<<'\n';
  return 0;
}
