// testLinfunction.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2025 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.
#include<string>
#include<cmath>
#include<iostream>
#include<iomanip>
#include "jbnumlib.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// tests the Li2,3,4,5 implementations for a few points
struct nnint{
  int nn;
  dcomplex zz;
};

dcomplex linint(const dcomplex z,void* ap){
 nnint nn = *((nnint*)(ap));
 return pow(log(z),nn.nn-1)/(z-1./nn.zz);
}

dcomplex jbdlinint(int n, dcomplex z){
  nnint nn={n,z};
  void* ap = &nn;
  int ifacn = 1;
  for(int i=2;i<n;i++) ifacn *= i;
  return pow(-1.,n)/double(ifacn)*jbwquad15(linint,dcomplex(0.),dcomplex(1.),1e-8,ap);
}

dcomplex jbdlinpow2(int n, dcomplex z){
dcomplex aa[57],zn=1.;
//std::cout << zn << '\n';
  const int nmax = 56;
  for(int k=1; k <= nmax; k++){
    zn *= z;
    aa[k]=zn/pow(double(k),n);
   }
  dcomplex result=0.;
  for(int k=nmax; k > 0; k--){
    result += aa[k];
  }
//for(int k=1; k<=nmax; k++){
//std:: cout << aa[k] <<'\n';
//}
  return result;
}

int main(void){
    std::string test[3]={"test around 0\n","test larger than 1\n","test around 1\n"};
    dcomplex x[3] = {dcomplex(-0.5,0.2),dcomplex(-1.5,.2),dcomplex(0.51,0.1)};

    for(int i=0;i<3;
        i++){
        std::cout << test[i]; 
        std::cout << "Li2: "<<jbdli2(x[i])<<' '<<jbdlinpow(2,x[i])<<' '<<abs(jbdli2(x[i])-jbdlinpow(2,x[i]))<<jbdlinint(2,x[i])<<'\n';
        std::cout << "Li3: "<<jbdli3(x[i])<<' '<<jbdlinpow(3,x[i])<<' '<<abs(jbdli3(x[i])-jbdlinpow(3,x[i]))<<jbdlinint(3,x[i])<<'\n';
        std::cout << "Li4: "<<jbdli4(x[i])<<' '<<jbdlinpow(4,x[i])<<' '<<abs(jbdli4(x[i])-jbdlinpow(4,x[i]))<<jbdlinint(4,x[i])<<'\n';
        std::cout << "Li5: "<<jbdli5(x[i])<<' '<<jbdlinpow(5,x[i])<<' '<<abs(jbdli5(x[i])-jbdlinpow(5,x[i]))<<jbdlinint(5,x[i])<<'\n';
        std::cout << "Li6: "<<jbdli6(x[i])<<' '<<jbdlinpow(6,x[i])<<' '<<abs(jbdli6(x[i])-jbdlinpow(6,x[i]))<<jbdlinint(6,x[i])<<'\n';
    //    std::cout << "Li5: "<<jbdlinpow2(5,x[i])<<'\n';// power series with fewer terms
    }
    return 0;
}