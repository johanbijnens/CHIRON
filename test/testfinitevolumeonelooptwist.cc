// testfinitevolumeonelooptwist.cc is part of the 
// CHIRON ChPT program collection
// Copyright (C) 2015-2016 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// note all in units of GeV
using namespace std;

#include<cmath>

#include "fourvector.h"
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumeonelooptwist.h"

int main(void){
  //const double hbarc = 0.1973269718;// hbarc in units of GeV fm  
  //double musq = 0.77*0.77;

  // setting an overall precision
  setprecisionfinitevolumeonelooptwistt(1e-12,1e-10);
  setprecisionfinitevolumeoneloopt(1e-12,1e-10);

  // printing out some numbers, pion and kaon masses
  double mpi = 0.1395;
  double mk = 0.495;
  double m1sq = mpi*mpi;
  double m2sq = mk*mk;
  fourvector thetazero = {0.,0.,0.,0.};
  fourvector theta1 = {0.,0.1,0.2,0.3};
  fourvector theta2 = {0.,M_PI/4.,-M_PI/3.,M_PI/12.};

  double xlm = 3.;
  double xl= xlm/sqrt(m1sq);

  cout << "mass m1 : "<<sqrt(m1sq)<<'\n';
  cout << "m1*L    : "<<xlm<<'\n';
  cout << "compare tadpoles for twist zero case with earlier\n";
  double x1 = AbVt(m1sq,xl);
  double x2 = AbVtwistt(1,0,11,m1sq,xl,thetazero);
  cout << "AbV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = BbVt(m1sq,xl);
  x2 = AbVtwistt(2,0,11,m1sq,xl,thetazero);
  cout << "BbV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = CbVt(m1sq,xl);
  x2 = AbVtwistt(3,0,11,m1sq,xl,thetazero);
  cout << "CbV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = DbVt(m1sq,xl);
  x2 = AbVtwistt(4,0,11,m1sq,xl,thetazero);
  cout << "DbV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = A22bVt(m1sq,xl);
  x2 = AbVtwistt(1,22,11,m1sq,xl,thetazero);
  cout << "A22bV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = B22bVt(m1sq,xl);
  x2 = AbVtwistt(2,22,11,m1sq,xl,thetazero);
  cout << "B22bV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = C22bVt(m1sq,xl);
  x2 = AbVtwistt(3,22,11,m1sq,xl,thetazero);
  cout << "C22bV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = D22bVt(m1sq,xl);
  x2 = AbVtwistt(4,22,11,m1sq,xl,thetazero);
  cout << "D22bV " <<' '<<x1<<' '<<x2<<'\n';
  // note the - sign is from t^{11}
  x1 = -A23bVt(m1sq,xl);
  x2 = AbVtwistt(1,23,11,m1sq,xl,thetazero);
  cout << "A23bV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = -B23bVt(m1sq,xl);
  x2 = AbVtwistt(2,23,11,m1sq,xl,thetazero);
  cout << "B23bV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = -C23bVt(m1sq,xl);
  x2 = AbVtwistt(3,23,11,m1sq,xl,thetazero);
  cout << "C23bV " <<' '<<x1<<' '<<x2<<'\n';
  x1 = -D23bVt(m1sq,xl);
  x2 = AbVtwistt(4,23,11,m1sq,xl,thetazero);
  cout << "D23bV " <<' '<<x1<<' '<<x2<<'\n';

  cout <<"\nComparing for some different twists and zero\n"
       <<"Numbers are for three twist vectors :\n"
       <<"theta   : "<<thetazero<<'\n' 
       <<"theta1  : "<<theta1<<'\n' 
       <<"theta2  : "<<theta2<<'\n';
  for(int i=1;i<=4;i++){
    cout << "AbVtwistt("<<i<<",0,0,m1sq,L,theta) :"
	 <<' '<<AbVtwistt(i,0,0,m1sq,xl,thetazero)
	 <<' '<<AbVtwistt(i,0,0,m1sq,xl,theta1)
	 <<' '<<AbVtwistt(i,0,0,m1sq,xl,theta2) << '\n';
  }
  for(int i=1;i<=4;i++){
    for (int j=1;j<=3;j++){
      cout << "AbVtwistt("<<i<<",2,"<<j<<",m1sq,L,theta) :"
	 <<' '<<AbVtwistt(i,2,j,m1sq,xl,thetazero)
	 <<' '<<AbVtwistt(i,2,j,m1sq,xl,theta1)
	 <<' '<<AbVtwistt(i,2,j,m1sq,xl,theta2) << '\n';
    }}
  for(int i=1;i<=4;i++){
    cout << "AbVtwistt("<<i<<",22,0,m1sq,L,theta) :"
	 <<' '<<AbVtwistt(i,22,0,m1sq,xl,thetazero)
	 <<' '<<AbVtwistt(i,22,0,m1sq,xl,theta1)
	 <<' '<<AbVtwistt(i,22,0,m1sq,xl,theta2) << '\n';
  }
  int comp[6]={11,12,13,22,23,33};
  for(int i=1;i<=4;i++){
    for (int j=0;j<6;j++){
      cout << "AbVtwistt("<<i<<",23,"<<comp[j]<<",m1sq,L,theta) :"
	 <<' '<<AbVtwistt(i,23,comp[j],m1sq,xl,thetazero)
	 <<' '<<AbVtwistt(i,23,comp[j],m1sq,xl,theta1)
	 <<' '<<AbVtwistt(i,23,comp[j],m1sq,xl,theta2) << '\n';
    }}

 ///////// testing relations ///////////////////////////////////////////
  cout <<"Testing 4 A_22^V+t^ii A_23ii^V-m1sq* A^V = A^V(one propagator less)\n";
  double x3,x4;
  for(int n=1;n <=4;n++){
    x1 = 4.*AbVtwistt(n,22,0,m1sq,xl,theta2);
    x2 = -AbVtwistt(n,23,11,m1sq,xl,theta2)-AbVtwistt(n,23,22,m1sq,xl,theta2)
         -AbVtwistt(n,23,33,m1sq,xl,theta2)-m1sq*AbVtwistt(n,0,0,m1sq,xl,theta2);
    x3 = 0.;
    if (n != 1){
    x3 = AbVtwistt(n-1,0,0,m1sq,xl,theta2);
    }
    x4 = x1+x2; 
    cout<<"n= "<<n<<" : "<<x3<<' '<<x4<<' '<<(x3/x4-1.)<<endl;
  }

  ////////////////// Bubble integrals /////////////////////////////////////
  cout <<"\nResults for bubble integrals\n";
  cout <<"mass m2 : "<<sqrt(m2sq)<<'\n';
  cout <<"twist: no twist, theta1 and second zero, theta1 and theta2\n";
  cout <<"p = {0.17,0.,0.,0.} + twist part\n";
  fourvector pext = {0.17,0.,0.,0.};
  fourvector p1 = pext+theta1/xl;
  fourvector p2 = pext+theta1/xl+theta2/xl;
  double psq = pext*pext;
  double p1sq = p1*p1;
  double p2sq = p2*p2;

  for(int i=1;i<=4;i++){
    cout << "BbVtwistt("<<i<<",0,0,m1sq,m2sq,psq,L,theta,p) :"
	 <<' '<<BbVtwistt(i,0,0,m1sq,m2sq,psq ,xl,thetazero,pext)
	 <<' '<<BbVtwistt(i,0,0,m1sq,m2sq,p1sq,xl,theta1,p1)
	 <<' '<<BbVtwistt(i,0,0,m1sq,m2sq,p2sq,xl,theta2,p2) << '\n';
  }
  for(int i=1;i<=4;i++){
    cout << "BbVtwistt("<<i<<",1,0,m1sq,m2sq,psq,L,theta,p) :"
	 <<' '<<BbVtwistt(i,1,0,m1sq,m2sq,psq ,xl,thetazero,pext)
	 <<' '<<BbVtwistt(i,1,0,m1sq,m2sq,p1sq,xl,theta1,p1)
	 <<' '<<BbVtwistt(i,1,0,m1sq,m2sq,p2sq,xl,theta2,p2) << '\n';
  }
  for(int i=1;i<=4;i++){
    for(int j=1;j<=3;j++){
      cout << "BbVtwistt("<<i<<",2,"<<j<<",m1sq,m2sq,psq,L,theta,p) :"
	 <<' '<<BbVtwistt(i,2,j,m1sq,m2sq,psq ,xl,thetazero,pext)
	 <<' '<<BbVtwistt(i,2,j,m1sq,m2sq,p1sq,xl,theta1,p1)
	 <<' '<<BbVtwistt(i,2,j,m1sq,m2sq,p2sq,xl,theta2,p2) << '\n';
    }}
  for(int i=1;i<=4;i++){
    cout << "BbVtwistt("<<i<<",21,0,m1sq,m2sq,psq,L,theta,p) :"
	 <<' '<<BbVtwistt(i,21,0,m1sq,m2sq,psq ,xl,thetazero,pext)
	 <<' '<<BbVtwistt(i,21,0,m1sq,m2sq,p1sq,xl,theta1,p1)
	 <<' '<<BbVtwistt(i,21,0,m1sq,m2sq,p2sq,xl,theta2,p2) << '\n';
  }
  for(int i=1;i<=4;i++){
    cout << "BbVtwistt("<<i<<",22,0,m1sq,m2sq,psq,L,theta,p) :"
	 <<' '<<BbVtwistt(i,22,0,m1sq,m2sq,psq ,xl,thetazero,pext)
	 <<' '<<BbVtwistt(i,22,0,m1sq,m2sq,p1sq,xl,theta1,p1)
	 <<' '<<BbVtwistt(i,22,0,m1sq,m2sq,p2sq,xl,theta2,p2) << '\n';
  }
  for(int i=1;i<=4;i++){
    for(int j=0;j<6;j++){
      cout << "BbVtwistt("<<i<<",23,"<<comp[j]<<",m1sq,m2sq,psq,L,theta,p) :"
	 <<' '<<BbVtwistt(i,23,comp[j],m1sq,m2sq,psq ,xl,thetazero,pext)
	 <<' '<<BbVtwistt(i,23,comp[j],m1sq,m2sq,p1sq,xl,theta1,p1)
	 <<' '<<BbVtwistt(i,23,comp[j],m1sq,m2sq,p2sq,xl,theta2,p2) << '\n';
    }}



  pext = p2;
  psq = pext*pext;
  double x5,x6;
  cout << "shift relations (k -> qext-k)\n";
  // shift relation BbV
  x1 = BbVtwistt(1,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = BbVtwistt(1,0,0,m2sq,m1sq,psq,xl,theta2,pext);
  cout << x1 <<' '<<x2<<' '<<(x1/x2-1.)<<endl;
  // shift relation B1bv
  x1 = BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = BbVtwistt(1,1,0,m2sq,m1sq,psq,xl,theta2,pext);
  x3 = BbVtwistt(1,0,0,m2sq,m1sq,psq,xl,theta2,pext);
  cout << x1+x2 << ' ' << x3 <<' '<<((x1+x2)/x3-1.)<<endl;
  // shift relations B2bv mu =1,2,3
  x1 =  BbVtwistt(1,2,1,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -BbVtwistt(1,2,1,m2sq,m1sq,psq,xl,theta2,pext);
  cout << x1 << ' ' << x2 <<' '<<(x1/x2-1.)<<endl;
  x1 =  BbVtwistt(1,2,2,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -BbVtwistt(1,2,2,m2sq,m1sq,psq,xl,theta2,pext);
  cout << x1 << ' ' << x2 <<' '<<(x1/x2-1.)<<endl;
  x1 =  BbVtwistt(1,2,3,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -BbVtwistt(1,2,3,m2sq,m1sq,psq,xl,theta2,pext);
  cout << x1 << ' ' << x2 <<' '<<(x1/x2-1.)<<endl;

  // one-loop bubbles
  cout << "relation B1, B2, B case propagators 1 \n";
  x1 =     psq*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -pext.out(1)*BbVtwistt(1,2,1,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(2)*BbVtwistt(1,2,2,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(3)*BbVtwistt(1,2,3,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*BbVtwistt(1,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 = 0.5*AbVtwistt(1,0,0,m2sq,xl,theta2)-0.5*AbVtwistt(1,0,0,m1sq,xl,theta1);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B1, B2, B case propagators 2 \n";
  x1 =     psq*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -pext.out(1)*BbVtwistt(2,2,1,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(2)*BbVtwistt(2,2,2,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(3)*BbVtwistt(2,2,3,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*BbVtwistt(2,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 = 0.5*BbVtwistt(1,0,0,m1sq,m2sq,psq,xl,theta1,pext)
      -0.5*AbVtwistt(2,0,0,m1sq,xl,theta1);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B1, B2, B case propagators 3 \n";
  x1 =     psq*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -pext.out(1)*BbVtwistt(3,2,1,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(2)*BbVtwistt(3,2,2,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(3)*BbVtwistt(3,2,3,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*BbVtwistt(3,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 = 0.5*AbVtwistt(2,0,0,m2sq,xl,theta2)
      -0.5*BbVtwistt(1,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B1, B2, B case propagators 4 \n";
  x1 =     psq*BbVtwistt(4,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -pext.out(1)*BbVtwistt(4,2,1,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(2)*BbVtwistt(4,2,2,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(3)*BbVtwistt(4,2,3,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*BbVtwistt(4,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 = 0.5*BbVtwistt(3,0,0,m1sq,m2sq,psq,xl,theta1,pext)
      -0.5*BbVtwistt(2,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B, propagators 1 \n";
  x1 = psq*BbVtwistt(1,21,0,m1sq,m2sq,psq,xl,theta1,pext)
       +4.*BbVtwistt(1,22,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -BbVtwistt(1,23,11,m1sq,m2sq,psq,xl,theta1,pext)
       -BbVtwistt(1,23,22,m1sq,m2sq,psq,xl,theta1,pext)
       -BbVtwistt(1,23,33,m1sq,m2sq,psq,xl,theta1,pext);
  x3 =  m1sq*BbVtwistt(1,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 =  AbVtwistt(1,0,0,m2sq,xl,theta2);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B, propagators 2 \n";
  x1 = psq*BbVtwistt(2,21,0,m1sq,m2sq,psq,xl,theta1,pext)
    +4.*BbVtwistt(2,22,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -BbVtwistt(2,23,11,m1sq,m2sq,psq,xl,theta1,pext)
       -BbVtwistt(2,23,22,m1sq,m2sq,psq,xl,theta1,pext)
       -BbVtwistt(2,23,33,m1sq,m2sq,psq,xl,theta1,pext);
  x3 =  m1sq*BbVtwistt(2,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 =  BbVtwistt(1,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B, propagators 3 \n";
  x1 = psq*BbVtwistt(3,21,0,m1sq,m2sq,psq,xl,theta1,pext)
    +4.*BbVtwistt(3,22,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -BbVtwistt(3,23,11,m1sq,m2sq,psq,xl,theta1,pext)
       -BbVtwistt(3,23,22,m1sq,m2sq,psq,xl,theta1,pext)
       -BbVtwistt(3,23,33,m1sq,m2sq,psq,xl,theta1,pext);
  x3 =  m1sq*BbVtwistt(3,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 =  AbVtwistt(2,0,0,m2sq,xl,theta2);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;


  cout << "relation B21, B22, B23, B, propagators 4 \n";
  x1 = psq*BbVtwistt(4,21,0,m1sq,m2sq,psq,xl,theta1,pext)
    +4.*BbVtwistt(4,22,0,m1sq,m2sq,psq,xl,theta1,pext);
  x2 = -BbVtwistt(4,23,11,m1sq,m2sq,psq,xl,theta1,pext)
       -BbVtwistt(4,23,22,m1sq,m2sq,psq,xl,theta1,pext)
       -BbVtwistt(4,23,33,m1sq,m2sq,psq,xl,theta1,pext);
  x3 =  m1sq*BbVtwistt(4,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 =  BbVtwistt(3,0,0,m1sq,m2sq,psq,xl,theta1,pext);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 0 propagators 1 \n";
  x1 = pext.out(0)*(psq*BbVtwistt(1,21,0,m1sq,m2sq,psq,xl,theta1,pext)
		+BbVtwistt(1,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = -pext.out(1)*BbVtwistt(1,23,10,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(2)*BbVtwistt(1,23,20,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(3)*BbVtwistt(1,23,30,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = pext.out(0)*
    0.5*(-m2sq+m1sq+psq)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 =  pext.out(0)*0.5*AbVtwistt(1,0,0,m2sq,xl,theta2);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;


  cout << "relation B21, B22, B23, B1, B2 mu = 0 propagators 2 \n";
  x1 = pext.out(0)*(psq*BbVtwistt(2,21,0,m1sq,m2sq,psq,xl,theta1,pext)
		+BbVtwistt(2,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = -pext.out(1)*BbVtwistt(2,23,10,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(2)*BbVtwistt(2,23,20,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(3)*BbVtwistt(2,23,30,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = pext.out(0)*
    0.5*(-m2sq+m1sq+psq)*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 =  0.5*pext.out(0)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;


  cout << "relation B21, B22, B23, B1, B2 mu = 0 propagators 3 \n";
  x1 = pext.out(0)*(psq*BbVtwistt(3,21,0,m1sq,m2sq,psq,xl,theta1,pext)
		+BbVtwistt(3,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = -pext.out(1)*BbVtwistt(3,23,10,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(2)*BbVtwistt(3,23,20,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(3)*BbVtwistt(3,23,30,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = pext.out(0)*
    0.5*(-m2sq+m1sq+psq)*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 =  pext.out(0)*0.5*AbVtwistt(2,0,0,m2sq,xl,theta2)
       -0.5*pext.out(0)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 0 propagators 4 \n";
  x1 = pext.out(0)*(psq*BbVtwistt(4,21,0,m1sq,m2sq,psq,xl,theta1,pext)
		+BbVtwistt(4,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = -pext.out(1)*BbVtwistt(4,23,10,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(2)*BbVtwistt(4,23,20,m1sq,m2sq,psq,xl,theta1,pext)
       -pext.out(3)*BbVtwistt(4,23,30,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = pext.out(0)*
    0.5*(-m2sq+m1sq+psq)*BbVtwistt(4,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x4 =  pext.out(0)*0.5*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext)
       -0.5*pext.out(0)*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 1 propagators 1 \n";
  x1 = pext.out(1)*(psq*BbVtwistt(1,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	       	+BbVtwistt(1,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(1,23,10,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(1,23,11,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(1,23,12,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(1,23,13,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(1)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
     +BbVtwistt(1,2,1,m1sq,m2sq,psq,xl,theta1,pext));
  x4 =  pext.out(1)*0.5*AbVtwistt(1,0,0,m2sq,xl,theta2)
    -0.5*AbVtwistt(1,2,1,m2sq,xl,theta2)-0.5*AbVtwistt(1,2,1,m1sq,xl,theta1);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 1 propagators 2 \n";
  x1 = pext.out(1)*(psq*BbVtwistt(2,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(2,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(2,23,10,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(2,23,11,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(2,23,12,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(2,23,13,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(1)*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(2,2,1,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = 0.5*(pext.out(1)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(1,2,1,m1sq,m2sq,psq,xl,theta1,pext))
       -0.5*AbVtwistt(2,2,1,m1sq,xl,theta1);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 1 propagators 3 \n";
  x1 = pext.out(1)*(psq*BbVtwistt(3,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(3,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(3,23,10,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(3,23,11,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(3,23,12,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(3,23,13,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(1)*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(3,2,1,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = pext.out(1)*0.5*AbVtwistt(2,0,0,m2sq,xl,theta2)
      -0.5*AbVtwistt(2,2,1,m2sq,xl,theta2)
      -0.5*(pext.out(1)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(1,2,1,m1sq,m2sq,psq,xl,theta1,pext));
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 1 propagators 4 \n";
  x1 = pext.out(1)*(psq*BbVtwistt(4,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(4,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(4,23,10,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(4,23,11,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(4,23,12,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(4,23,13,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(1)*BbVtwistt(4,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(4,2,1,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = 0.5*(pext.out(1)*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(3,2,1,m1sq,m2sq,psq,xl,theta1,pext))
      -0.5*(pext.out(1)*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(2,2,1,m1sq,m2sq,psq,xl,theta1,pext));
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 2 propagators 1 \n";
  x1 = pext.out(2)*(psq*BbVtwistt(1,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	       	+BbVtwistt(1,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(1,23,20,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(1,23,21,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(1,23,22,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(1,23,23,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(2)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
     +BbVtwistt(1,2,2,m1sq,m2sq,psq,xl,theta1,pext));
  x4 =  pext.out(2)*0.5*AbVtwistt(1,0,0,m2sq,xl,theta2)
    -0.5*AbVtwistt(1,2,2,m2sq,xl,theta2)-0.5*AbVtwistt(1,2,2,m1sq,xl,theta1);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 2 propagators 2 \n";
  x1 = pext.out(2)*(psq*BbVtwistt(2,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(2,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(2,23,20,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(2,23,21,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(2,23,22,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(2,23,23,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(2)*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(2,2,2,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = 0.5*(pext.out(2)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(1,2,2,m1sq,m2sq,psq,xl,theta1,pext))
       -0.5*AbVtwistt(2,2,2,m1sq,xl,theta1);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 2 propagators 3 \n";
  x1 = pext.out(2)*(psq*BbVtwistt(3,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(3,22,2,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(3,23,20,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(3,23,21,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(3,23,22,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(3,23,23,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(2)*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(3,2,2,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = pext.out(2)*0.5*AbVtwistt(2,0,0,m2sq,xl,theta2)
      -0.5*AbVtwistt(2,2,2,m2sq,xl,theta2)
      -0.5*(pext.out(2)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(1,2,2,m1sq,m2sq,psq,xl,theta1,pext));
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 2 propagators 4 \n";
  x1 = pext.out(2)*(psq*BbVtwistt(4,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(4,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(4,23,20,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(4,23,21,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(4,23,22,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(4,23,23,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(2)*BbVtwistt(4,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(4,2,2,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = 0.5*(pext.out(2)*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(3,2,2,m1sq,m2sq,psq,xl,theta1,pext))
      -0.5*(pext.out(2)*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(2,2,2,m1sq,m2sq,psq,xl,theta1,pext));
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 3 propagators 1 \n";
  x1 = pext.out(3)*(psq*BbVtwistt(1,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	       	+BbVtwistt(1,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(1,23,30,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(1,23,31,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(1,23,32,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(1,23,33,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(3)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
     +BbVtwistt(1,2,3,m1sq,m2sq,psq,xl,theta1,pext));
  x4 =  pext.out(3)*0.5*AbVtwistt(1,0,0,m2sq,xl,theta2)
    -0.5*AbVtwistt(1,2,3,m2sq,xl,theta2)-0.5*AbVtwistt(1,2,3,m1sq,xl,theta1);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 3 propagators 2 \n";
  x1 = pext.out(3)*(psq*BbVtwistt(2,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(2,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(2,23,30,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(2,23,31,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(2,23,32,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(2,23,33,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(3)*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(2,2,3,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = 0.5*(pext.out(3)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(1,2,3,m1sq,m2sq,psq,xl,theta1,pext))
       -0.5*AbVtwistt(2,2,3,m1sq,xl,theta1);
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 3 propagators 3 \n";
  x1 = pext.out(3)*(psq*BbVtwistt(3,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(3,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(3,23,30,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(3,23,31,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(3,23,32,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(3,23,33,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(3)*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(3,2,3,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = pext.out(3)*0.5*AbVtwistt(2,0,0,m2sq,xl,theta2)
      -0.5*AbVtwistt(2,2,3,m2sq,xl,theta2)
      -0.5*(pext.out(3)*BbVtwistt(1,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(1,2,3,m1sq,m2sq,psq,xl,theta1,pext));
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  cout << "relation B21, B22, B23, B1, B2 mu = 3 propagators 4 \n";
  x1 = pext.out(3)*(psq*BbVtwistt(4,21,0,m1sq,m2sq,psq,xl,theta1,pext)
	           +BbVtwistt(4,22,0,m1sq,m2sq,psq,xl,theta1,pext));
  x2 = pext.out(0)*BbVtwistt(4,23,30,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(1)*BbVtwistt(4,23,31,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(2)*BbVtwistt(4,23,32,m1sq,m2sq,psq,xl,theta1,pext)
      -pext.out(3)*BbVtwistt(4,23,33,m1sq,m2sq,psq,xl,theta1,pext);
  x3 = 0.5*(-m2sq+m1sq+psq)*
    (pext.out(3)*BbVtwistt(4,1,0,m1sq,m2sq,psq,xl,theta1,pext)
            +BbVtwistt(4,2,3,m1sq,m2sq,psq,xl,theta1,pext));
  x4 = 0.5*(pext.out(3)*BbVtwistt(3,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(3,2,3,m1sq,m2sq,psq,xl,theta1,pext))
      -0.5*(pext.out(3)*BbVtwistt(2,1,0,m1sq,m2sq,psq,xl,theta1,pext)
                   +BbVtwistt(2,2,3,m1sq,m2sq,psq,xl,theta1,pext));
  x5 = x1+x2;
  x6 = x3+x4;
  cout<< x1<<' '<<x2<<' '<<x3<<' '<<x4<<' '
      <<x5<<' '<<x6<<' '<<(x5/x6-1.)<<endl;

  return 0;
}
