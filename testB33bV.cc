// note all in units of GeV
// test AbV 33 case
//      BbV
// g++ -std=c++11 testB33bV.cc  -lchiron -ljbnumlib
using namespace std;

#include<cmath>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<string>

#include "fourvector.h"
#include "finitevolumeonelooptwist.h"

int intstring(string ain){
  int i;
  istringstream(ain) >> i;
  return i;
}

int main(void){
  //const double hbarc = 0.1973269718;// hbarc in units of GeV fm  
  //double musq = 0.77*0.77;

  // setting an overall precision
  setprecisionfinitevolumeonelooptwistt(1e-12,1e-9);

  // printing out some numbers, pion and kaon masses
  double mpi = 0.1395;
  double mk = 0.495;
  double m1sq = mpi*mpi;
  double m2sq = mk*mk;
  fourvector thetazero = {0.,0.,0.,0.};
  fourvector theta1 = {0.,0.1,0.2,0.3};
  fourvector theta2 = {0.,M_PI/4.,-M_PI/3.,M_PI/12.};
  //fourvector theta1 = thetazero;
  //fourvector theta2 = thetazero;

  double xlm = 1.5;
  double xl= xlm/sqrt(m1sq);
  fourvector qext = (theta2+theta1)/xl+fourvector(0.,2.*M_PI,0.,0.)/xl;
  double qsq = qext*qext+0.1*m1sq;
  double q0 = sqrt(qsq-qext*qext);
  cout << q0<<' '<<qext<<'\n';
  qext.set(0,q0);
  vector<double> qex = qext;

  cout <<"++++++++++++++ Input+++++++++++++++++++++++++++++\n";
  cout << "mass m1 : "<<sqrt(m1sq)<<'\n';
  cout << "mass m2 : "<<sqrt(m2sq)<<'\n';
  cout << "m1*L    : "<<xlm<<'\n';
  cout << "L       : "<<xl<<'\n';
  cout << "qsq     : "<<qsq<<' '<<qext*qext<<'\n';
  cout << "qext    : "<<qext<<'\n';
  cout << "qex : "<< qex[0]<<' '<<qex[1]<<' '<< qex[2]<<' '<< qex[3]<<' '<<'\n';


  string lmu[]={"0","1","2","3"};

  cout << "Output for A33bV: component  value\n";
  for(int i=0;i < 4;i++){
  for(int j=0;j < 4;j++){
  for(int k=0;k < 4;k++){
    int ncomponent;
    string incomp = lmu[i]+lmu[j]+lmu[k];
    istringstream(incomp) >> ncomponent;
    cout << incomp <<' '<<setw(3)<<ncomponent
	 <<' '<<AbVtwistt(2,33,ncomponent,m1sq,xl,theta2)<<'\n';
  }}}

  double d1,d2;

  cout << "Relations from multiplying with g_munu :\n";
  d1 =        AbVtwistt(1,33,100,m1sq,xl,theta2)
             -AbVtwistt(1,33,111,m1sq,xl,theta2)
             -AbVtwistt(1,33,122,m1sq,xl,theta2)
             -AbVtwistt(1,33,133,m1sq,xl,theta2);
  d2 = m1sq*AbVtwistt(1,2,1,m1sq,xl,theta2);
  cout << d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  d1 =        AbVtwistt(1,33,200,m1sq,xl,theta2)
             -AbVtwistt(1,33,211,m1sq,xl,theta2)
             -AbVtwistt(1,33,222,m1sq,xl,theta2)
             -AbVtwistt(1,33,233,m1sq,xl,theta2);
  d2 = m1sq*AbVtwistt(1,2,2,m1sq,xl,theta2);
  cout << d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  d1 =        AbVtwistt(1,33,300,m1sq,xl,theta2)
             -AbVtwistt(1,33,311,m1sq,xl,theta2)
             -AbVtwistt(1,33,322,m1sq,xl,theta2)
             -AbVtwistt(1,33,333,m1sq,xl,theta2);
  d2 = m1sq*AbVtwistt(1,2,3,m1sq,xl,theta2);
  cout << d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  for(int i=2;i<5;i++){
  d1 =        AbVtwistt(i,33,100,m1sq,xl,theta2)
             -AbVtwistt(i,33,111,m1sq,xl,theta2)
             -AbVtwistt(i,33,122,m1sq,xl,theta2)
             -AbVtwistt(i,33,133,m1sq,xl,theta2);
  d2 = AbVtwistt(i-1,2,1,m1sq,xl,theta2)+m1sq*AbVtwistt(i,2,1,m1sq,xl,theta2);
  cout << d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  d1 =        AbVtwistt(i,33,200,m1sq,xl,theta2)
             -AbVtwistt(i,33,211,m1sq,xl,theta2)
             -AbVtwistt(i,33,222,m1sq,xl,theta2)
             -AbVtwistt(i,33,233,m1sq,xl,theta2);
  d2 = AbVtwistt(i-1,2,2,m1sq,xl,theta2)+m1sq*AbVtwistt(i,2,2,m1sq,xl,theta2);
  cout << d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  d1 =        AbVtwistt(i,33,300,m1sq,xl,theta2)
             -AbVtwistt(i,33,311,m1sq,xl,theta2)
             -AbVtwistt(i,33,322,m1sq,xl,theta2)
             -AbVtwistt(i,33,333,m1sq,xl,theta2);
  d2 = AbVtwistt(i-1,2,3,m1sq,xl,theta2)+m1sq*AbVtwistt(i,2,3,m1sq,xl,theta2);
  cout << d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  }

  cout << "\n++++++++++++ Now for the Bubbles ++++++++++++++++++++++++\n";
  cout << "B31bV"<<BbVtwistt(4,31,0,m1sq,m2sq,qsq,xl,theta1,qext)<<'\n';
  cout << "B32bV"<<BbVtwistt(4,32,0,m1sq,m2sq,qsq,xl,theta1,qext)<<'\n';
  cout << "Components of B33bV: component nprop value\n";
  for(int i=0;i < 4;i++){
  for(int j=0;j < 4;j++){
  for(int k=0;k < 4;k++){
    int ncomp;
    istringstream(lmu[i]+lmu[j]+lmu[k]) >> ncomp;
    for(int nprop=1;nprop<5;nprop++){
    cout << lmu[i]+lmu[j]+lmu[k]
	 <<' '
	 <<BbVtwistt(nprop,33,ncomp,m1sq,m2sq,qsq,xl,theta1,qext)<<'\n';}
  }}}


  cout << "Relations due to shifting momenta and exchanging masses\n";
  int ichange[] = {1,3,2,4};
  cout <<"B31bV and nprop=1,...,4\n";
  for(int nprop=0;nprop<4;nprop++){
   d1 = BbVtwistt(nprop+1,31,0,m1sq,m2sq,qsq,xl,theta1,qext);
  d2 =    -BbVtwistt(ichange[nprop],31,0,m2sq,m1sq,qsq,xl,theta2,qext)
       +3.*BbVtwistt(ichange[nprop],21,0,m2sq,m1sq,qsq,xl,theta2,qext)
       -3.*BbVtwistt(ichange[nprop],1,0,m2sq,m1sq,qsq,xl,theta2,qext)
          +BbVtwistt(ichange[nprop],0,0,m2sq,m1sq,qsq,xl,theta2,qext);
  cout << d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  }
  cout <<"B32bV and nprop=1,...,4\n";
   for(int nprop=0;nprop<4;nprop++){
  d1 = BbVtwistt(nprop+1,32,0,m1sq,m2sq,qsq,xl,theta1,qext);
  d2 = BbVtwistt(ichange[nprop],22,0,m2sq,m1sq,qsq,xl,theta2,qext)
      -BbVtwistt(ichange[nprop],32,0,m2sq,m1sq,qsq,xl,theta2,qext);
  cout << d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  }

   cout <<"Components of B33bV\n and nprop=1,..,4\n";

  for(int i=0;i<4;i++){
  for(int j=0;j<4;j++){
  for(int k=0;k<4;k++){
    for(int nprop = 0; nprop<4;nprop++){
    int ij =  intstring(lmu[i]+lmu[j]);
    int ik =  intstring(lmu[i]+lmu[k]);
    int jk =  intstring(lmu[j]+lmu[k]);
    int ijk = intstring(lmu[i]+lmu[j]+lmu[k]);
    d1 =  BbVtwistt(nprop+1,33,ijk,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 =        - BbVtwistt(ichange[nprop],33,ijk,m2sq,m1sq,qsq,xl,theta2,qext)
          +qex[i]*BbVtwistt(ichange[nprop],23, jk,m2sq,m1sq,qsq,xl,theta2,qext)
          +qex[j]*BbVtwistt(ichange[nprop],23, ik,m2sq,m1sq,qsq,xl,theta2,qext)
          +qex[k]*BbVtwistt(ichange[nprop],23, ij,m2sq,m1sq,qsq,xl,theta2,qext)
   -qex[i]*qex[j]*BbVtwistt(ichange[nprop], 2,  k,m2sq,m1sq,qsq,xl,theta2,qext)
   -qex[i]*qex[k]*BbVtwistt(ichange[nprop], 2,  j,m2sq,m1sq,qsq,xl,theta2,qext)
   -qex[j]*qex[k]*BbVtwistt(ichange[nprop], 2,  i,m2sq,m1sq,qsq,xl,theta2,qext);
    cout << lmu[i]+lmu[j]+lmu[k]
	 <<' '<<d1<<' '<<d2<<' '<<d1/d2-1.
	 <<'\n';
    }
  }}}

  cout <<"B33bV reduces at qsq=0 to the correct tadpoles at equal masses\n";
  for(int i=0;i<4;i++){
  for(int j=0;j<4;j++){
  for(int k=0;k<4;k++){
    int ijk = intstring(lmu[i]+lmu[j]+lmu[k]);
    fourvector qzero;
    d1 =  BbVtwistt(1,33,ijk,m1sq,m1sq,0.,xl,theta1,qzero);
    d2 = AbVtwistt(2,33,ijk,m1sq,xl,theta1);
    cout << lmu[i]+lmu[j]+lmu[k]
	 <<' '<<d1<<' '<<d2<<' '<<d1/d2-1.
	 <<'\n';
  }}}

  cout <<"B33bV reduces at qsq=0 to the correct tadpoles at different masses\n";
  for(int i=0;i<4;i++){
  for(int j=0;j<4;j++){
  for(int k=0;k<4;k++){
    int ijk = intstring(lmu[i]+lmu[j]+lmu[k]);
    fourvector qzero;
    d1 =  BbVtwistt(1,33,ijk,m1sq,m2sq,0.,xl,theta1,qzero);
    // note the + sign since I need to do r -> -r in the second term
    d2 =( AbVtwistt(1,33,ijk,m1sq,xl,theta1)
	 +AbVtwistt(1,33,ijk,m2sq,xl,-theta1))/(m1sq-m2sq);
    cout << lmu[i]+lmu[j]+lmu[k]
	 <<' '<<d1<<' '<<d2<<' '<<d1/d2-1.
	 <<'\n';
  }}}

  cout < "Checking Passarino-Veltman like g_munu relations\n";
  for(int i=0;i<4;i++){
    int i00 = intstring(lmu[0]+lmu[0]+lmu[i]);
    int i11 = intstring(lmu[1]+lmu[1]+lmu[i]);
    int i22 = intstring(lmu[2]+lmu[2]+lmu[i]);
    int i33 = intstring(lmu[3]+lmu[3]+lmu[i]);

    d1 =   qsq*qex[i]*BbVtwistt(1,31,  0,m1sq,m2sq,qsq,xl,theta1,qext)
           +6.*qex[i]*BbVtwistt(1,32,  0,m1sq,m2sq,qsq,xl,theta1,qext)
                     +BbVtwistt(1,33,i00,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(1,33,i11,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(1,33,i22,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(1,33,i33,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 =  m1sq*qex[i]*BbVtwistt(1,1,0,m1sq,m2sq,qsq,xl,theta1,qext)
                +m1sq*BbVtwistt(1,2,i,m1sq,m2sq,qsq,xl,theta1,qext)
        +qex[i]*AbVtwistt(1,0,0,m2sq,xl,theta2)
               -AbVtwistt(1,2,i,m2sq,xl,theta2);
    cout << i<<' '<<d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';

    d1 =   qsq*qex[i]*BbVtwistt(2,31,  0,m1sq,m2sq,qsq,xl,theta1,qext)
           +6.*qex[i]*BbVtwistt(2,32,  0,m1sq,m2sq,qsq,xl,theta1,qext)
                     +BbVtwistt(2,33,i00,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(2,33,i11,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(2,33,i22,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(2,33,i33,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 =  m1sq*qex[i]*BbVtwistt(2,1,0,m1sq,m2sq,qsq,xl,theta1,qext)
                +m1sq*BbVtwistt(2,2,i,m1sq,m2sq,qsq,xl,theta1,qext)
         +qex[i]*BbVtwistt(1,1,0,m1sq,m2sq,qsq,xl,theta1,qext)
                +BbVtwistt(1,2,i,m1sq,m2sq,qsq,xl,theta1,qext);
    cout << i<<' '<<d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';


    d1 =   qsq*qex[i]*BbVtwistt(3,31,  0,m1sq,m2sq,qsq,xl,theta1,qext)
           +6.*qex[i]*BbVtwistt(3,32,  0,m1sq,m2sq,qsq,xl,theta1,qext)
                     +BbVtwistt(3,33,i00,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(3,33,i11,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(3,33,i22,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(3,33,i33,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 =  m1sq*qex[i]*BbVtwistt(3,1,0,m1sq,m2sq,qsq,xl,theta1,qext)
                +m1sq*BbVtwistt(3,2,i,m1sq,m2sq,qsq,xl,theta1,qext)
        +qex[i]*AbVtwistt(2,0,0,m2sq,xl,theta2)
               -AbVtwistt(2,2,i,m2sq,xl,theta2);
    cout << i<<' '<<d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';

    d1 =   qsq*qex[i]*BbVtwistt(4,31,  0,m1sq,m2sq,qsq,xl,theta1,qext)
           +6.*qex[i]*BbVtwistt(4,32,  0,m1sq,m2sq,qsq,xl,theta1,qext)
                     +BbVtwistt(4,33,i00,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(4,33,i11,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(4,33,i22,m1sq,m2sq,qsq,xl,theta1,qext)
                     -BbVtwistt(4,33,i33,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 =  m1sq*qex[i]*BbVtwistt(4,1,0,m1sq,m2sq,qsq,xl,theta1,qext)
                +m1sq*BbVtwistt(4,2,i,m1sq,m2sq,qsq,xl,theta1,qext)
         +qex[i]*BbVtwistt(3,1,0,m1sq,m2sq,qsq,xl,theta1,qext)
                +BbVtwistt(3,2,i,m1sq,m2sq,qsq,xl,theta1,qext);
    cout << i<<' '<<d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  }


  cout << "Checking the Passarino-Veltman-like  q_mu relation\n";
  // define g_{munu}
  int gmunu[4][4] ={{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
  for(int i=0;i<4;i++){
  for(int j=0;j<4;j++){
    int ij = intstring(lmu[i]+lmu[j]);
    int ij0 = intstring(lmu[0]+lmu[i]+lmu[j]);
    int ij1 = intstring(lmu[1]+lmu[i]+lmu[j]);
    int ij2 = intstring(lmu[2]+lmu[i]+lmu[j]);
    int ij3 = intstring(lmu[3]+lmu[i]+lmu[j]);

    d1 =   qsq*qex[i]*qex[j]*BbVtwistt(1,31,  0,m1sq,m2sq,qsq,xl,theta1,qext)
      +(gmunu[i][j]*qsq+2.*qex[i]*qex[j])*
                             BbVtwistt(1,32,  0,m1sq,m2sq,qsq,xl,theta1,qext)
      +qex[0]*BbVtwistt(1,33,ij0,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[1]*BbVtwistt(1,33,ij1,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[2]*BbVtwistt(1,33,ij2,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[3]*BbVtwistt(1,33,ij3,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 = 0.5*(m1sq-m2sq+qsq)*(
                qex[i]*qex[j]*BbVtwistt(1,21,0,m1sq,m2sq,qsq,xl,theta1,qext)
                 +gmunu[i][j]*BbVtwistt(1,22,0,m1sq,m2sq,qsq,xl,theta1,qext)
		            +BbVtwistt(1,23,ij,m1sq,m2sq,qsq,xl,theta1,qext))
      +0.5*gmunu[i][j]*AbVtwistt(1,22,0,m2sq,xl,theta2)
      +0.5*AbVtwistt(1,23,ij,m2sq,xl,theta2)
      -0.5*qex[i]*AbVtwistt(1,2,j,m2sq,xl,theta2)
      -0.5*qex[j]*AbVtwistt(1,2,i,m2sq,xl,theta2)
      +0.5*qex[i]*qex[j]*AbVtwistt(1,0,0,m2sq,xl,theta2)
      -0.5*gmunu[i][j]*AbVtwistt(1,22,0,m1sq,xl,theta1)
      -0.5*AbVtwistt(1,23,ij,m1sq,xl,theta1);
    cout << i<<j<<' '<<d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';

    d1 =   qsq*qex[i]*qex[j]*BbVtwistt(2,31,  0,m1sq,m2sq,qsq,xl,theta1,qext)
      +(gmunu[i][j]*qsq+2.*qex[i]*qex[j])*
                             BbVtwistt(2,32,  0,m1sq,m2sq,qsq,xl,theta1,qext)
      +qex[0]*BbVtwistt(2,33,ij0,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[1]*BbVtwistt(2,33,ij1,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[2]*BbVtwistt(2,33,ij2,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[3]*BbVtwistt(2,33,ij3,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 = 0.5*(m1sq-m2sq+qsq)*(
                qex[i]*qex[j]*BbVtwistt(2,21,0,m1sq,m2sq,qsq,xl,theta1,qext)
                 +gmunu[i][j]*BbVtwistt(2,22,0,m1sq,m2sq,qsq,xl,theta1,qext)
		            +BbVtwistt(2,23,ij,m1sq,m2sq,qsq,xl,theta1,qext))
        +0.5*qex[i]*qex[j]*BbVtwistt(1,21, 0,m1sq,m2sq,qsq,xl,theta1,qext)
          +0.5*gmunu[i][j]*BbVtwistt(1,22, 0,m1sq,m2sq,qsq,xl,theta1,qext)
                      +0.5*BbVtwistt(1,23,ij,m1sq,m2sq,qsq,xl,theta1,qext)
      -0.5*gmunu[i][j]*AbVtwistt(2,22,0,m1sq,xl,theta1)
      -0.5*AbVtwistt(2,23,ij,m1sq,xl,theta1);
    cout << i<<j<<' '<<d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';

    d1 =   qsq*qex[i]*qex[j]*BbVtwistt(3,31,  0,m1sq,m2sq,qsq,xl,theta1,qext)
      +(gmunu[i][j]*qsq+2.*qex[i]*qex[j])*
                             BbVtwistt(3,32,  0,m1sq,m2sq,qsq,xl,theta1,qext)
      +qex[0]*BbVtwistt(3,33,ij0,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[1]*BbVtwistt(3,33,ij1,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[2]*BbVtwistt(3,33,ij2,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[3]*BbVtwistt(3,33,ij3,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 = 0.5*(m1sq-m2sq+qsq)*(
                qex[i]*qex[j]*BbVtwistt(3,21,0,m1sq,m2sq,qsq,xl,theta1,qext)
                 +gmunu[i][j]*BbVtwistt(3,22,0,m1sq,m2sq,qsq,xl,theta1,qext)
		            +BbVtwistt(3,23,ij,m1sq,m2sq,qsq,xl,theta1,qext))
      +0.5*gmunu[i][j]*AbVtwistt(2,22,0,m2sq,xl,theta2)
      +0.5*AbVtwistt(2,23,ij,m2sq,xl,theta2)
      -0.5*qex[i]*AbVtwistt(2,2,j,m2sq,xl,theta2)
      -0.5*qex[j]*AbVtwistt(2,2,i,m2sq,xl,theta2)
      +0.5*qex[i]*qex[j]*AbVtwistt(2,0,0,m2sq,xl,theta2)
        -0.5*qex[i]*qex[j]*BbVtwistt(1,21, 0,m1sq,m2sq,qsq,xl,theta1,qext)
          -0.5*gmunu[i][j]*BbVtwistt(1,22, 0,m1sq,m2sq,qsq,xl,theta1,qext)
                      -0.5*BbVtwistt(1,23,ij,m1sq,m2sq,qsq,xl,theta1,qext);
    cout << i<<j<<' '<<d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';

    d1 =   qsq*qex[i]*qex[j]*BbVtwistt(4,31,  0,m1sq,m2sq,qsq,xl,theta1,qext)
      +(gmunu[i][j]*qsq+2.*qex[i]*qex[j])*
                             BbVtwistt(4,32,  0,m1sq,m2sq,qsq,xl,theta1,qext)
      +qex[0]*BbVtwistt(4,33,ij0,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[1]*BbVtwistt(4,33,ij1,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[2]*BbVtwistt(4,33,ij2,m1sq,m2sq,qsq,xl,theta1,qext)
      -qex[3]*BbVtwistt(4,33,ij3,m1sq,m2sq,qsq,xl,theta1,qext);
    d2 = 0.5*(m1sq-m2sq+qsq)*(
                qex[i]*qex[j]*BbVtwistt(4,21,0,m1sq,m2sq,qsq,xl,theta1,qext)
                 +gmunu[i][j]*BbVtwistt(4,22,0,m1sq,m2sq,qsq,xl,theta1,qext)
		            +BbVtwistt(4,23,ij,m1sq,m2sq,qsq,xl,theta1,qext))
        +0.5*qex[i]*qex[j]*BbVtwistt(3,21, 0,m1sq,m2sq,qsq,xl,theta1,qext)
          +0.5*gmunu[i][j]*BbVtwistt(3,22, 0,m1sq,m2sq,qsq,xl,theta1,qext)
                      +0.5*BbVtwistt(3,23,ij,m1sq,m2sq,qsq,xl,theta1,qext)
        -0.5*qex[i]*qex[j]*BbVtwistt(2,21, 0,m1sq,m2sq,qsq,xl,theta1,qext)
          -0.5*gmunu[i][j]*BbVtwistt(2,22, 0,m1sq,m2sq,qsq,xl,theta1,qext)
                      -0.5*BbVtwistt(2,23,ij,m1sq,m2sq,qsq,xl,theta1,qext);
    cout << i<<j<<' '<<d1<<' '<<d2<<' '<<d1/d2-1.<<'\n';
  }}

  return 0;
}
