// testvectorformlo.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// testing the vector form-factors as expressed in lowest order masses

#include<iostream>
#include<complex>

#include "vectorformlo.h"

#include "oneloopintegrals.h"
int main(void){
  using namespace std;
  // electromagnetic form-factor
  // first some more physical values to compare with the results of
  // Bijnens-Talavera 2002, Fig 4 and Fig 5 there
  lomass mass1(0.13956995,0.493677,0.0924,0.77);
  Li Li1;
  Li1.setmu(0.77);
  Li1.setli(9,0.0069);
  cout << "Inputs used for comparing Bijnens-Talavera 2002\n";
  cout << mass1;
  cout << Li1;
  cout << "fvpipp4Rlo : " << fvpipp4Rlo(-0.25,mass1)
       <<' '<< fvpipp4Rlo(0.,mass1)
       <<' '<< fvpipp4Rlo(0.25,mass1)<<'\n';
  cout << "fvpipp4Llo : " << fvpipp4Llo(-0.25,mass1,Li1)
       <<' '<< fvpipp4Llo(0.,mass1,Li1)
       <<' '<< fvpipp4Llo(0.25,mass1,Li1)<<'\n';
  cout << "fvkpp4Rlo  : " << fvkpp4Rlo(-0.25,mass1)
       <<' '<< fvkpp4Rlo(0.,mass1)
       <<' '<< fvkpp4Rlo(0.25,mass1)<<'\n';
  cout << "fvkpp4Llo  : " << fvkpp4Llo(-0.25,mass1,Li1)
       <<' '<< fvkpp4Llo(0.,mass1,Li1)
       <<' '<< fvkpp4Llo(0.25,mass1,Li1)<<'\n';
  cout << "fvkop4Rlo  : " << fvkop4Rlo(-0.25,mass1)
       <<' '<< fvkop4Rlo(0.,mass1)
       <<' '<< fvkop4Rlo(0.25,mass1)<<'\n';
  cout << "fvkop4Llo  : " << fvkop4Llo(-0.25,mass1,Li1)
       <<' '<< fvkop4Llo(0.,mass1,Li1)
       <<' '<< fvkop4Llo(0.25,mass1,Li1)<<'\n';

  // Kto pi or Kell3 form-factors
  // compare with Bijnens-Talavera Fig. 6 and 7
  cout << "fvpkpip4Rlo: " << fvpkpip4Rlo(-0.25,mass1)
       <<' '<< fvpkpip4Rlo(0.,mass1)
       <<' '<< fvpkpip4Rlo(0.25,mass1)<<'\n';
  cout << "fvmkpip4Rlo: " << fvmkpip4Rlo(-0.25,mass1)
       <<' '<< fvmkpip4Rlo(0.,mass1)
       <<' '<< fvmkpip4Rlo(0.25,mass1)<<'\n';
  cout << "fvmkpip4Rlo (BT): " << fvmkpip4RloBT(-0.25,mass1)
       <<' '<< fvmkpip4RloBT(0.,mass1)
       <<' '<< fvmkpip4RloBT(0.25,mass1)<<'\n';

  // for equal masses all formfactors should be the same and f_- vanish
  lomass mass2(0.21,0.21,0.0924,0.77);
  cout << "All f_+ formfactors equal and f_- vanish for equal masses:\n" ;

  cout << "fvpipp4Rlo : " << fvpipp4Rlo(-0.25,mass2)
       <<' '<< fvpipp4Rlo(0.,mass2)
       <<' '<< fvpipp4Rlo(0.25,mass2)<<'\n';
  cout << "fvkpp4Rlo  : " << fvkpp4Rlo(-0.25,mass2)
       <<' '<< fvkpp4Rlo(0.,mass2)
       <<' '<< fvkpp4Rlo(0.25,mass2)<<'\n';
  cout << "fvpkpip4Rlo: " << fvpkpip4Rlo(-0.25,mass2)
       <<' '<< fvpkpip4Rlo(0.,mass2)
       <<' '<< fvpkpip4Rlo(0.25,mass2)<<'\n';
  cout << "fvmkpip4Rlo: " << fvmkpip4Rlo(-0.25,mass2)
       <<' '<< fvmkpip4Rlo(0.,mass2)
       <<' '<< fvmkpip4Rlo(0.25,mass2)<<'\n';
 return 0;
}
