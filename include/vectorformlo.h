// vectorformlo.h is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#ifndef VECTORFORMLO_H
#define VECTORFORMLO_H


#include<complex>
#include "Li.h"
#include "inputs.h"

typedef std::complex<double> dcomplex;

// electromagnetic form-factors
dcomplex fvpipp4lo(const double t, const lomass mass,const Li Liin);
dcomplex fvpipp4Rlo(const double t, const lomass mass);
double fvpipp4Llo(const double t, const lomass mass, const Li Liiin);
dcomplex fvkpp4lo(const double t, const lomass mass,const Li Liin);
dcomplex fvkpp4Rlo(const double t, const lomass mass);
double fvkpp4Llo(const double t, const lomass mass, const Li Liiin);
dcomplex fvkop4lo(const double t, const lomass mass,const Li Liin);
dcomplex fvkop4Rlo(const double t, const lomass mass);
double fvkop4Llo(const double t, const lomass mass, const Li Liiin);
// K to pi or Kl3 form-factors
dcomplex fvpkpip4lo(const double t, const lomass mass,const Li Liin);
dcomplex fvpkpip4Rlo(const double t, const lomass mass);
double fvpkpip4Llo(const double t, const lomass mass, const Li Liiin);
dcomplex fvmkpip4lo(const double t, const lomass mass,const Li Liin);
dcomplex fvmkpip4Rlo(const double t, const lomass mass);
dcomplex fvmkpip4RloBT(const double t, const lomass mass);
double fvmkpip4Llo(const double t, const lomass mass, const Li Liiin);


#endif // VECTORFORMLO_H
