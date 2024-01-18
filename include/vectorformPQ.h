// vectorformPQ.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2016 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include<complex>
#include "inputsnf.h"
#include "Linf.h"
typedef std::complex<double> dcomplex;

double fvpnf3p4L(const double qsq, const quarkmassnf mass,const Linf Liin);

dcomplex fvpv2s2nf3p4R(const double qsq, const quarkmassnf mass);
dcomplex fvpv2s2nf3p4Rp(const double qsq, const quarkmassnf mass);
dcomplex fvpv2s3nf3p4R(const double qsq, const quarkmassnf mass);
dcomplex fvpv2s3nf3p4Rp(const double qsq, const quarkmassnf mass);
dcomplex fvpv3s2nf3p4R(const double qsq, const quarkmassnf mass);
dcomplex fvpv3s3nf3p4R(const double qsq, const quarkmassnf mass);
