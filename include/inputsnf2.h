// inputsnf2.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2019 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of classes physmassnf2, quarkmassn2,
// lomassnf2 as well as associated functions


#ifndef INPUTSNF2_H
#define INPUTSNF2_H

#include <string>
#include <iostream>
using namespace std;

class physmassnf2{
 private:
  double mpi,fpi,mu;
 public:
  physmassnf2(const double mpiin =0.135, const double fpiin=0.0922,
	      const double muin =0.77);
  void setmpi(const double mpiin=0.135);
  void setfpi(const double fpiin=0.0922);
  void setmu(const double muiin=0.77);
  void out(double &mpiout, double &fpiout, double &muout) const;
  double getmpi(void) const;
  double getfpi(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const physmassnf2 & bb); // output
  friend istream & operator>>(istream & is,physmassnf2 & bb); // input
  friend bool operator==(const physmassnf2 mass1,const physmassnf2 mass2);
  ~physmassnf2(void);
};

class quarkmassnf2;
class lomassnf2{
 private:
  double mp,f,mu;
 public:
  lomassnf2(const double mp =0.135, const double f=0.090,
	    const double muin =0.77);
  lomassnf2(const quarkmassnf2 mass);
  void setmp(const double mpin=0.135);
  void setf(const double fin=0.09);
  void setmu(const double muin=0.77);
  void out(double &mpout, double &fout, double &muout) const;
  double getmp(void) const;
  double getf(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const lomassnf2 & bb); // output
  friend istream & operator>>(istream & is,lomassnf2 & bb); // input
  friend bool operator==(const lomassnf2 mass1,const lomassnf2 mass2);
  ~lomassnf2(void);
};

class quarkmassnf2{
 private:
  double Bmhat,f,mu;
 public:
  quarkmassnf2(const double Bmhat =0.01, const double f0=0.090,
	       const double muin =0.77);
  quarkmassnf2(const lomassnf2 mass);
  void setBmhat(const double Bmhatin=0.01);
  void setf(const double fin=0.09);
  void setmu(const double muin=0.77);
  void out(double &B0mhatout, double &f0out, double &muout) const;
  double getBmhat(void) const;
  double getf(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const quarkmassnf2 & bb); // output
  friend istream & operator>>(istream & is,quarkmassnf2 & bb); // input
  friend bool operator==(const quarkmassnf2 mass1,const quarkmassnf2 mass2);
  ~quarkmassnf2(void);
};
#endif // INPUTSNF2_H
