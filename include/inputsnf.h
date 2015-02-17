// inputsnf.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of classes physmass, Li and Ci as well
// as associated functions


#ifndef INPUTSNF_H
#define INPUTSNF_H

#include<string>
#include<vector>
#include<iostream>
using namespace std;

class quarkmassnf{
 private:
  int nq;// holds the number of quark species
  vector<double> B0mq;// length should be nq if initialized correctly
  double f0,mu;
 public:
  quarkmassnf(const double f0=0.090, const double muin =0.77, const int nqin=3);
  quarkmassnf(const vector<double> B0mqin,const double f0in=0.090,
	      const double muin =0.77);
  void setf0(const double f0in=0.09);
  void setmu(const double muin=0.77);
  void setB0mq(const vector<double> B0mqin);
  void setB0mq(const double B0mqin, const int i);
  void setB0mq(const int i, const double B0mqin=0.);
  void out(vector<double> &B0mq) const;
  void out(vector<double> &B0mq, double &f0out, double &muout) const;
  void out(vector<double> &B0mq, double &f0out, double &muout, int &nq) const;
  int getnq(void) const;
  vector<double> getB0mq(void) const;
  double getB0mq(const int i) const;
  double getf0(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const quarkmassnf & bb); // output
  friend istream & operator>>(istream & is,quarkmassnf & bb); // input
  friend bool operator==(const quarkmassnf mass1,const quarkmassnf mass2);
  ~quarkmassnf(void);
};
#endif // INPUTSNF_H
