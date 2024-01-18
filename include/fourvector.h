// fourvector.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of the class fourvector as well
// as a few associated functions


#ifndef FOURVECTOR_H
#define FOURVECTOR_H

#include <iostream>
#include <vector>
#include <iomanip>

class fourvector{
private:
  std::vector<double> v={0.,0.,0.,0.};
  //  std::vector<double> v (4,0.);// initializes four zero elements
public:
  fourvector(const double a0=0.,const double a1=0.,const double a2=0.,
	     const double a3=0.);
  fourvector(const std::vector<double> vv);
  fourvector operator-(const fourvector vv) const;
  fourvector operator+(const fourvector vv) const;
  fourvector& operator-=(const fourvector vv);
  fourvector& operator+=(const fourvector vv);
  fourvector operator-(void) const;
  fourvector operator*(const double aa) const;
  double operator*(const fourvector vv) const;
  fourvector operator/(const double aa) const;
  friend fourvector operator*(const double aa, const fourvector vv);
  // allows a fourvector to be used as a vector when accessing values
  operator std::vector<double>() const;
  double out(const int n) const;
  void set(const int, const double xn);
  void set(const double xn, const int n);
  void lower(void);

  friend std::ostream & operator<<(std::ostream & os,const fourvector & vv);
  friend std::istream & operator>>(std::istream & is, fourvector & vv);

  ~fourvector(void){}
};

#endif // FOURVECTOR_H
