// li.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2017 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.
#ifndef LINF2_H
#define LINF2_H

class libar;

class li{
protected:
double l1r,l2r,l3r,l4r,l5r,l6r,l7r,h1r,h2r,h3r;
  double mu;
  std::string name;
public:
  li(const double l1r = 0.,const double l2r = 0.,const double l3r = 0.,
     const double l4r = 0.,const double l5r = 0.,const double l6r = 0.,
     const double l7r = 0.,const double h1r = 0.,const double h2r = 0.,
     const double h3r = 0.,
     const double mu = 0.77,
     const std::string Name = "Nameless li");
  li(const libar libarin, const double mu = 0.77, const double mpi = 0.13957061);
  ~li(void);
  void setli(const int n, const double lin);
  void setli(const double lin, const int n);
  void setmu(const double muin);
  void setname(const std::string inputname);
friend std::ostream & operator<<(std::ostream & os,const li & bb); // does the output
friend std::istream & operator>>(std::istream & is, li & liout);// reads output in again
  li operator+(const li & bb) const; // defines the sum of two sets of li
  li operator-(const li & bb) const; // defines the difference of two sets of li
  li operator*(const double & xx) const; // multiplies li with a number
  friend li operator*(const double & aa,const  li & bb);// aa*li
  void changescale(const double newmu);
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,
	     double & H1t,double & H2t,double & H3t,
	     double & mut,std::string nameout) const;
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,
	     double & H1t,double & H2t,double &H3t,
	     double & mut) const;
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & H1t,double & H2t,double & H3t) const;
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,
	     double & mut) const;
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t) const;
  double out(const int n) const;
  double getmu(void) const;
};

class libar{
protected:
  double l1bar,l2bar,l3bar,l4bar,l5bar,l6bar,l7r;
  double h1bar,h2bar,h3r;
  std::string name;
public:
  libar(const double l1bar = 0.,const double l2bar = 0.,const double l3bar = 0.,
     const double l4bar = 0.,const double l5bar = 0.,const double l6bar = 0.,
     const double l7r = 0.,const double h1bar = 0.,const double h2bar = 0.,
     const double h3r = 0.,
     const std::string Name = "Nameless libar");
  libar(const li liin, const double mpi = 0.13957061);
  ~libar(void);
  void setlibar(const int n, const double lin);
  void setlibar(const double lin, const int n);
  void setname(const std::string inputname);
  void print(void) const;
friend std::ostream & operator<<(std::ostream & os,const libar & bb); // does the output
friend std::istream & operator>>(std::istream & is, libar & liout);//reads output in again
  void out(double & l1t,double & l2t,double & l3t,
		double & l4t,double & l5t,double & l6t,
		double & l7t,double & h1t,double & h2t,
		  double & h3t, std::string nameout) const;
  void out(double & l1t,double & l2t,double & l3t,
		double & l4t,double & l5t,double & l6t,
		double & l7t,double & h1t,double & h2t,
	   double & h3t) const;
  void out(double & l1t,double & l2t,double & l3t,
		double & l4t,double & l5t,double & l6t,
		double & l7t) const;
  double out(const int n) const;
};
//
//libar libarfromLip4(const physmass massin, const Li Liin);
//libar libarfromLip6(const physmass massin, const Li Liin, const Ci Ciin);


#endif
