#ifndef COSHFUNC_H
#define COSHFUNC_H

#include<iostream>
#include<cmath>
#include<fstream>
#include<Eigen/Dense>
#include<random>
#include <cstdlib>
#include <cstdio>
#include <complex>


#include "parametermap.h"
#define PI 4.0*atan(1.0)

class CCosh{
public:
  int nmax;
  Eigen::MatrixXd A;
  Eigen::VectorXd I;
  Eigen::VectorXd Z;

  CCosh();
  CCosh(int nmaxset);

  void Set_NMax(int nmaxset);
  double GetG(int n,double x);
  double GetRandomX();
  void TestOverlap();
  void CalcOverlap(int n,int nprime);
  void CalcOverlapNumerical(int n,int nprime);
  void Integrate(double x);
  
  int seed;
  std::mt19937 rng;
  std::uniform_real_distribution<double> uniform_dist;
};

class CDistribution{
 public:
  CDistribution();
  string parsfilename;
  CParameterMap parmap;
  virtual void GenX(int ab,double &x,double &weight);
};

class CDistCosh : public CDistribution{
 public:
  std::vector<double> width;
  std::vector<std::vector<double> > g;
  CCosh coshcalc;  

  //  CDistCosh(string parsfilename);
 CDistCosh(std::vector<double> Width, std::vector<std::vector<double> > G);

  void Set_WG(std::vector<double> Width, std::vector<std::vector<double> > G);
  void GenX(int ab,double &x,double &weight);
  void Functions(int n, double x, Eigen::MatrixXd &M);
  double Integrate(int ab);
  double IntegrateNumerical(int ab);
};


#endif
