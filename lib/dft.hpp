#ifndef DFT_H
#define DFT_H

#include<iostream>
#include<Eigen/Dense>
#include<ctime>
#include<vector>
#include<complex>
#include<fstream>

std::complex<double> f_function(double x);
void project_function(Eigen::VectorXcd &f, int &nx, double &L);
void write(Eigen::VectorXcd &f, Eigen::VectorXcd &ftest, Eigen::VectorXcd &F, double &L);

class DFT{
public:
  int nx, no;
  std::vector<double> L;
  Eigen::MatrixXcd D,Dstar;

  DFT();
  DFT(int N);
  DFT(int NX, int NO);

  void set_limit(double start, double end);
};

#endif
