#include<iostream>
#include<Eigen/Dense>
#include<ctime>
#include<vector>
#include<complex>
#include<fstream>
#include "lib/dft.cpp"


std::ofstream ofile;

int main(int argc, char* argv[]){
  double time,x,dx,L,omega;
  std::complex<double> z;
  int nx,no;
  std::clock_t start, end;
  
  //Assign discretization
  Eigen::MatrixXcd D,Dstar;
  Eigen::VectorXcd F,f,ftest;  
  nx = 500;
  no = 200;
  L = 1000;
  dx = L/(double) nx;

  //Construct f(x), D, and F(omega)
  project_function(f, nx, L);
  DFT dft(nx,no);
  dft.set_limit(0,1000);
  F = dft.D*f;
  ftest = dft.Dstar*F/(double) nx;

  //print to file
  start = clock();
  write(f,ftest,F,L);
  end = clock();
  time = (end - start)/(double) CLOCKS_PER_SEC;
  printf("time: %f\n",time);  


  return 0;
}
