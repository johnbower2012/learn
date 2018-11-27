#include "coshfunc.cpp"
#include<iostream>
#include<fstream>
#include<Eigen/Dense>

std::ofstream ofile;



int main(int argc, char* argv[]){
  int ab=4;
  int nmax=10;
  std::vector<double> width(ab);
  std::vector<std::vector<double> > g(ab,std::vector<double>(nmax+1));
  CCosh cosh(nmax);

  for(int i=0;i<ab;i++){
    width[i] = i+1;
    for(int j=0;j<nmax+1;j++){
      g[i][j] = j*pow(-1,j)+2;
    }
    g[i][0] = 1.0;
    for(int j=1;j<nmax+1;j++){
      g[i][0] -= g[i][j]*cosh.Z(j);
    }
    g[i][0] /= cosh.Z(0);
  }
    
  CDistCosh dist(width,g);
  double 
    x=0.0,
    sumx=0.0,
    sumw=0.0,
    weight=0.0;
  int 
    n=100000;

  n=1000;
  x=50;
  Eigen::MatrixXd M;
  dist.Functions(n,x,M);
  ofile.open("wave.dat");
  for(int j=0;j<n;j++){
    for(int i=0;i<ab+1;i++){
      ofile << M(j,i) << " ";
    }
    ofile << std::endl;
  }
  ofile.close();
    
  
  /*
  for(int i=0;i<n;i++){
    dist.GenX(0,x,weight);
    sumx += x;
    sumw += weight;
  }
  std::cout << sumx/double(n) << std::endl;
  std::cout << sumw/double(n) << std::endl;
  */
  return 0;
}
