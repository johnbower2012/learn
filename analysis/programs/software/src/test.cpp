#include<iostream>
#include<Eigen/Dense>
#include "emulator.h"
#include "system.h"

int main(int argc, char* argv[]){
  int train=10,
    test=200,
    param=1,
    obs=1,
    hpa=3;
  double xinit=0.0,
    xfin=6.28,
    dx;
  dx = (xfin-xinit)/(double) train;
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(train,param);
  Eigen::MatrixXd XS = Eigen::MatrixXd::Zero(test,param);  
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(train,obs);
  Eigen::MatrixXd hp = Eigen::MatrixXd::Zero(obs,hpa);
  Eigen::MatrixXd beta = Eigen::MatrixXd::Zero(obs,1+param);
  Eigen::MatrixXd in = Eigen::MatrixXd::Zero(train,param+obs);
  Eigen::MatrixXd output;
  for(int i=0;i<train;i++){
    X(i,0) = xinit + dx*i;
    Y(i,0) = cos(X(i,0));
    in(i,0) = X(i,0);
    in(i,1) = Y(i,0);
  }
  dx = (xfin-xinit)/(double) (test);
  for(int i=0;i<test;i++){
    XS(i,0) = xinit + dx*i;
  }
  beta(0,0) = 0.0;
  beta(0,1) = 0.0001;
  hp(0,0) = 1.0;
  hp(0,1) = 0.5;
  hp(0,2) = 0.01;
  
  emulator emul(X,hp,beta);
  emul.Emulate(XS,Y,output);
  std::cout << in << std::endl;
  std::cout << output << std::endl;
  WriteFile("in.dat",in," ");
  WriteFile("out.dat",output," ");

  return 0;
}
