#include<iostream>
#include<Eigen/Dense>
#include "emulator.h"
#include "system.h"
#include "mcmc.h"

int main(int argc, char* argv[]){
  int train=100,
    test=200,
    param=2,
    obs=1,
    hpa=3;
  double xinit=0.0,
    xfin=10,
    dx;
  dx = (xfin-xinit)/(double) train;
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(train,param);
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(train,obs);
  Eigen::MatrixXd hp = Eigen::MatrixXd::Zero(obs,hpa);
  Eigen::MatrixXd in = Eigen::MatrixXd::Zero(train,param+obs);
  Eigen::MatrixXd output;
  for(int i=0;i<train;i++){
    X(i,0) = xinit + dx*i;
    X(i,1) = xinit + dx*((i*i*i)%train);
    Y(i,0) = (X(i,0)-5)*(X(i,0)-5) + fabs(X(i,1) - 1);
    in(i,0) = X(i,0);
    in(i,1) = X(i,1);
    in(i,2) = Y(i,0);
  }

  hp(0,0) = 2.0;
  hp(0,1) = 2.0;
  hp(0,2) = 0.01;

  emulator emul(X,hp);

  Eigen::MatrixXd Range(param,2);
  Range(0,0) =Range(1,0)= xinit;
  Range(0,1) =Range(1,1)= xfin;
  Eigen::MatrixXd targetValue(1,1);
  targetValue(0,0) = 0;
  Eigen::VectorXd Widths(param);
  Widths(0) = (xfin - xinit)/10.0;
  Widths(1) = (xfin - xinit)/10.0;
  int Samples=100000;

  Eigen::MatrixXd History;
  Eigen::MatrixXd position(1,param);
  position(0,0) = 5;
   position(0,1) = 5;
  MCMC mcmc(targetValue,Range,Widths);
  mcmc.setPosition(position);
  std::vector<std::string> header(param);
  header[0] = "x1";
  header[1] = "x2";
  mcmc.Run(Samples,History,emul,Y);
  Eigen::MatrixXd HisTory = History.block(0,0,Samples,param);  
  WriteCSVFile("mcmctrace.csv",header,HisTory,",");
  WriteFile("in.dat",in," ");

  return 0;
}
