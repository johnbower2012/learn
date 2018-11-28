#ifndef __EMULATOR_H__
#define __EMULATOR_H__

#include<Eigen/Dense>
#include<chrono>
#include<random>
#include<cmath>
#include<iostream>

class emulator{
public:
  int trainPoints;
  int paramCount;
  int obsCount;
  int hyperparamCount;
  
  double epsilon=1e-8;
  Eigen::VectorXd noise;

  Eigen::MatrixXd X;
  Eigen::MatrixXd Hyperparam;
  Eigen::MatrixXd Beta;
  Eigen::MatrixXd Noise;

  std::vector<Eigen::MatrixXd> Kernel;
  std::vector<Eigen::MatrixXd> KernelInv;
  std::vector<Eigen::MatrixXd> HMatrix;

  emulator(Eigen::MatrixXd newX, Eigen::MatrixXd newHyperparam, Eigen::MatrixXd newBeta);

  void kernelFunction(Eigen::MatrixXd A, Eigen::MatrixXd B, int obsIndex, Eigen::MatrixXd &kernlMatrix);
  void regressionLinearFunction(Eigen::MatrixXd testX, int obsIndex, Eigen::MatrixXd &HMatrix);
  void Emulate(Eigen::MatrixXd testX, Eigen::MatrixXd Y, Eigen::MatrixXd &outMatrix);
};

#endif
