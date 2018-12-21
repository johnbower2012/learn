#ifndef __MCMC_H__
#define __MCMC_H__

#include<Eigen/Dense>
#include<random>
#include<chrono>
#include "emulator.h"

class MCMC{
 public:
  double maxLogLikelihood;
  double Likelihood;
  int obsCount;
  int paramCount;
  int Samples;
  Eigen::MatrixXd targetValue;
  Eigen::MatrixXd Range;

  Eigen::MatrixXd Position;
  Eigen::MatrixXd testPosition;
  Eigen::VectorXd Widths;

  unsigned seed;
  std::default_random_engine generator;
  std::normal_distribution<double> normal_dist;
  std::uniform_real_distribution<double> uniform_dist;

  bool NR;
  MCMC();
  MCMC(Eigen::MatrixXd newTargetValue, Eigen::MatrixXd newRange, Eigen::VectorXd newWidths, bool NR);

  void setTargetValue(Eigen::MatrixXd newTargetValue);
  void setRange(Eigen::MatrixXd newRange);
  void setPosition();
  void setPosition(Eigen::MatrixXd newPosition);
  void setWidths(Eigen::VectorXd newWidths);

  void setSeed(unsigned seed);
  void setSeedClock();

  double normal();
  double uniform();

  void step();
  double getLogLikelihood();
  double getLogLikelihood(Eigen::MatrixXd Z);
  bool decide(Eigen::MatrixXd Z);

  void Run(int Samples, Eigen::MatrixXd &History, emulator obsEmulator, Eigen::MatrixXd Y);
};

#endif
