#ifndef __EMULATOR_CLASS_H__
#define __EMULATOR_CLASS_H__

#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include<Eigen/Dense>
#include "time.h"

class emulator{
public:
  int train_points;
  int parameter_count;
  int observable_count;
  int hyperparameter_count;

  double epsilon;
  Eigen::VectorXd noise;

  Eigen::MatrixXd X;
  std::vector<Eigen::MatrixXd> kernel;
  std::vector<Eigen::MatrixXd> kernel_inverse;
  Eigen::MatrixXd hyperparameters;
  std::vector<Eigen::MatrixXd> H;
  Eigen::MatrixXd beta;

  emulator(Eigen::MatrixXd X, Eigen::MatrixXd hyperparameters, Eigen::MatrixXd beta, double epsilon);

  Eigen::MatrixXd kernel_function(Eigen::MatrixXd A, Eigen::MatrixXd B, int obs_index);
  Eigen::MatrixXd regression_linear_function(Eigen::MatrixXd X_mat, int obs_index);
  Eigen::MatrixXd emulate(arma::mat X_s, Eigen::MatrixXd Y_mat);
  Eigen::MatrixXd emulate_nr(Eigen::MatrixXd X_s, Eigen::VectorXd y);
};

#endif
