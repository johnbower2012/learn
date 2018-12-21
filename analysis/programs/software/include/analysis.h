#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include<Eigen/Dense>
#include<vector>
#include<iostream>

void RemoveRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void RemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
void AddOnesRow(Eigen::MatrixXd matrix, Eigen::MatrixXd &outMatrix);
void AddOnesColumn(Eigen::MatrixXd matrix, Eigen::MatrixXd &outMatrix);

void ScaleMatrixColumns(Eigen::MatrixXd Matrix,Eigen::VectorXd &Mean, Eigen::VectorXd &Std, Eigen::MatrixXd &Scaled);
void ScaleMatrixRows(Eigen::MatrixXd Matrix,Eigen::VectorXd &Mean, Eigen::VectorXd &Stdd, Eigen::MatrixXd &Scaled);

void AverageColumns(Eigen::VectorXd &average, Eigen::MatrixXd matrix);
void AverageRows(Eigen::VectorXd &average, Eigen::MatrixXd matrix);
void TildeFunction(Eigen::MatrixXd &tilde, Eigen::VectorXd mean, Eigen::VectorXd error, Eigen::MatrixXd matrix);
void CovarianceFunction(Eigen::MatrixXd &cov, Eigen::MatrixXd matrix);
void EigenSolve(Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors, Eigen::MatrixXd matrix);
void FindMax(Eigen::VectorXd vector, int &index, double &max);
void FindMax(Eigen::MatrixXd vector, int &index, double &max);
void EigenSort(Eigen::VectorXd &eigval, Eigen::MatrixXd &eigvec);

void ZerothMoment(Eigen::MatrixXd Function, double &Moment);
void FirstMoment(Eigen::MatrixXd Function, double &Moment);
void SecondMoment(Eigen::MatrixXd Function, double &Moment);
void MatrixMoments(std::vector<Eigen::MatrixXd> Matrix, Eigen::VectorXd DelX, Eigen::MatrixXd &Obs);
void linearRegressionLeastSquares(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd &Beta);
void ComputeFit(Eigen::MatrixXd ModelZ, Eigen::MatrixXd EmulatorZ, Eigen::MatrixXd &outMatrix);

#endif
