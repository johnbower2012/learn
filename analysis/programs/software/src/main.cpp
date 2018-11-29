#include<iostream>
#include<Eigen/Dense>
#include<string>
#include<vector>
#include<stdlib.h>
#include "system.h"
#include "analysis.h"
#include "coshfunc.h"
#include "emulator.h"
#include "mcmc.h"

int main(int argc, char* argv[]){
  std::vector<std::string> modelfilenames={"I211_J211.dat","I2212_J2212.dat","I321_J2212.dat","I321_J321.dat"},
    expfilenames={"star_pipi.dat","star_ppbar.dat","star_pK.dat","star_KK.dat"},
      paramNames;

  std::string foldername="model_output",
    writefilename="parameters.dat",
    rangename=foldername+"/parameter_priors.dat",
    delimiter=" ",
    infilename=foldername+"/moments_parameters.dat",
    outfilename=foldername+"/wave.dat";

  int start=atoi(argv[1]),
    finish=atoi(argv[2]),
    column=1,
    parameters,
    observables,
    ab=4,
    n=500,
    nmax=4,
    modelfiles = modelfilenames.size(),
    expfiles = expfilenames.size();
  int samples=finish-start;
  double x=8;

  std::vector<Eigen::MatrixXd> ModelMatrix(modelfiles),
    ExpMatrix(expfiles);
  Eigen::MatrixXd Dy,ModelObs,ExpObs,Parameters,GAB,range;
  Eigen::VectorXd modeldy, expdy;

  //Load Data
  LoadDataFile(foldername, modelfilenames[0], delimiter, 0, 1, 0, Dy);
  modeldy = Dy.col(0);
  LoadDataFile(foldername, expfilenames[0], delimiter, 0, 1, 0, Dy);
  expdy = Dy.col(0);

  //Write Parameters files
  WriteParameterFiles(rangename, foldername, writefilename, delimiter, start, finish, ab,Parameters);
  LoadParamFile(rangename,paramNames,range,delimiter);
  WriteFile(infilename,Parameters,delimiter);

  //Construct gab functions
  CDistCosh
    dist;
  LoadFile(infilename, Parameters, delimiter);
  dist.FunctionSet(n,x,samples,ab,nmax,Parameters,GAB);
  WriteFile(outfilename,GAB,delimiter);

  //Run Scott's program
  //system("cd ..;pwd; cd src/");

  //Conduct Observables Analysis
  double er=0.1;
  Eigen::VectorXd Error,
    Mean,
    EigenValues;
  Eigen::MatrixXd ModelTilde,
    ExpTilde,
    Covariance,
    EigenVectors,
    ModelZ,
    ExpZ;

  LoadDataFiles(foldername, modelfilenames, delimiter, start, finish, column, ModelMatrix);
  LoadDataFiles(foldername, expfilenames, delimiter, 0, 1, column, ExpMatrix);

  MatrixMoments(ModelMatrix,modeldy,ModelObs);
  MatrixMoments(ExpMatrix,expdy,ExpObs);
  Error = er*ExpObs.col(0);

  AverageRows(Mean,ModelObs);
  TildeFunction(ModelTilde,Mean,Error,ModelObs);
  TildeFunction(ExpTilde,Mean,Error,ExpObs);

  CovarianceFunction(Covariance,ModelObs);
  EigenSolve(EigenValues,EigenVectors,Covariance);
  EigenSort(EigenValues,EigenVectors);
  ModelZ = ModelObs.transpose()*EigenVectors;
  ExpZ = ExpObs.transpose()*EigenVectors;
  //std::cout << "Parameters:\n" << Parameters << std::endl;  
  //std::cout << "ModelObs:\n" << ModelObs << std::endl;
  //std::cout << "ExpObs:\n" << ExpObs << std::endl;
  //std::cout << "Error:\n" << Error << std::endl;
  //std::cout << "Mean:\n" << Mean << std::endl;
  //std::cout << "ModelTilde:\n" << ModelTilde << std::endl;
  //std::cout << "ExpTilde:\n" << ExpTilde << std::endl;
  //std::cout << "Covariance:\n" << Covariance << std::endl;
  //std::cout << "EV:\n" << EigenValues << std::endl;
  //std::cout << "EV:\n" << EigenVectors << std::endl;
  //std::cout << "----------------\n" << std::endl;
  //std::cout << "ModelZ:\n" << ModelZ << std::endl;
  //AverageColumns(Mean,ModelZ);
  //std::cout << "MMean:\n" << Mean.transpose() << std::endl;
  //std::cout << "ExpZ:\n" << ExpZ << std::endl;

  observables = ModelZ.cols();
  parameters = Parameters.cols();

  Eigen::MatrixXd Hyperparameters,
    outMatrix,
    plot(finish-start,parameters+observables);
  Eigen::VectorXd Width(parameters);

  for(int i=0;i<parameters;i++){
    Width(i) = (range(i,1) - range(i,0))/10.0;
  }

  for(int i=0;i<finish-start;i++){
    for(int j=0;j<parameters;j++){
      plot(i,j) = Parameters(i,j);
    }
    for(int j=0;j<observables;j++){
      plot(i,j+parameters) = ModelZ(i,j);
    }
  }
  WriteFile("plot.dat",plot," ");
  LoadFile("hyperparameters.dat",Hyperparameters," ");

  emulator emulation(Parameters, Hyperparameters);
  MCMC mcmc(ExpZ,range,Width);
  mcmc.setPosition();
  int Samples=10000;
  
  Eigen::MatrixXd History;
  mcmc.Run(Samples,History,emulation,ModelZ);
  Eigen::MatrixXd trace = History.block(0,0,Samples,parameters);
  WriteCSVFile("mcmctrace.csv",paramNames,trace,",");

  return 0;
}