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

double func(double x)
{
  //return sin(6.28*x);
  return sin(3.14*x);
  //return sin(x);
  //return (6.0*x - 2.0)*(6.0*x - 2.0)*sin(12.0*x - 4.0);
}

int main(int argc, char* argv[]){


  /***********testing*********/
  bool testing=true;
  if(testing){
  int points = 23,
    test = 50;
  double start=-7.5,
    end=7.5,
    dx=(end-start)/(double) (points),
    x=0;

  Eigen::MatrixXd function = Eigen::MatrixXd::Zero(points,2);
  Eigen::MatrixXd hyperparam = Eigen::MatrixXd::Zero(1,3);
  Eigen::MatrixXd testX = Eigen::MatrixXd::Zero(test,1);
  Eigen::MatrixXd testY;
  Eigen::MatrixXd testOut = Eigen::MatrixXd::Zero(test,4);
  if( argc != 4 )
    {
      printf("Usage: r l points\n");
      return 1;
    }

  hyperparam(0,0) = atof(argv[1]);
  hyperparam(0,1) = atof(argv[2]);
  hyperparam(0,2) = 0.1;
  points = atoi(argv[3]);
  dx=(end-start)/(double) (points);

 for(int i = 0; i < points; i++){
    function(i,0) = x = start + dx*(double) i;
    //function(i,1) = func(x);

  }
  dx = (end - start)/(double) (test);
 for(int i = 0; i < test; i++){
    testX(i,0) = start + dx*(double) i;
  }

  emulator emul(function.col(0),hyperparam);
  emul.Emulate_NoPrior(function.col(0),testY);
  function.col(1) = testY.col(0);
  emul.setPrint(true);
  emul.Emulate_NR(testX,function.col(1),testY);

  testOut.block(0,0,test,1) = testX;
  testOut.block(0,1,test,1) = testY.col(0);

  for(int i=0;i<test;i++){
    testOut(i,2) = testY(i,0) - 2.0*sqrt(testY(i,1));
    testOut(i,3) = testY(i,0) + 2.0*sqrt(testY(i,1));
  }
  WriteFile("test.dat",testOut," ");
  WriteFile("train.dat",function," ");

  emul.setPrint(true);
  int Row=10.0/0.2;
  int Col=5.0/0.1;
  Eigen::MatrixXd logprob = Eigen::MatrixXd::Zero(Row,Col);
  std::fstream ofile("log.dat",std::ios::out);
  for( double r = 0.1; r < 10.0; r += 0.2 )
    {
      int R = r/0.2;
      hyperparam(0,0) = r;
      for( double l = 0.05; l < 5.0; l += 0.1 )
	{
	  int L = l/0.1;
	  hyperparam(0,1) = l;
	  for( double n = 0.1; n < 0.2; n += 0.5)
	    {
	      //int N = n/0.05;
	      hyperparam(0,2) = n;
	      emul.Construct(function.col(0),hyperparam);
	      emul.setPrint(true);
	      logprob(R,L) = emul.Emulate_NR(testX,function.col(1),testY);
	      ofile << r << " " << l << " " << n << " " << logprob(R,L) << '\n';

	    }
	}
    }
  ofile.close();

  Eigen::MatrixXf::Index maxRow, maxCol;
  float max = logprob.maxCoeff(&maxRow, &maxCol);
  printf("%f %f %f\n",max,0.1+0.2*(double)maxRow,0.05 + 0.1*(double)maxCol);
}



  /********COMPUTATION*********/
  if(!testing){
  std::vector<std::string> modelfilenames={"I211_J211.dat","I2212_J2212.dat","I321_J2212.dat","I321_J321.dat"},
    expfilenames={"star_pipi.dat","star_ppbar.dat","star_pK.dat","star_KK.dat"},
      paramNames;

  std::string foldername="model_output",
    writefilename="parameters.dat",
    rangename=foldername+"/parameter_priors.dat",
    delimiter=" ",
    infilename=foldername+"/moments_parameters.dat",
    outfilename=foldername+"/wave.dat";

  if(argc != 3){
    printf("Usage: ./main start finish\n");
    return 1;
  }

  int start=atoi(argv[1]),
    finish=atoi(argv[2]),
    column=1,
    parameters,
    observables,
    ab=4,
    n=500,
    nmax=0,
    modelfiles = modelfilenames.size(),
    expfiles = expfilenames.size();
  int samples=finish-start;
  double x=8;

  std::vector<Eigen::MatrixXd> ModelMatrix(modelfiles),
    ExpMatrix(expfiles);
  Eigen::MatrixXd Dy,ModelObs,ExpObs,Parameters,GAB,range;
  Eigen::VectorXd modeldy, expdy;

  //Write Parameters files
  WriteParameterFiles(rangename, foldername, writefilename, delimiter, start, finish, ab, Parameters);
  LoadParamFile(rangename,paramNames,range,delimiter);
  WriteFile(infilename,Parameters,delimiter);


  //Construct gab functions
  /*
  CDistCosh
    dist;
  LoadFile(infilename, Parameters, delimiter);
  dist.FunctionSet(n,x,samples,ab,nmax,Parameters,GAB);
  WriteFile(outfilename,GAB,delimiter);
  */
  //Run Scott's program
  //system("cd ..;pwd; cd src/");
  

  //Load Data
  LoadDataFile(foldername, modelfilenames[0], delimiter, 0, 1, 0, Dy);
  modeldy = Dy.col(0);
  LoadDataFile(foldername, expfilenames[0], delimiter, 0, 1, 0, Dy);
  expdy = Dy.col(0);
  
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

  printf("Removing first row of model data and first row of exp data...\n");
  for(int i=0;i<ModelMatrix.size();i++){
    RemoveRow(ModelMatrix[i],0);
  }
  for(int i=0;i<ExpMatrix.size();i++){
    RemoveRow(ExpMatrix[i],0);
  }


  MatrixMoments(ModelMatrix,modeldy,ModelObs);
  MatrixMoments(ExpMatrix,expdy,ExpObs);
  Error = er*ExpObs.col(0);

  AverageRows(Mean,ModelObs);
  TildeFunction(ModelTilde,Mean,Error,ModelObs);
  TildeFunction(ExpTilde,Mean,Error,ExpObs);

  CovarianceFunction(Covariance,ModelTilde);
  EigenSolve(EigenValues,EigenVectors,Covariance);
  EigenSort(EigenValues,EigenVectors);
  ModelZ = ModelTilde.transpose()*EigenVectors;
  ExpZ = ExpTilde.transpose()*EigenVectors;

  WriteFile(foldername+"/expz.dat",ExpZ," ");

  /*
  std::cout << "Parameters:\n" << Parameters << std::endl;  
  std::cout << "ModelObs:\n" << ModelObs << std::endl;
  std::cout << "ExpObs:\n" << ExpObs << std::endl;
  std::cout << "Error:\n" << Error << std::endl;
  std::cout << "Mean:\n" << Mean << std::endl;
  std::cout << "ModelTilde:\n" << ModelTilde << std::endl;
  std::cout << "ExpTilde:\n" << ExpTilde << std::endl;
  std::cout << "Covariance:\n" << Covariance << std::endl;
  std::cout << "EV:\n" << EigenValues.transpose() << std::endl;
  std::cout << "EV:\n" << EigenVectors << std::endl;
  std::cout << "----------------\n" << std::endl;
  std::cout << "ModelZ:\n" << ModelZ << std::endl;
  AverageColumns(Mean,ModelZ);
  std::cout << "MMean:\n" << Mean.transpose() << std::endl;
  std::cout << "ExpZ:\n" << ExpZ << std::endl;
  */
  observables = ModelZ.cols();
  parameters = Parameters.cols();
  int test=finish-start;
  Eigen::MatrixXd Hyperparameters,
    outMatrix(test,parameters+observables),
    plot(finish-start,parameters+observables),
    testX = Parameters,
    testY;

  plot.block(0,0,finish-start,parameters) = Parameters;
  plot.block(0,parameters,finish-start,observables) = ModelZ;
  WriteFile("trainplot.dat",plot," ");

  Eigen::VectorXd Width(parameters);
  for(int i=0;i<parameters;i++){
    Width(i) = (range(i,1) - range(i,0))/50.0;
  }
  LoadFile("hyperparameters.dat",Hyperparameters," ");


  Eigen::MatrixXd Beta;
  linearRegressionLeastSquares(ModelZ,Parameters,Beta);
  WriteFile("beta.dat",Beta," ");


  /*  
  int row = atoi(argv[3]);
  testX = Parameters.row(row);
  Eigen::MatrixXd testY_ = ModelZ.row(row);
  RemoveRow(Parameters,row);
  RemoveRow(ModelZ,row);
  std::cout << "testX\n" << testX << std::endl;
  std::cout << "ModelY\n" << testY_ << std::endl;

  Eigen::MatrixXd Fit;
  emulator emulation(Parameters, Hyperparameters, Beta);
  outMatrix = Eigen::MatrixXd::Zero(1,parameters + observables);
  emulation.Emulate(testX, ModelZ, testY);
  test = 1;
  outMatrix.block(0,0,test,parameters) = testX;
  outMatrix.block(0,parameters,test,observables) = testY;
  std::cout << "WR_testY\n" << testY << std::endl;
  ComputeFit(testY_,testY,Fit);
  std::cout << Fit << std::endl;
  */




  test = 27;
  LHCSampling(testX, test, 4, range);
  emulator emulation(Parameters, Hyperparameters.row(0), Beta.row(0));

  emulation.setPrint(true);

  int rInt = 80.0/5.0-2;
  int lInt = 0.5/0.05;
  Eigen::MatrixXd logprob = Eigen::MatrixXd::Zero(rInt,lInt);
  Eigen::MatrixXd hyperparam = Eigen::MatrixXd::Zero(1,3);

  std::fstream ofile("log.dat",std::ios::out);
  for( double r = 10; r < 80.0; r += 5 )
    {
      int R = r/5.0 - 2;
      hyperparam(0,0) = r;
      for( double l = 0.01; l < 0.5; l += 0.05 )
	{
	  int L = l/0.05;
	  hyperparam(0,1) = l;
	  for( double n = 0.1; n < 0.5; n += 0.5)
	    {
	      hyperparam(0,2) = n;
	      emulation.Construct(Parameters,hyperparam.row(0));
	      emulation.setPrint(true);
	      logprob(R,L) = emulation.Emulate_NR(testX,ModelZ.col(0),testY);
	      ofile << r << " " << l << " " << n << " " << logprob(R,L) << '\n';
	      std::cout << r << " " << l << " " << n << " " << logprob(R,L) << '\n';
	    }
	}
    }
  ofile.close();

  Eigen::MatrixXf::Index maxRow, maxCol;
  float max = logprob.maxCoeff(&maxRow, &maxCol);
  printf("%f %f %f\n",max,10+5*(double)maxRow,0.01+ 0.05*(double)maxCol);


  /*
  std::cout << "NR_testY\n" << testY << std::endl;
  ComputeFit(testY_,testY,Fit);
  std::cout << Fit << std::endl;
  */

  /*  


  ComputeFit(ModelZ,testY,Fit);
  WriteFile("fit.dat",Fit," ");
  */
  /*
  MCMC mcmc(ExpZ,range,Width,false);
  mcmc.setPosition();
  int Samples=10000;
  
  Eigen::MatrixXd History;
  mcmc.Run(Samples,History,emulation,ModelZ);
  Eigen::MatrixXd trace = History.block(0,0,Samples,parameters);
  WriteCSVFile("mcmctrace_50.csv",paramNames,trace,",");
  */
  }
  return 0;
}
