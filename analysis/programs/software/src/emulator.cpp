#include "emulator.h"

emulator::emulator(Eigen::MatrixXd newX, Eigen::MatrixXd newHyperparam, Eigen::MatrixXd newBeta){
  this->X = newX;
  this->Hyperparam = newHyperparam;
  this->Beta = newBeta;

  this->trainPoints = this->X.rows();
  this->paramCount = this->X.cols();
  this->hyperparamCount = this->Hyperparam.cols();
  this->obsCount = this->Beta.rows();

  this->Noise = this->Hyperparam.col(this->hyperparamCount - 1);

  Eigen::MatrixXd trainI = Eigen::MatrixXd::Identity(trainPoints, trainPoints);

  this->Kernel = std::vector<Eigen::MatrixXd>(obsCount);
  this->KernelInv = std::vector<Eigen::MatrixXd>(obsCount);
  this->HMatrix = std::vector<Eigen::MatrixXd>(obsCount);
  for(int ob=0;ob<obsCount;ob++){
    this->kernelFunction(this->X,this->X,ob,this->Kernel[ob]);
    this->Kernel[ob] += trainI*(this->Noise(ob)*Noise(ob) + this->epsilon);
    this->KernelInv[ob] = Kernel[ob].inverse();
    regressionLinearFunction(this->X, ob,  this->HMatrix[ob]);
  }
}
void emulator::kernelFunction(Eigen::MatrixXd A, Eigen::MatrixXd B, int obsIndex, Eigen::MatrixXd &kernelMatrix){
  double xi, xj, sum, diff,
    theta1 = this->Hyperparam(obsIndex,0)*this->Hyperparam(obsIndex,0),
    theta2 = 0.5/(this->Hyperparam(obsIndex,1)*this->Hyperparam(obsIndex,1));
  int ai = A.rows(),
    bi = B.rows();
  kernelMatrix = Eigen::MatrixXd::Zero(ai,bi);

  for(int a=0;a<ai;a++){
    for(int b=0;b<bi;b++){
      sum = 0.0;
      for(int pc=0;pc<paramCount;pc++){
	xi = A(a,pc);
	xj = B(b,pc);
	diff = xi - xj;
	sum += diff*diff;
      }
      kernelMatrix(a,b) = theta1*exp(-sum*theta2);
    }
  }
}
void emulator::regressionLinearFunction(Eigen::MatrixXd testX, int obsIndex, Eigen::MatrixXd &HMatrix){
  int points = testX.rows();
  HMatrix = Eigen::MatrixXd::Zero(1,points);
  for(int i=0;i<points;i++){
    HMatrix(0,i) = this->Beta(obsIndex,0);
    for(int j=0;j<this->paramCount;j++){
      HMatrix(0,i) += Beta(obsIndex,j+1)*testX(i,j);
    }
  }
}
void emulator::Emulate(Eigen::MatrixXd testX, Eigen::MatrixXd Y, Eigen::MatrixXd &outMatrix){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> dist(0.0,1.0);

  int testPoints = testX.rows();
  
  //matrices and vectors for computational use
  Eigen::MatrixXd KernelS = Eigen::MatrixXd::Zero(trainPoints,testPoints),
    KernelSS = Eigen::MatrixXd::Zero(testPoints,testPoints),
    
    meanMatrix = Eigen::MatrixXd::Zero(testPoints,1),
    varMatrix = Eigen::MatrixXd::Zero(testPoints, testPoints),
    L = Eigen::MatrixXd::Zero(testPoints, testPoints),

    HS = Eigen::MatrixXd::Zero(1, testPoints),
    R = Eigen::MatrixXd::Zero(trainPoints, testPoints),
    betaMatrix = Eigen::MatrixXd::Zero(trainPoints, 1),

    tempPP = Eigen::MatrixXd::Zero(1,1),
    tempPPInv = Eigen::MatrixXd::Zero(1,1),

    postFunc = Eigen::MatrixXd::Zero(testPoints, 1),
    randomSample = Eigen::MatrixXd::Zero(testPoints, 1),

    testI = Eigen::MatrixXd::Identity(testPoints,testPoints);

  outMatrix = Eigen::MatrixXd::Zero(testPoints,this->paramCount+this->obsCount*5);

  for(int ob=0;ob<obsCount;ob++){
    kernelFunction(this->X, testX, ob, KernelS);
    kernelFunction(testX, testX, ob, KernelSS);
    KernelSS += testI*this->epsilon;

    regressionLinearFunction(testX, ob, HS);
    R = HS - this->HMatrix[ob]*this->KernelInv[ob]*KernelS;
    tempPP = this->HMatrix[ob]*KernelInv[ob]*this->HMatrix[ob].transpose();
    tempPPInv = tempPP.inverse();
    betaMatrix = tempPPInv*this->HMatrix[ob]*this->KernelInv[ob]*Y.col(ob);

    meanMatrix = KernelS.transpose()*this->KernelInv[ob]*Y.col(ob) + R.transpose()*betaMatrix;

    varMatrix = KernelSS - KernelS.transpose()*this->KernelInv[ob]*KernelS + R.transpose()*tempPPInv*R;
    L = varMatrix.llt().matrixL();

    for(int tp=0;tp<testPoints;tp++){
      randomSample(tp,0) = dist(generator);
    }
    postFunc = meanMatrix + L*randomSample;

    outMatrix.col(paramCount + ob) = meanMatrix.col(0);
    outMatrix.col(paramCount + this->obsCount*3 + ob) = postFunc.col(0);
    for(int tp=0;tp<testPoints;tp++){
      outMatrix(tp,paramCount + this->obsCount + ob) = meanMatrix(tp,0) - 2.0*sqrt(varMatrix(tp,tp));
      outMatrix(tp,paramCount + this->obsCount*2 + ob) = meanMatrix(tp,0) + 2.0*sqrt(varMatrix(tp,tp));
      outMatrix(tp,paramCount + this->obsCount*4 + ob) = sqrt(varMatrix(tp,tp));
    }
  }

  for(int param=0;param<paramCount;param++){
    outMatrix.col(param) = testX.col(param);
  }
}
