#include "mcmc.h"

/*SETUP--
***************/
MCMC::MCMC(){
  this->maxLogLikelihood = -std::numeric_limits<double>::infinity();
  this->Likelihood = maxLogLikelihood;
  this->normal_dist = std::normal_distribution<double>(0.0,1.0);
  this->uniform_dist = std::uniform_real_distribution<double>(0.0,1.0);

  this->setSeed(1);
}
MCMC::MCMC(Eigen::MatrixXd newTargetValue, Eigen::MatrixXd newRange, Eigen::VectorXd newWidths, bool nR){
  this->NR = nR;
  this->maxLogLikelihood = -std::numeric_limits<double>::infinity();
  this->Likelihood = maxLogLikelihood;
  this->setSeed(1);
  this->targetValue = newTargetValue;
  this->obsCount = this->targetValue.cols();
  this->Range = newRange;
  this->paramCount = this->Range.rows();
  this->Position = Eigen::MatrixXd::Zero(1,paramCount);
  this->testPosition = Eigen::MatrixXd::Zero(1,paramCount);
  this->Widths = newWidths;
}
void MCMC::setTargetValue(Eigen::MatrixXd newTargetValue){
  this->targetValue = newTargetValue;
  this->obsCount = this->targetValue.cols();
}
void MCMC::setRange(Eigen::MatrixXd newRange){
  this->Range = newRange;
  this->paramCount = this->Range.rows();
  this->Position = Eigen::MatrixXd::Zero(1,paramCount);
  this->testPosition = Eigen::MatrixXd::Zero(1,paramCount);
}
void MCMC::setPosition(Eigen::MatrixXd newPosition){
  this->Position = newPosition;
}
void MCMC::setPosition(){
  for(int i=0;i<this->paramCount;i++){
    this->Position(0,i) = (this->Range(i,0) + this->Range(i,1))/2.0;
  }
}
void MCMC::setWidths(Eigen::VectorXd newWidths){
  this->Widths = newWidths;
}

/*RANDOM--
***************/
void MCMC::setSeed(unsigned seed){
  this->generator.seed(seed);
}
void MCMC::setSeedClock(){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator.seed(seed);
}
double MCMC::normal(){
  return this->normal_dist(generator);
}
double MCMC::uniform(){
  return this->uniform_dist(generator);
}

/*STEP--
**************/
double MCMC::getLogLikelihood(Eigen::MatrixXd Z){
  double z=0.0, diff=0.0;
  for(int i=0;i<this->obsCount;i++){
    diff = (Z(0,i) - this->targetValue(0,i));
    z += diff*diff;
  }
  z *= 0.5;
  return -z;
}
void MCMC::step(){
  for(int i=0;i<this->paramCount;i++){
    testPosition(0,i) = Position(0,i) + normal()*this->Widths(i);
  }
}
bool MCMC::decide(Eigen::MatrixXd Z){
  double random, LH, ratio;
  for(int i=0;i<paramCount;i++){
    if(this->testPosition(0,i) < this->Range(i,0)){
      return false;
    } else if(this->testPosition(0,i) > this->Range(i,1)){
      return false;
    }
  }
  LH = this->getLogLikelihood(Z);

  ratio = LH - this->Likelihood;
  if(ratio>0){
    this->Position = this->testPosition;
    this->Likelihood = LH;
    if(this->Likelihood > this->maxLogLikelihood){
      this->maxLogLikelihood = this->Likelihood;
      printf("MaxLogLH: %f\n",this->maxLogLikelihood);
    }
    return true;
  }else{
    random = uniform();
    if(log(random)<ratio){
      this->Position = this->testPosition;
      this->Likelihood = LH;
      if(this->Likelihood > this->maxLogLikelihood){
	this->maxLogLikelihood = this->Likelihood;
	printf("MaxLogLH: %f\n",this->maxLogLikelihood);
      }
      return true;
    }else{
      return false;
    }
  }
}
void MCMC::Run(int Samples, Eigen::MatrixXd &History, emulator obsEmulator, Eigen::MatrixXd Y){
  Eigen::MatrixXd testValue;
  bool takeStep=false;
  int print=Samples/10,
    percent=0;
  double acceptedSteps = 0.0;

  History = Eigen::MatrixXd::Zero(Samples,this->paramCount+this->obsCount);

  for(int step=0;step<Samples;step++){
    this->step();
    if(this->NR==true){
      obsEmulator.Emulate_NR(this->testPosition,Y,testValue);
    }else{
      obsEmulator.Emulate(this->testPosition,Y,testValue);
    }      
    takeStep = this->decide(testValue);
    History.block(step,0,1,this->paramCount) = this->Position;
    History.block(step,this->paramCount,1,this->obsCount) = testValue;

    if(takeStep==true){
      acceptedSteps++;
    }
    if((step+1)%print==0){
      percent+=10;
      printf("%d%% completed...%f%% steps accepted\n",percent,acceptedSteps*100/((double)Samples*percent/100));
    }
  }
  acceptedSteps *= 100.0/(double) Samples;
  printf("Percent Acceptance: %f%%...\n",acceptedSteps);
  printf("MaxLLH: %f...\n",this->maxLogLikelihood);
}
    

    
  
  
