#include "dft.hpp"

DFT::DFT(){
}
DFT::DFT(int N){
  this->nx = N;
  this->no = N;

  D = Eigen::MatrixXcd::Zero(this->no+1,this->nx+1);
  Dstar = Eigen::MatrixXcd::Zero(this->nx+1,this->no+1);
}
DFT::DFT(int NX, int NO){
  this->nx = NX;
  this->no = NO;

  D = Eigen::MatrixXcd::Zero(this->no+1,this->nx+1);
  Dstar = Eigen::MatrixXcd::Zero(this->nx+1,this->no+1);
}
void DFT::set_limit(double start, double end){
  this->L.resize(2);
  this->L[0] = start;
  this->L[1] = end;
  
  this->D.resize(this->no+1,this->nx+1);
  this->Dstar.resize(this->nx+1,this->no+1);

  double
    twopi_nx = 8.0*atan(1)/(double) this->nx,
    omega;
  std::complex<double> 
    z;

  for(int row=0; row<this->no+1; row++){
    omega = -0.5*(this->no) + row;
    for(int col=0; col<this->nx+1; col++){
      z = std::complex<double>(0.0,twopi_nx*omega*(double) col);
      D(row,col) = exp(z);
      Dstar(col,row) = exp(-z);
    }
  }
}
std::complex<double> f_function(double x){
  double sigma = (x-700.0)/10.0;
  sigma *= sigma;
  return (1 - sigma)*exp(-0.5*sigma);
}
void project_function(Eigen::VectorXcd &f, int &nx, double &L){
  double 
    sigma=0, 
    x=0, 
    dx=L/(double) (nx);
  f.resize(nx+1);

  for(int row=0;row<nx+1;row++){
    x = dx*(double) row;
    if(x<600){
    } else if(x>800){
    } else{
      sigma = (x-700.0)/10.0;
      sigma *= sigma;
      f(row) = (1 - sigma)*exp(-0.5*sigma);
    }
  }
}
void write(Eigen::VectorXcd &f, Eigen::VectorXcd &ftest, Eigen::VectorXcd &F, double &L, std::string func_name, std::string freq_name){
  int nx = f.size();
  int no = F.size();
  double dx = L/(double) nx;

  std::ofstream ofile;
  ofile.open(func_name);
  for(int row=0;row<nx;row++){
    ofile << dx*(double)row << " " <<  f(row).real() << " " << ftest(row).real() << " " << ftest(row).imag() << std::endl;
  }
  ofile.close();

  ofile.open(freq_name);
  for(int row=0;row<no;row++){
    ofile << (-no/2 + row)/L << " " << std::abs(F(row)) << std::endl;
  }
  ofile.close();
}
void load(Eigen::MatrixXd &M, int Nx, int Nt, std::string name){
  M.resize(Nx+1,Nt+2);
  std::ifstream ifile;
  ifile.open(name);
  for(int row=0;row<Nx+1;row++){
    for(int col=0;col<Nt+2;col++){
      ifile >> M(row,col);
    }
  }
  ifile.close();
}
