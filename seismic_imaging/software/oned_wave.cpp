#include<iostream>
#include<Eigen/Dense>
#include<cmath>
#include<fstream>
#include<string>
#include<ctime>
#include "../lib/dft.cpp"

void project_function(Eigen::VectorXd &F, int &Nx, double &L);
void project_velocity(Eigen::VectorXd &F, int &Nx, double &L);
void project_B(Eigen::MatrixXd &B, Eigen::VectorXd &C, double Dx, double Dt);
void write_matrix(Eigen::MatrixXd &A,std::string name);

int main(int argc, char* argv[]){
  int
    Nx=2000,
    Nt;

  double
    L=1000.0,
    dx=L/(double) Nx,
    //dt=dx/3000.0, #currently evaulates to 0.0001666667
    dt=0.0001,
    time;

  std::clock_t start, finish;

  Nt = 5.0/dt;

  Eigen::VectorXd 
    u0,
    dudt0,
    c,
    b;

  Eigen::MatrixXd B,U;

  U = Eigen::MatrixXd::Zero(Nx+1,Nt+2);

  project_function(u0,Nx,L);
  dudt0=u0;
  project_velocity(c,Nx,L);
  project_B(B,c,dx,dt);

  for(int xi=0;xi<Nx+1;xi++){
    U(xi,0) = dx*(double) xi;
  }
  U.col(1) = u0;
  U.col(2) = 0.5*B*u0 + dt*dudt0;

  start = clock();
  for(int ti=2;ti<Nt+1;ti++){
    //Matrix method, slow:
    //    U.col(ti+1) = -U.col(ti-1) + B*U.col(ti);

    U(0,ti+1) = 0.0; //D boundary x=0
    //    U(0,ti+1) = -U(0,ti-1) + 2.0*U(0,ti) + 2.0*(c(0)*c(0)*dt*dt/dx/dx)*(U(1,ti) - U(0,ti)); //VN boundary x=0
    for(int xi=1;xi<Nx;xi++){
      U(xi,ti+1) = -U(xi,ti-1) + B(xi,xi-1)*U(xi-1,ti) + B(xi,xi)*U(xi,ti) + B(xi,xi+1)*U(xi+1,ti);
    }
    //    U(Nx,ti+1)=0.0; //D boundary x=L
    U(Nx,ti+1) = -U(Nx,ti-1) + 2.0*U(Nx,ti) + 2.0*(c(Nx)*c(Nx)*dt*dt/dx/dx)*(U(Nx-1,ti) - U(Nx,ti));//VN boundary x=L
    if(ti%1000==0){
      printf("%d / %d\n",ti,Nt);
    }
  }
  finish = clock();
  time = (finish - start)/(double) CLOCKS_PER_SEC;
  printf("time: %f\n",time);

  std::string name = "../data/wave.dat";
  write_matrix(U,name);
  name = "../data/u(5).dat";
  Eigen::MatrixXd output(Nx+1,2);
  output.col(0) = U.col(0);
  output.col(1) = U.col(Nt+1);
  write_matrix(output,name);

  DFT dft(Nx);
  std::string t0_freq = "../data/freq_t0.dat";
  std::string tNt_freq = "../data/freq_tNt.dat";
  std::string t0_func = "../data/func_t0.dat";
  std::string tNt_func = "../data/func_tNt.dat";
  Eigen::VectorXcd f,F,ftest;
  Eigen::MatrixXcd D,Dstar;
  project_function(f,Nx,L);
  dft.set_limit(0,L);
  F = dft.D*f;
  ftest = dft.Dstar*F/(double) (Nx+1);

  //print to file
  start = clock();
  write(f,ftest,F,L,t0_func,t0_freq);
  finish = clock();
  time = (finish - start)/(double) CLOCKS_PER_SEC;
  printf("time: %f\n",time);  

  f = U.col(Nt+1);
  F = dft.D*f;
  ftest = dft.Dstar*F/(double) (Nx+1);

  //print to file
  start = clock();
  write(f,ftest,F,L,tNt_func,tNt_freq);
  finish = clock();
  time = (finish - start)/(double) CLOCKS_PER_SEC;
  printf("time: %f\n",time);  


  return 0;
}

void project_function(Eigen::VectorXd &F, int &Nx, double &L){
  double 
    sigma=0, 
    x=0, 
    dx=L/(double) Nx;
  F.resize(Nx+1);

  for(int row=0;row<Nx+1;row++){
    x = dx*(double) row;
    if(x<600){
    } else if(x>800){
    } else{
      sigma = (x-700.0)/10.0;
      sigma *= sigma;
      F(row) = (1 - sigma)*exp(-0.5*sigma);
    }
  }
}
void project_velocity(Eigen::VectorXd &F, int &Nx, double &L){
  double 
    x=0, 
    dx=L/(double) Nx;
  F.resize(Nx+1);

  for(int row=0;row<Nx+1;row++){
    x = dx*(double) row;
    if(x<500){
      F(row) = 3000;
    } else{
      F(row) = 1000;
    }
  }
}
void project_B(Eigen::MatrixXd &B, Eigen::VectorXd &C, double Dx, double Dt){
  int Nx=C.size()-1;
  double sigma=0;
  B=Eigen::MatrixXd::Zero(Nx+1,Nx+1);
  //Will leave B(0,0) = B(Nx,Nx) = 0 for D boundary conditions
  for(int ix=1;ix<Nx;ix++){
    sigma=C(ix-1)*Dt/Dx;
    sigma*=sigma;
    B(ix,ix-1) = sigma;

    sigma=C(ix)*Dt/Dx;
    sigma*=sigma;
    B(ix,ix) = 2.0*(1.0-sigma);

    sigma=C(ix+1)*Dt/Dx;
    sigma*=sigma;
    B(ix,ix+1) = sigma;
  }
  //VN condition x=0
  /*
  sigma=C(0)*Dt/Dx;
  sigma*=sigma;
  B(0,0) = 2.0*(1.0-sigma);
  sigma=C(1)*Dt/Dx;
  sigma*=sigma;
  B(0,1) = 2.0*sigma;
  */

  //VN condition x=L
  sigma=C(Nx-1)*Dt/Dx;
  sigma*=sigma;
  B(Nx,Nx-1) = 2.0*sigma;
  sigma=C(Nx)*Dt/Dx;
  sigma*=sigma;
  B(Nx,Nx) = 2.0*(1.0-sigma);
}
void write_matrix(Eigen::MatrixXd &A,std::string name){
  int ni=A.rows();
  int nj=A.cols();
  std::cout << nj << std::endl;
  std::ofstream ofile;
  ofile.open(name.c_str());
  for(int row=0;row<ni;row++){
    for(int col=0;col<nj;col++){
      ofile << " " << A(row,col);
    }
    ofile << std::endl;
  }
  ofile << std::endl;
  ofile.close();

}
