#include "coshfunc.h"

#define PI 4.0*atan(1.0)
CDistribution::CDistribution(){
}
void CDistribution::GenX(int ab, double &x, double &weight){
}

/*
CDistCosh::CDistCosh(string parsfilename){
  int ab=4,ng=10;
  width.resize(ab);
  g.resize(ab,std::vector<double> (ng+1));
  parameter::ReadParsFromFile(parmap,parsfilename.c_str());

  width[0]=parameter::getD(parmap,"InvCoshPar_uu_width",1.0);
  width[1]=parameter::getD(parmap,"InvCoshPar_ud_width",1.0);
  width[2]=parameter::getD(parmap,"InvCoshPar_us_width",1.0);
  width[3]=parameter::getD(parmap,"InvCoshPar_ss_width",1.0);

  g[0][0]=parameter::getD(parmap,"InvCoshPar_uu_g0",1.0);
  g[0][1]=parameter::getD(parmap,"InvCoshPar_uu_g1",0.0);
  g[0][2]=parameter::getD(parmap,"InvCoshPar_uu_g2",0.0);
  g[0][3]=parameter::getD(parmap,"InvCoshPar_uu_g3",0.0);
  g[0][4]=parameter::getD(parmap,"InvCoshPar_uu_g4",0.0);
  g[0][5]=parameter::getD(parmap,"InvCoshPar_uu_g5",0.0);
  g[0][6]=parameter::getD(parmap,"InvCoshPar_uu_g6",0.0);
  g[0][7]=parameter::getD(parmap,"InvCoshPar_uu_g7",0.0);
  g[0][8]=parameter::getD(parmap,"InvCoshPar_uu_g8",0.0);
  g[0][9]=parameter::getD(parmap,"InvCoshPar_uu_g9",0.0);
  g[0][10]=parameter::getD(parmap,"InvCoshPar_uu_g10",0.0);

  g[1][0]=parameter::getD(parmap,"InvCoshPar_ud_g0",1.0);
  g[1][1]=parameter::getD(parmap,"InvCoshPar_ud_g1",0.0);
  g[1][2]=parameter::getD(parmap,"InvCoshPar_ud_g2",0.0);
  g[1][3]=parameter::getD(parmap,"InvCoshPar_ud_g3",0.0);
  g[1][4]=parameter::getD(parmap,"InvCoshPar_ud_g4",0.0);
  g[1][5]=parameter::getD(parmap,"InvCoshPar_ud_g5",0.0);
  g[1][6]=parameter::getD(parmap,"InvCoshPar_ud_g6",0.0);
  g[1][7]=parameter::getD(parmap,"InvCoshPar_ud_g7",0.0);
  g[1][8]=parameter::getD(parmap,"InvCoshPar_ud_g8",0.0);
  g[1][9]=parameter::getD(parmap,"InvCoshPar_ud_g9",0.0);
  g[1][10]=parameter::getD(parmap,"InvCoshPar_ud_g10",0.0);

  g[2][0]=parameter::getD(parmap,"InvCoshPar_us_g0",1.0);
  g[2][1]=parameter::getD(parmap,"InvCoshPar_us_g1",0.0);
  g[2][2]=parameter::getD(parmap,"InvCoshPar_us_g2",0.0);
  g[2][3]=parameter::getD(parmap,"InvCoshPar_us_g3",0.0);
  g[2][4]=parameter::getD(parmap,"InvCoshPar_us_g4",0.0);

  g[2][5]=parameter::getD(parmap,"InvCoshPar_us_g5",0.0);
  g[2][6]=parameter::getD(parmap,"InvCoshPar_us_g6",0.0);
  g[2][7]=parameter::getD(parmap,"InvCoshPar_us_g7",0.0);
  g[2][8]=parameter::getD(parmap,"InvCoshPar_us_g8",0.0);
  g[2][9]=parameter::getD(parmap,"InvCoshPar_us_g9",0.0);
  g[2][10]=parameter::getD(parmap,"InvCoshPar_us_g10",0.0);

  g[3][0]=parameter::getD(parmap,"InvCoshPar_ss_g0",1.0);
  g[3][1]=parameter::getD(parmap,"InvCoshPar_ss_g1",0.0);
  g[3][2]=parameter::getD(parmap,"InvCoshPar_ss_g2",0.0);
  g[3][3]=parameter::getD(parmap,"InvCoshPar_ss_g3",0.0);
  g[3][4]=parameter::getD(parmap,"InvCoshPar_ss_g4",0.0);
  g[3][5]=parameter::getD(parmap,"InvCoshPar_ss_g5",0.0);
  g[3][6]=parameter::getD(parmap,"InvCoshPar_ss_g6",0.0);
  g[3][7]=parameter::getD(parmap,"InvCoshPar_ss_g7",0.0);
  g[3][8]=parameter::getD(parmap,"InvCoshPar_ss_g8",0.0);
  g[3][9]=parameter::getD(parmap,"InvCoshPar_ss_g9",0.0);
  g[3][10]=parameter::getD(parmap,"InvCoshPar_ss_g10",0.0);

  coshcalc=new CCosh(ng);
}
*/
CDistCosh::CDistCosh(std::vector<double> Width, std::vector<std::vector<double> > G) : coshcalc(){
  int ab=Width.size(),ng=G[0].size();
  this->width = Width;
  this->g = G;
  
  coshcalc.Set_NMax(ng);
}
void CDistCosh::Set_WG(std::vector<double> Width, std::vector<std::vector<double> > G){
  this->width = Width;
  this->g = G;
  int ng=G[0].size();
  coshcalc.Set_NMax(ng);
}
void CDistCosh::GenX(int ab, double &x, double &weight){
  int n,ng=coshcalc.nmax;
  x=coshcalc.GetRandomX();
  double Z=0.0;
  weight=0.0;
  for(n=0;n<ng+1;n++){
    weight+=coshcalc.GetG(n,x)*g[ab][n];
    Z += coshcalc.Z(n)*g[ab][n];
  }
  //weight /= Z;
  weight=weight/(coshcalc->GetG(0,x)*g[ab][0]);
  x*=width[ab];
}
void CDistCosh::Functions(int n, double x, Eigen::MatrixXd &M){
  int ab = width.size();
  int nmax = coshcalc.nmax;
  double dx=x/double(n), xx=0.0;
  M = Eigen::MatrixXd::Zero(n+1,ab+1);
  for(int xi=0;xi<n+1;xi++){
    M(xi,0) = xx;
    for(int species=0;species<ab;species++){
      for(int m=0;m<nmax+1;m++){
	M(xi,species+1) += coshcalc.GetG(m,xx/width[species])*g[species][m];
      }
    }
    xx += dx;
  }
}
double CDistCosh::Integrate(int ab){
  int m;
  int nmax = coshcalc.nmax;
  double sum=0.0;
  for(m=0;m<nmax+1;m++){
    sum += coshcalc.Z(m)*g[ab][m];
  }
  sum *= width[ab];
  return sum;
}
double CDistCosh::IntegrateNumerical(int ab){
  double
    x=100.0,
    xx=0.0,
    dx=0.0001,
    sum=0.0;
  int
    m;
  int nmax = coshcalc.nmax;
  for(m=0;m<nmax+1;m++){
    for(xx=0.0;xx<x;xx+=dx){
      sum += this->coshcalc.GetG(m,xx)*g[ab][m];
    }
  }
  sum *= dx;
  sum *= width[ab];
  return sum;
}
CCosh::CCosh(){
}
CCosh::CCosh(int nmaxset){
  this->nmax=nmaxset;
  Eigen::VectorXd x,y,ytest;
  int n,nprime,m,mprime;
  double norm,overlap;
  I.resize(2*nmax+2);
  I(0)=0.5*PI;
  I(1)=1.0;

  for(n=2;n<=2*nmax+1;n++){
    I(n)=I(n-2)*double(n-1)/double(n);
  }
  
  A = Eigen::MatrixXd::Zero(nmax+1,nmax+1);
  Z = Eigen::VectorXd::Zero(nmax+1);
  Eigen::MatrixXd K;

  A(0,0)=1.0;
  Z(0) = A(0,0)*I(0);
  for(n=1;n<=nmax;n++){
    A(n,n)=1.0;
    K = Eigen::MatrixXd::Zero(n,n);
    x = Eigen::VectorXd::Zero(n);
    y = Eigen::VectorXd::Zero(n);
    ytest = Eigen::VectorXd::Zero(n);
    for(nprime=0;nprime<n;nprime++){
      for(m=0;m<n;m++){
	K(nprime,m)=0.0;
	for(mprime=0;mprime<=nprime;mprime++){
	  K(nprime,m)+=I(m+mprime+1)*A(nprime,mprime);
	}
      }
    }
    for(nprime=0;nprime<n;nprime++){
      y(nprime)=0.0;
      for(mprime=0;mprime<=nprime;mprime++)
	y(nprime)+=I(n+mprime+1)*A(nprime,mprime);
    }
    x = K.colPivHouseholderQr().solve(y);
    for(m=0;m<n;m++)
      A(n,m)=-x(m);
    norm=0.0;
    for(m=0;m<=n;m++){
      for(mprime=0;mprime<=n;mprime++){
	norm+=A(n,m)*A(n,mprime)*I(m+mprime+1);
      }
    }
    if(norm<0.0){
      printf("norm(%d)=%g\n",n,norm);
      printf("negative norm, try making nmax smaller\n");
      exit(1);
    }
    x=x/sqrt(norm);
    y=y/sqrt(norm);
    //ytest=K*x-y;                                                                                                              
    //ytest.print("ytest:");                                                                                                    
    for(m=0;m<=n;m++){
      A(n,m)=A(n,m)/sqrt(norm);
      Z(n) += A(n,m)*I(m);
    }
   }

  seed = 123456;
  rng = std::mt19937(seed);
  uniform_dist = std::uniform_real_distribution<double>(0,1);
}
void CCosh::Set_NMax(int nmaxset){
  this->nmax=nmaxset;
  Eigen::VectorXd x,y,ytest;
  int n,nprime,m,mprime;
  double norm,overlap;
  I.resize(2*nmax+2);
  I(0)=0.5*PI;
  I(1)=1.0;

  for(n=2;n<=2*nmax+1;n++){
    I(n)=I(n-2)*double(n-1)/double(n);
  }
  
  A = Eigen::MatrixXd::Zero(nmax+1,nmax+1);
  Z = Eigen::VectorXd::Zero(nmax+1);
  Eigen::MatrixXd K;

  A(0,0)=1.0;
  Z(0) = A(0,0)*I(0);
  for(n=1;n<=nmax;n++){
    A(n,n)=1.0;
    K = Eigen::MatrixXd::Zero(n,n);
    x = Eigen::VectorXd::Zero(n);
    y = Eigen::VectorXd::Zero(n);
    ytest = Eigen::VectorXd::Zero(n);
    for(nprime=0;nprime<n;nprime++){
      for(m=0;m<n;m++){
	K(nprime,m)=0.0;
	for(mprime=0;mprime<=nprime;mprime++){
	  K(nprime,m)+=I(m+mprime+1)*A(nprime,mprime);
	}
      }
    }
    for(nprime=0;nprime<n;nprime++){
      y(nprime)=0.0;
      for(mprime=0;mprime<=nprime;mprime++)
	y(nprime)+=I(n+mprime+1)*A(nprime,mprime);
    }
    x = K.colPivHouseholderQr().solve(y);
    for(m=0;m<n;m++)
      A(n,m)=-x(m);
    norm=0.0;
    for(m=0;m<=n;m++){
      for(mprime=0;mprime<=n;mprime++){
	norm+=A(n,m)*A(n,mprime)*I(m+mprime+1);
      }
    }
    if(norm<0.0){
      printf("norm(%d)=%g\n",n,norm);
      printf("negative norm, try making nmax smaller\n");
      exit(1);
    }
    x=x/sqrt(norm);
    y=y/sqrt(norm);
    //ytest=K*x-y;                                                                                                              
    //ytest.print("ytest:");                                                                                                    
    for(m=0;m<=n;m++){
      A(n,m)=A(n,m)/sqrt(norm);
      Z(n) += A(n,m)*I(m);
    }
   }

  seed = 123456;
  rng = std::mt19937(seed);
  uniform_dist = std::uniform_real_distribution<double>(0,1);
}
double CCosh::GetG(int n,double x){
  double g,coshx;
  int m;
  coshx=cosh(x);
  g=0.0;
  for(m=0;m<n+1;m++){
    g+=this->A(n,m)*pow(coshx,-m-1);
  }
  return g;
}
double CCosh::GetRandomX(){
  //random number from (2/pi)/cosh(x)
  return asinh(tan(0.5*PI*uniform_dist(rng)));
}
void CCosh::TestOverlap(){
  int n,nprime,nmax=this->nmax;
  for(n=0;n<nmax+1;n++){
    for(nprime=0;nprime<nmax+1;nprime++){
      this->CalcOverlap(n,nprime);
      this->CalcOverlapNumerical(n,nprime);
      printf("----------------------\n");
    }
  }
}
void CCosh::CalcOverlap(int n,int nprime){
  double overlap=0.0;
  int m,mprime;
  for(m=0;m<=n;m++){
    for(mprime=0;mprime<=nprime;mprime++){
      overlap+=A(n,m)*I(m+mprime+1)*A(nprime,mprime);
    }
  }
  printf("overlap(%d,%d)=%g\n",n,nprime,overlap);
}
void CCosh::CalcOverlapNumerical(int n,int nprime){
  double overlap=0.0,g,gprime;
  double x,dx=0.0002,xmax=35.0;
  for(x=0.5*dx;x<xmax;x+=dx){
    g=this->GetG(n,x);
    gprime=this->GetG(nprime,x);
    overlap+=g*gprime*dx;
  }
  if(abs(overlap)<5.0E-9)
    overlap=0.0;
  printf("numerical: overlap(%d,%d)=%g\n",n,nprime,overlap);
}
void CCosh::Integrate(double x){
  printf("Testing getting random no. from dist (2/pi)/cosh(x)\n");
  double answer=0,dx=0.0001,xx;
  for(xx=0.5*dx;xx<x;xx+=dx){
    answer+=dx/cosh(xx);
  }
  printf("%g =? %g\n",answer,atan(sinh(x)));

  double y,y2bar,dy=0.01;
  int ibin,nbins=1000;
  double *prob=new double[nbins];
  int n,nmax=1000000;
  for(n=0;n<nmax;n++){
    y=GetRandomX();
    ibin=lrint(floor(y/dy));
    if(ibin<nbins)
      prob[ibin]+=1.0/(nmax*dy);
  }
  for(ibin=0;ibin<nbins;ibin++){
    y=(ibin+0.5)*dy;
    printf("%5d %6.3f %10.7f =? %10.7f\n",ibin,y,prob[ibin],(2.0/PI)/cosh(y));
  }
}
