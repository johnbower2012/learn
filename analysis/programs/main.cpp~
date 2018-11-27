#include "coshfunc.h"
#include "file.h"

int main(int argc, char* argv[]){
  int 
    ab=4,
    nmax=4,
    n=500,
    samples=1000;
  double 
    x=8;

  Eigen::VectorXd 
    width(ab);
  Eigen::MatrixXd 
    g(ab,nmax+1),
    G,
    M(n+1,ab*samples+1),
    m;

  std::string 
    infilename = "moments_parameters.dat",
    outfilename = "wave.dat";

  CDistCosh 
    dist;

  load_file(infilename, G, 24, 1000);
  dist.FunctionSet(n,x,samples,ab,nmax,G,M);
  print_file(outfilename,M);

  return 0;
}

