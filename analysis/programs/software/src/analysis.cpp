#include "analysis.h"

void RemoveRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();

  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

  matrix.conservativeResize(numRows,numCols);
}

void RemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;

  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

  matrix.conservativeResize(numRows,numCols);
}
void AverageColumns(Eigen::VectorXd &average, Eigen::MatrixXd matrix){
  int rows = matrix.rows(),
    cols = matrix.cols();
  average = Eigen::VectorXd::Zero(cols);
  for(int col=0;col<cols;col++){
    for(int row=0;row<rows;row++){
      average(col) += matrix(row,col);
    }
    average(col) /= rows;
  }
}
void AddOnesRow(Eigen::MatrixXd matrix, Eigen::MatrixXd &outMatrix){
  int rows = matrix.rows(),
    cols = matrix.cols();
  outMatrix = Eigen::MatrixXd::Zero(rows+1,cols);
  for(int col=0;col<cols;col++){
    outMatrix(0,col) = 1.0;
    for(int row=0;row<rows;row++){
      outMatrix(row+1,col) = matrix(row,col);
    }
  }
}
void AddOnesColumn(Eigen::MatrixXd matrix, Eigen::MatrixXd &outMatrix){
  int rows = matrix.rows(),
    cols = matrix.cols();
  outMatrix = Eigen::MatrixXd::Zero(rows,cols+1);
  for(int row=0;row<rows;row++){
    outMatrix(row,0) = 1.0;
    for(int col=0;col<cols;col++){
      outMatrix(row,col+1) = matrix(row,col);
    }
  }
}


void ScaleMatrixColumns(Eigen::MatrixXd Matrix,Eigen::VectorXd &Mean, Eigen::VectorXd &Std, Eigen::MatrixXd &Scaled)
{
  int rows = Matrix.rows(),
    cols = Matrix.cols();
  Mean = Eigen::VectorXd::Zero(cols);
  Std = Eigen::VectorXd::Zero(cols);
  Scaled = Eigen::MatrixXd::Zero(rows,cols);
  for(int col=0;col<cols;col++)
    {
      for(int row=0;row<rows;row++)
	{
	  Mean(col) += Matrix(row,col);
	}
      Mean(col) /= (double) rows;
      for(int row=0;row<rows;row++)
	{
	  Std(col) += (Matrix(row,col) - Mean(col))*(Matrix(row,col) - Mean(col));
	}
      Std(col) /= (double)rows;
      Std(col) = sqrt(Std(col));
      for(int row=0;row<rows;row++)
	{
	  Scaled(row,col) = (Matrix(row,col) - Mean(col))/Std(col);
	}
    }
}
void ScaleMatrixRows(Eigen::MatrixXd Matrix,Eigen::VectorXd &Mean, Eigen::VectorXd &Std, Eigen::MatrixXd &Scaled)
{
  int rows = Matrix.rows(),
    cols = Matrix.cols();
  Mean = Eigen::VectorXd::Zero(rows);
  Std = Eigen::VectorXd::Zero(rows);
  Scaled = Eigen::MatrixXd::Zero(rows,cols);
  for(int row=0;row<rows;row++)
    {
      for(int col=0;col<cols;col++)
	{
	  Mean(row) += Matrix(row,col);
	}
      Mean(row) /= (double) cols;
      for(int col=0;col<cols;col++)
	{
	  Std(row) += (Matrix(row,col) - Mean(row))*(Matrix(row,col) - Mean(row));
	}
      Std(row) /= (double)cols;
      Std(row) = sqrt(Std(row));
      for(int col=0;col<cols;col++)
	{
	  Scaled(row,col) = (Matrix(row,col) - Mean(row))/Std(row);
	}
    }
}
void AverageRows(Eigen::VectorXd &average, Eigen::MatrixXd matrix){
  int rows = matrix.rows(),
    cols = matrix.cols();
  average = Eigen::VectorXd::Zero(rows);
  for(int row=0;row<rows;row++){
    for(int col=0;col<cols;col++){
      average(row) += matrix(row,col);
    }
    average(row) /= cols;
  }
}
void TildeFunction(Eigen::MatrixXd &tilde, Eigen::VectorXd mean, Eigen::VectorXd error, Eigen::MatrixXd matrix){
  int rows = matrix.rows(),
    cols = matrix.cols();
  tilde = Eigen::MatrixXd::Zero(rows,cols);
  for(int row=0;row<rows;row++){
    for(int col=0;col<cols;col++){
      tilde(row,col) = (matrix(row,col) - mean(row))/error(row);
    }
  }
}
void CovarianceFunction(Eigen::MatrixXd &cov, Eigen::MatrixXd matrix){
  int rows=matrix.rows(),
    cols=matrix.cols();
  cov = Eigen::MatrixXd::Zero(rows,rows);
  for(int row1=0;row1<rows;row1++){
    for(int row2=0;row2<rows;row2++){
      for(int col=0;col<cols;col++){
	cov(row1,row2) += matrix(row1,col)*matrix(row2,col);
      }
      cov(row1,row2) /= (double) cols;
    }
  }
}
void EigenSolve(Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors, Eigen::MatrixXd matrix){
  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(matrix);
  if(eigensolver.info() != Eigen::Success) abort();
  eigenvalues = eigensolver.eigenvalues().real();
  eigenvectors = eigensolver.eigenvectors().real();
}
void FindMax(Eigen::VectorXd vector, int &index, double &max){
  int elems = vector.size();
  index=0;
  max=vector(0);
  for(int elem=0;elem<elems;elem++){
    if(vector(elem) > max){
      max = vector(elem);
      index = elem;
    }
  }
}
void FindMax(Eigen::MatrixXd vector, int &index, double &max){
  int elems = vector.rows();
  index=0;
  max=vector(0,0);
  for(int elem=0;elem<elems;elem++){
    if(vector(elem,0) > max){
      max = vector(elem,0);
      index = elem;
    }
  }
}
void EigenSort(Eigen::VectorXd &eigval, Eigen::MatrixXd &eigvec){
  int index,
    vals = eigval.size();
  double max;
  bool tick=false;
  Eigen::MatrixXd eigsort = Eigen::MatrixXd::Zero(vals,vals),
    eigval_matrix = Eigen::MatrixXd::Zero(vals,2);
  for(int val=0;val<vals;val++){
    eigval_matrix(val,0) = 0;
    eigval_matrix(val,1) = -1;
  } 
  for(int val=0;val<vals;val++){
    index=0;
    for(int search2=0;search2<val;search2++){
      if(eigval_matrix(search2,1)==index){
	index++;
	search2=-1;
      }
    }
    max = eigval(index,0);
    for(int search1=index;search1<vals;search1++){
      tick=false;
      for(int search2=0;search2<val;search2++){
	if(eigval_matrix(search2,1)==search1){
	  tick=true;
	  break;
	}
      }
      if(tick==true){
	continue;
      }
      if(eigval(search1)>max){
	max = eigval(search1);
	index = search1;
      }
      eigval_matrix(val,1) = index;
      eigval_matrix(val,0) = max;
    }
  }
  for(int col=0;col<vals;col++){
    eigval(col) = eigval_matrix(col,0);
    index=eigval_matrix(col,1);
    for(int row=0;row<vals;row++){
      eigsort(row,col) = eigvec(row,index);
    }
  }
  eigvec = eigsort;
}


void ZerothMoment(Eigen::MatrixXd Function, double &Moment){
  int points=Function.rows()-1;
  double f=0,dx=0;
  Moment=0;
  for(int point=0;point<points;point++){
    f = (Function(point+1,1) + Function(point,1))/2.0;
    dx = Function(point+1,0) - Function(point,0);
    Moment += f*dx;
  }
}
void FirstMoment(Eigen::MatrixXd Function, double &Moment){
  int points=Function.rows()-1;
  double f=0,x=0,dx=0,
    zero=0;
  Moment=0;
  for(int point=0;point<points;point++){
    f = (Function(point+1,1) + Function(point,1))/2.0;
    x = (Function(point+1,0) + Function(point,0))/2.0;    
    dx = Function(point+1,0) - Function(point,0);
    Moment += f*x*dx;
  }
  ZerothMoment(Function,zero);
  //  Moment /= zero;
}
void SecondMoment(Eigen::MatrixXd Function, double &Moment){
  int points=Function.rows()-1;
  double f=0,x=0,dx=0,
    zero=0, first=0;
  ZerothMoment(Function,zero);
  FirstMoment(Function,first);
  Moment=0;
  for(int point=0;point<points;point++){
    f = (Function(point+1,1) + Function(point,1))/2.0;
    x = (Function(point+1,0) + Function(point,0))/2.0;    
    dx = Function(point+1,0) - Function(point,0);
    Moment += f*(x-first)*(x-first)*dx;
  }
  //  Moment /= zero;
}
void MatrixMoments(std::vector<Eigen::MatrixXd> matrix, Eigen::VectorXd DelX, Eigen::MatrixXd &Obs){
  int files = matrix.size(),
    points = matrix[0].rows(),
    obs_file = 3,
    obs = obs_file*matrix.size(),
    runs = matrix[0].cols();
  double zero, first, second;
  Eigen::MatrixXd function = Eigen::MatrixXd::Zero(points,2);
  Obs = Eigen::MatrixXd::Zero(obs,runs);
  for(int file=0;file<files;file++){
    for(int run=0;run<runs;run++){
      for(int point=0;point<points;point++){
	function(point,0) = DelX(point);
	function(point,1) = matrix[file](point,run);
      }
      ZerothMoment(function,zero);
      FirstMoment(function,first);
      SecondMoment(function,second);
      Obs(file*obs_file, run) = zero;
      Obs(file*obs_file +1, run) = first;
      Obs(file*obs_file +2, run) = second;
    }
  }    
}
void linearRegressionLeastSquares(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd &Beta){
  int points = Y.rows();
  int parameters = X.cols();
  int observables = Y.cols();
  Beta = Eigen::MatrixXd::Zero(observables,parameters+1); 

  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(parameters+1,parameters+1),
    X_;
  AddOnesColumn(X,X_);
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(parameters);
  Eigen::VectorXd y = Eigen::VectorXd::Zero(points);

  temp = X_.transpose()*X_;
  temp = temp.inverse();
  for(int j=0;j<observables;j++){
    y = Y.col(j);
    beta = temp*X_.transpose()*y;
    Beta.row(j) = beta;
  }
}
void ComputeFit(Eigen::MatrixXd ModelZ, Eigen::MatrixXd EmulatorZ, Eigen::MatrixXd &outMatrix){
  int samples = ModelZ.rows(),
    obs = ModelZ.cols();

  outMatrix = Eigen::MatrixXd::Zero(samples,obs+1);

  double sum = 0.0;
  for(int i=0;i<samples;i++){
    sum = 0.0;
    for(int j=0;j<obs;j++){
      outMatrix(i,j+1) = ModelZ(i,j) - EmulatorZ(i,j);
      sum += outMatrix(i,j+1)*outMatrix(i,j+1);
    }
    outMatrix(i,0) = -0.5*sum;
  }
}
  
