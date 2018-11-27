#include<iostream>
#include<string>
#include<Eigen/Dense>
#include "file.h"

int main(int argc, char* argv[]){
  std::string filename="numbers",
    delimiter=" ";
  Eigen::MatrixXd Matrix;
  printf("Print File:\n----------------\n");
  print_file(filename);
  printf("\nFormatted Print File:\n----------------\n");
  print_formatted_file(filename,delimiter);
  printf("\nSum File:\n----------------\n");
  sum_file(filename,delimiter);
  printf("\nData from File:\n----------------\n");
  load_file(filename,Matrix,delimiter);

  std::cout << Matrix << std::endl;

  return 0;  
}

