#ifndef __FILE_H__
#define __FILE_H__

#include<string>
#include<fstream>
#include<Eigen/Dense>

bool BothAreSpaces(char lhs, char rhs);
void RemoveSpaces(std::string& str);
void print_file(std::string filename);
void print_formatted_file(std::string filename,std::string delimiter);
void load_file(std::string filename, Eigen::MatrixXd &Matrix,std::string delimiter);
void sum_file(std::string filename,std::string delimiter);

#endif
