#include "file.h"

bool BothAreSpaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }
void RemoveSpaces(std::string& str){
  std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
  str.erase(new_end, str.end()); 
}
void print_file(std::string filename){
  std::string line;
  std::fstream myfile(filename,std::ios::in);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      printf("%s\n",line.c_str());
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}
void print_formatted_file(std::string filename,std::string delimiter){
  std::string line;
  std::fstream myfile(filename,std::ios::in);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      if(line[0]!='#' && line!="\0"){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	printf("%s\n",line.c_str());
      }
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}
void load_file(std::string filename, Eigen::MatrixXd &Matrix,std::string delimiter){
  std::string line,
    number;
  std::size_t pos,length;
  std::ifstream myfile(filename);
  int rows=0,row=0,
    cols=0,col=0,
    max=0;

  if(myfile.is_open()){
    while(getline(myfile,line)){
      cols=0;
      if( line[0]!='#' && line!="\0" ){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	while( (pos=line.find(delimiter)) != std::string::npos){
	  line.erase(0,pos+delimiter.length());
	  cols++;
	}
	if(cols>max){
	  max=cols;
	}
	rows++;
      }
      cols=max;
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }

  Matrix = Eigen::MatrixXd::Zero(rows,cols);
  myfile.open(filename);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      if(line[0]!='#' && line!="\0"){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	col=0;
	while( (pos=line.find(delimiter)) != std::string::npos){
	  number=line.substr(0,pos+delimiter.length());
	  Matrix(row,col) = std::stof( number );
	  line.erase(0,pos+delimiter.length());
	  col++;
	}
	row++;
      }
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }

}
void sum_file(std::string filename,std::string delimiter){
  std::string line,
    number;
  std::size_t pos,length;
  std::ifstream myfile(filename);
  double sum;
  if(myfile.is_open()){
    while(getline(myfile,line)){
      if(line[0]!='#' && line!="\0"){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	sum=0;
	while( (pos=line.find(delimiter)) != std::string::npos){
	  number=line.substr(0,pos+delimiter.length());
	  sum += std::stof( number );
	  printf("%s",number.c_str());
	  line.erase(0,pos+delimiter.length());
	}
	printf("\n");
	printf("sum: %f\n",sum);
      }
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}
