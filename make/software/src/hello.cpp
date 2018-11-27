#include <string>
#include "hello.h"

void helloworld(){
  printf("Hello, world!\n");
}
void hello(std::string world){
  printf("Hello, %s\n",world.c_str());
}
