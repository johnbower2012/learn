## output file name
MAIN=main

## software directories
#### include: *.h  obj: *.o  lib: *.a  src: *.cpp
SOFT=../software
IDIR=$(SOFT)/include
ODIR=$(SOFT)/obj
LDIR=$(SOFT)/lib
SDIR=$(SOFT)/src

## compiler and compiler flags
CC=g++
CFLAGS=-I$(IDIR)

## additional libraries to link
LIBS=-larmadillo


## dependecies
#### pattern substition adds folder location to *.h files
_DEPS=hello.h
DEPS=$(patsubst %, $(IDIR)/%, $(_DEPS))

## object files
#### pattern substition adds folder location to *.o files
_OBJ=main.o hello.o
OBJ=$(patsubst %, $(ODIR)/%, $(_OBJ))