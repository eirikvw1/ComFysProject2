CC= g++ -std=c++11 -Wall -Wextra -Wpedantic -larmadillo -lblas -llapack
CFLAGS=-I.
DEPS = func.h


all: pro2_b 

pro2_b: pro2_b.cpp
	${CC} -o pro2.exe pro2_b.cpp  -larmadillo -llapack -lblas

