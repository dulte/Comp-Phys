all: gaussian LUDecomp

gaussian:
	g++ gaussian.cpp -std=c++11 -O3 -o gauss.o
LUDecomp:
		g++ -g -std=c++11 lu_decomp.cpp -llapack -lblas -larmadillo -O3 -o lu.o
