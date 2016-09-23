all: gaussian gaussian_benchmark LUDecomp

gaussian:
	g++ gaussian.cpp -std=c++11 -O3 -o gauss.o

gaussian_benchmark:
		g++ gaussian_benchmark.cpp -std=c++11 -O3 -o gauss_benchmark.o

LUDecomp:
		g++ -g -std=c++11 lu_decomp.cpp -llapack -lblas -larmadillo -O3 -o lu.o
