#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <armadillo>
#include <string>

using namespace std;

double exact_function(double x);
void lu_decompose(int n, string filename);

int main(int argc, char *argv[]){

    clock_t start, finish;
    double timeused;
    string filename;
    int i;
    int n;
    int exponent;

    if (argc > 1) {
        exponent = atoi(argv[1]);
    } else {
        cout << "Bad usage. Program " << argv[0] << " should read max "
            << "exponent of 10 for n" << endl;
        exit(1);
    }

    string argument;
    for (i = 1; i <= exponent; ++i) {
        filename = "data/luDecomp";
        argument = to_string(i);
        filename.append(to_string(i));
        filename.append(".txt");
        n = pow(10,i);
        if (i == 4) {
            cout << "warning, n = " << n << " may take long time"<< endl;
        }
        start = clock();
        lu_decompose(n, filename);
        finish = clock();
        timeused = (double) (start - finish) / CLOCKS_PER_SEC;
        cout << "Time used for lu decomp with n = "<< n <<"(seconds): " << timeused << endl;
    }

    return 0;
}
void lu_decompose(int n, string filename){
    using arma::mat;
    using arma::zeros;
    using arma::vec;
    using arma::solve;
    mat A = zeros(n,n);

    vec f = zeros(n);
    vec x = zeros(n);
    vec w = zeros(n);
    mat L,U;

    double h = 1.0/(n+1);
    int i;
    for (i = 0; i < n; ++i) {
        A(i,i) = 2;
        f(i) = h*h*exact_function(i*h);
        if (i-1 >= 0) {
            A(i,i-1) = -1;
        }
        if (i+1 <= n-1) {
            A(i,i+1) = -1;
        }
    }

    //arma::vec u = arma::solve(A,f);
    arma::lu(L,U,A);
    w = solve(L,f);
    x = solve(U,w);
    x.save(filename , arma::raw_ascii);
}

double exact_function(double x){
    return 100 * exp(-10*x);
}
