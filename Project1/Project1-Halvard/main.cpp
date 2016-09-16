#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <armadillo>
#include <string>

using namespace std;

double exact_sol(double x);
double function(double x);

void writeTextToFile(double * array, ofstream & outfile, int N);

void lu_decompose(int n);
void initialize(int n, double *a, double *b, double *c, double * u, double * f , double * f_tilde, double * exact);
void gaussian_general(int n, double *a, double *b, double *c, double * u, double * f , double * f_tilde, double * exact);


int main(int argc, char *argv[])
{
    string gaussian_general = ("general");
    string gaussian_specific = ("specific");

    int i;

    int n = 1000;
    double h = 1/(n+1.); // step length

    double *a = new double[n+2];
    double *b = new double[n+2];
    double *c= new double[n+2];
    double *f = new double[n+2];
    double *u = new double[n+2];
    double *f_tilde = new double[n+2];
    double *exact = new double[n+2];

    for (i = 0; i < n; ++i) {
        a[i] = 2;
        b[i] = -1;
        c[i] = -1;
        f[i] = h*h*function(i * h);
        exact[i] = exact_sol(i * h);
    }



    // Solving with the general gaussian elimination
    clock_t start, finish;

    if (chosen_method.compare(gaussian_general)!= 0) {
        cout << "Using general algorithm for solving matrix problem" << endl;
        start = clock();
        gaussian_general(n, a, b, c, u, f, f_tilde, exact);
    } else if (chosen_method.compare(gaussian_specific)!= 0){
        cout << "Using general algorithm for solving matrix problem" << endl;
        start = clock();
        gaussian_specific(n, a, u, f);
    } else {
        cout << chosen_method << " is not an implemented method" << endl;
        sys.exit();
    }
    finish = clock();
    double timeused = (double) (start - finish) / CLOCKS_PER_SEC;
    cout << "Time used (seconds): " << timeused << endl;



    ofstream outFile1("ComputedSolution.txt");
    ofstream outFile2("ExactSolution.txt");
    writeTextToFile(u,outFile1, n+1);
    writeTextToFile(exact,outFile2, n+1);
    

    start = clock();
    lu_decompose(n);
    finish = clock();
    timeused = (double) (start - finish) / CLOCKS_PER_SEC;
    cout << "Time used for lu decomp(seconds): " << timeused << endl;
    return 0;
}


void writeTextToFile(double * array, ofstream & outFile, int N) {
    for (int i = 0; i < N; ++i) {
        outFile << array[i] << endl;
        
    }
}

void lu_decompose(int n){
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
        f(i) = h*h*function(i*h);
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
    x.save("ArmaSolution.txt", arma::raw_ascii);
}

void gaussian_general(int n, double *a, double *b, double *c, double * u
                        , double * f , double * f_tilde, double * exact)
{
    for (int i = 2; i < n+1; ++i) {
        a[i] = a[i] - b[i-1] * c[i-1]/a[i-1];
        f_tilde[i] =  f[i] - f_tilde[i-1] * c[i-1] / a[i-1];
    }

    u[0] = 0;
    u[n+1] = 0;

    for (int i = n; i >= 1; --i) {
        u[i] = (f_tilde[i] - b[i]*u[i+1])/a[i];
        //cout << "i = " << i << "  u_i = " << u[i] << "   Exact:  " << exact[i] << endl;
    }
}
void gaussian_specific(int n, double *a, double* u, double * f){
    for (int i = 2; i < N+1 ; i++){
        f[i] += f[i-1]*a-inv[i-1];
    }
    u[n] = f[n]*a_inv[n];

    for (int i = n-1; i > 0; i--){
        u[i] = (f[i] + u[i+1])*a_inv[i];
    }
}

//The specific algorithm
void gauss_elemination_spesific(int N, double *a_inv, double* u, double* f){
    for(int i = 2; i < N+1; i++){
        f[i] += f[i-1]*a_inv[i-1];
    }
    u[N] = f[N]*a_inv[N];

    for(int i = N-1; i > 0; i--){
        u[i] = (f[i] + u[i+1])*a_inv[i];
    }
}
double exact_sol(double x){
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}
double function(double x){
    return 100 * exp(-10*x);
}
