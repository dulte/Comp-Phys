#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <armadillo>

using namespace std;


//Function decleration
double funct(double x);
void gauss_elemination_general(int N, double *a, double* b, double* c, double* u, double* f);
void gauss_elemination_spesific(int N, double *a_inv, double* u, double* f);
void writeArrayToFile(ofstream & outFile, double * array, int numBlocks);
double calculate_error(double,double);
void compute_solution(int N, double * u, double *exact);
void LUDecomp(int N);


//Main
int main(int argc, char *argv[])
{

    //Div variables
    int N = 3000;
    int iterations = 1; // 25; //If this variable is set to 1, you will get the numerical approx. with the N above, else you wil get the log error with different N's



    //Open files to write results
    ofstream outFile_u("numericalSolution.bin");
    ofstream outFile_error("logerror.bin");
    ofstream outFile_h("logh.bin");

    LUDecomp(3002);

    cout << "hei" << endl;



    double* max_errors = new double[iterations];
    double* log_h = new double[iterations];

    for(int j = 0; j < iterations; j++){
        if (iterations > 1){
            N = pow(10,1+j/4.0);
            log_h[j] = -(1+j/4.0);
            cout << N << endl;
        }

        //Array that need updating depending on N
        double* u_numericalSolution = new double[N+2];
        double* exactSol = new double[N+2];
        double* error = new double[N+2];

        //Calculates the exact and numerical solution
        compute_solution(N,u_numericalSolution, exactSol);

        for(int k = 0; k < N+1; k++){
            error[k] = calculate_error(u_numericalSolution[k],exactSol[k]);
        }

        max_errors[j] = *max_element(error + 2, error + N-1);
        cout << max_errors[j] << endl;

        if (iterations > 1){
            delete [] u_numericalSolution;
            delete [] exactSol;
            delete [] error;
        }
        else{
            cout << u_numericalSolution[5] << "  " << exactSol[5] << endl;
            writeArrayToFile(outFile_u, u_numericalSolution, N+2);
        }
    }

    if (iterations > 1){
        writeArrayToFile(outFile_error, max_errors, iterations);
        writeArrayToFile(outFile_h,log_h,iterations);
    }
    outFile_u.close();
    outFile_h.close();
    outFile_error.close();

    return 0;
}



//Function definitions:


//Function we want to integrate
double funct(double x){
    return 100*exp(-10*x);
}

void compute_solution(int N, double * u, double * exact){

    double dt = 1.0/(N+1);

    double* a = new double[N+2];
    double* b = new double[N+2];
    double* c = new double[N+2];

    double* f = new double[N+2];

    double* a_special_inv = new double[N+2]; //A special version calculated for the algoritem for this special kind of matrix

    //Loop to fill the array needed to solve the problem
    for (int i = 0; i  < N+1; i++){
        a[i] = 2;
        b[i] = -1;
        c[i] = -1;

        f[i] = dt*dt*funct(dt*i);
        exact[i] = 1 - (1-exp(-10))*(dt*i) - exp(-10*dt*i);

        a_special_inv[i] = (double)i/(i+1);

    }

    //gauss_elemination_general(N,a,b,c,u_numericalSolution,f);
    gauss_elemination_spesific(N,a_special_inv,u,f);

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] f;
    delete [] a_special_inv;
}

//The general algorithm
void gauss_elemination_general(int N, double *a, double* b, double* c, double* u, double* f) {

    //Forward loop
    for(int i = 2; i < N+1; i++){

        double d = c[i-1]/a[i-1];
        a[i] -= b[i-1]*d;
        f[i] -= f[i-1]*d;
    }

    u[N] = f[N]/a[N];

    //Backward loop
    for(int i = N-1; i > 0; i--){
        u[i] = (f[i] - b[i]*u[i+1])/a[i];
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

double calculate_error(double computed,double exact){
    if (computed == 0  || exact == 0){
        return 0.0;
    }
    else{
        return log10(abs((computed-exact)/exact));
    }
}

void LUDecomp(int N){
    arma::mat A = arma::zeros(N,N);
    arma::mat L,U;
    arma::vec w,x;

    arma::vec func = arma::zeros(N);

    double dt = 1.0/(N+1);

    for (int i = 0; i < N; i++){
        func(i) = dt*dt*funct(dt*i);
        A(i,i) = 2;
        if (i-1 >= 0){
            A(i,i-1) = -1;
        }
        if (i+1 <= N-1){
            A(i,i+1) = -1;
        }
    }

    arma::lu(L,U,A);
    w = arma::solve(L,func);
    x = arma::solve(U,w);
    x.save("num.txt", arma::raw_ascii);
}

//Writes arrays to Files
void writeArrayToFile(ofstream & outFile, double * array, int numBlocks)
{
    outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(double));
}
