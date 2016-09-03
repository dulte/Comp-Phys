#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;


//Function decleration
double funct(double x);
void gauss_elemination_general(int N, double *a, double* b, double* c, double* u, double* f);
void gauss_elemination_spesific(int N, double *a_inv, double* u, double* f);
void writeArrayToFile(ofstream & outFile, double * array, int numBlocks);
double calculate_error(double,double);
void compute_solition(int N, double * u, double *exact);


//Main
int main(int argc, char *argv[])
{

    //Div variables
    int N = 10000;
    int iterations = 30;


    ofstream outFile_u("numericalSolution.bin");
    ofstream outFile_h("logerror.bin");



    double* max_errors = new double[iterations];

    for(int j = 0; j < iterations; j++){
        if (iterations > 1){
            N = pow(10,1+j/4.0);
            cout << N << endl;
        }
        else{
            N = 10000;
        }

        //Array that need updating depending on N
        double* u_numericalSolution = new double[N+2];
        double* exactSol = new double[N+2];
        double* error = new double[N+2];

        //Calculates the exact and numerical solution
        compute_solition(N,u_numericalSolution, exactSol);

        for(int k = 0; k < N+1; k++){
            error[k] = calculate_error(u_numericalSolution[k],exactSol[k]);
        }

        max_errors[j] = *max_element(error + 2, error + N-1);
        cout << max_errors[j] << endl;

        if (iterations > 1){
            delete [] u_numericalSolution;
            delete [] exactSol;
        }
        else{
            cout << u_numericalSolution[5] << "  " << exactSol[5] << endl;
            writeArrayToFile(outFile_u, u_numericalSolution, N+2);
        }
    }


    writeArrayToFile(outFile_h, max_errors, iterations);
    outFile_u.close();
    outFile_h.close();

    return 0;
}



//Function definitions:


//Function we want to integrate
double funct(double x){
    return 100*exp(-10*x);
}

void compute_solition(int N, double * u, double * exact){

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

//Writes arrays to Files
void writeArrayToFile(ofstream & outFile, double * array, int numBlocks)
{
    outFile.write(reinterpret_cast<char*>(array), numBlocks*sizeof(double));
}
