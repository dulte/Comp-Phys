#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


double funct(double x){
    return 100*exp(-10*x);
}

void gauss_elemination_general(int N, double *a, double* b, double* c, double* u, double* f) {
    for(int i = 2; i < N+1; i++){

        double d = c[i-1]/a[i-1];
        a[i] -= b[i-1]*d;
        f[i] -= f[i-1]*d;
    }

    u[N] = f[N]/a[N];

    for(int i = N-1; i > 0; i--){
        u[i] = (f[i] - b[i]*u[i+1])/a[i];
    }

}

int main(int argc, char *argv[])
{
    int N = 1000;
    double dt = 1.0/(N+1);

    double* a = new double[N+2];
    double* b = new double[N+2];
    double* c = new double[N+2];
    double* u = new double[N+2];
    double* f = new double[N+2];
    double* exactSol = new double[N+2];

    for (int i = 0; i  < N+1; i++){
        a[i] = 2;
        b[i] = -1;
        c[i] = -1;

        f[i] = dt*dt*funct(dt*i);
        exactSol[i] = 1 - (1-exp(-10))*(dt*i) - exp(-10*dt*i);

    }

    a[0] = 0;
    b[0] = 0;
    c[0] = 0;

    gauss_elemination_general(N,a,b,c,u,f);

    cout << u[5] << "  " << exactSol[5] << endl;
    return 0;
}
