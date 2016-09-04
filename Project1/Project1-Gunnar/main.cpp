#include <iostream>
#include <cmath>

using namespace std;

double funct(double x){
    return 100*exp(10*x);
}

double exact_func(double x){
   return 1-(1-exp(-10))*x-exp(-10*x);
}

void gauss_elimination_general(int N, double *a, double *b, double *c, double *u, double *f){
    for (int i=0; i < N+2; i++){
        a[i]=a[i]-(b[i-1]*c[i-1])/(float(a[i-1]));
        f[i]=f[i]-f[i-1]*c[i-1]/float(a[i-1]);
        u[i]=(f[i]-b[i]u[i+1])/float(a[i]);
    }
}

int main()
{
    int N = 10000;
    double dt = 1./N;
    double *a = new double(N+2);
    double *b = new double(N+2);
    double *c = new double(N+2);
    double *u = new double(N+2);
    double *f = new double(N+2);
    double *exact_solution=new double[N+2]
    u[0]=0;
    u[N+1]=0;
    for (int i=0; i < N+2; i++){
        a[i]=-2;
        b[i]=1;
        c[i]=1;

        f[i]=funct(dt*i);
        exact_solution[i]=exact_func(dt*i)
    }

    return 0;

}
