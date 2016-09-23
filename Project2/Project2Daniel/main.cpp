#include <iostream>
#include <armadillo>
#include "eigensolver.h"
#include <schrodinger.h>
#include <harmonicoscillator.h>

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{

    int n = 200;
    mat A;

    Schrodinger* s = new Schrodinger(new HarmonicOscillator());

    A = s->setupMatrix(n,10);


    cout << "Found matrix" << endl;

    eigensolver* solver = new eigensolver(A,n);

    solver->jacobi(0.000001);
    vec eigen = solver->get_eigenvalues();



    sort(eigen).print();
    cout << "Hello World!" << endl;
    return 0;
}
