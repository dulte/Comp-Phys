#include <iostream>
#include <armadillo>
#include "eigensolver.h"
#include <schrodinger.h>
#include <harmonicoscillator.h>

using namespace std;
using namespace arma;


int index_min(vec& V, int n);

int main(int argc, char *argv[])
{

    int n = 100;
    mat A;

    Schrodinger* s = new Schrodinger(new HarmonicOscillator());

    A = s->setupMatrix(n,10);


    cout << "Found matrix" << endl;

    eigensolver* solver = new eigensolver(A,n);

    solver->jacobi(0.000001);
    vec eigen = solver->get_eigenvalues();
    mat eigen_vec = solver->get_eingenvectors();

    //uword min_index = eigen.index_min();

    vec func = eigen_vec.col(index_min(eigen,n));

    func.save("eigenvec.txt",raw_ascii);



    sort(eigen).print();
    cout << "Hello World!" << endl;
    return 0;
}


int index_min(vec& V, int n){
    int min_i = 0;
    double min_val = 1000;

    for (int i = 0; i < n; i++){
        if(abs(V(i)) < min_val){
            min_val = abs(V(i));
            min_i = i;
        }
    }

    return min_i;
}
