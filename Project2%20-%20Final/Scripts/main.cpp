#include <iostream>
#include <armadillo>
#include "eigensolver.h"
#include <schrodinger.h>
#include <harmonicoscillator.h>
#include <coloumb.h>
#include <string>

using namespace std;
using namespace arma;


int index_min(vec& V, int n);

int main(int argc, char *argv[])
{

    int n = 200;
    mat A_harmonic;
    mat A_coloumb;

    Schrodinger* s_harmonic = new Schrodinger(new HarmonicOscillator(5));
    Schrodinger* s_coloumb = new Schrodinger(new Coloumb(5));


    A_harmonic = s_harmonic->setupMatrix(n,2);
    A_coloumb = s_coloumb->setupMatrix(n,2);


    cout << "Found matrix" << endl;

    eigensolver* solver_harmonic = new eigensolver(A_harmonic,n);
    eigensolver* solver_coloumb = new eigensolver(A_coloumb,n);

    solver_harmonic->jacobi(0.000001);

    cout << "Done with the first calculation." << endl;

    solver_coloumb->jacobi(0.000001);


    vec eigen_harmonic = solver_harmonic->get_eigenvalues();
    mat eigen_vec_harmonic = solver_harmonic->get_eingenvectors();


    vec eigen_coloumb = solver_coloumb->get_eigenvalues();
    mat eigen_vec_coloumb = solver_coloumb->get_eingenvectors();

    //uword min_index = eigen.index_min();

    //vec func = eigen_vec.col(index_min(eigen,n));

    eigen_vec_harmonic.save("eigenvec_harmonic.txt",raw_ascii);
    eigen_harmonic.save("eigenval_harmonic.txt",raw_ascii);

    eigen_coloumb.save("eigenval_coloumb.txt",raw_ascii);
    eigen_vec_coloumb.save("eigenvec_coloumb.txt",raw_ascii);


    //sort(eigen_harmonic).print();
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
