#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <armadillo>

using namespace arma;
using namespace std;

class eigensolver
{
public:
    eigensolver(mat& matrix_to_be_solved,int n);

    void jacobi(double eps, bool test_enabled);
    void jacobi(double eps) {jacobi(eps, false);}
    vec get_eigenvalues();
    mat get_eingenvectors();
    void test_max_value(double max_val, mat& M);


private:
    int col;
    mat A;
    mat R;
    vec eigenvalues;


    double off_mat(mat& M);
    vec largest_non_diag(mat& M);
};

#endif // EIGENSOLVER_H
