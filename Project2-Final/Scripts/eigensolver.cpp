#include "eigensolver.h"
#include <armadillo>
#include <cmath>

using namespace arma;
using namespace std;


eigensolver::eigensolver(mat& matrix_to_be_solved,int n)
{
    A = matrix_to_be_solved;
    col = n;
    R = eye<mat>(col,col);
}

void eigensolver::jacobi(double eps, bool test_enabled){

    vec test_vec1 = zeros(col);
    test_vec1(0) = 1;
    vec test_vec2 = zeros(col);
    test_vec2(1) = 1;


    int times_looped = 0;
    bool done = false;

    mat B = A;
    eigenvalues = zeros(col);

    while (!done){
        vec max_index = largest_non_diag(B);
        int k = max_index[0];
        int l = max_index[1];

        double t;
        double tau = (B(l,l)-B(k,k))/(2*B(k,l));

        if( tau >= 0 ){
           t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }else{
           t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }

        double c = 1.0/(sqrt(1+t*t));
        double s = c*t;



        double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
        a_kk = B(k,k);
        a_ll = B(l,l);
        B(k,k) = c*c*a_kk - 2.0*c*s*B(k,l) + s*s*a_ll;
        B(l,l) = s*s*a_kk + 2.0*c*s*B(k,l) + c*c*a_ll;
        B(k,l) = 0.0;
        B(l,k) = 0.0;
        for ( int i = 0; i < col; i++ ) {
          if ( i != k && i != l ) {
            a_ik = B(i,k);
            a_il = B(i,l);
            B(i,k) = c*a_ik - s*a_il;
            B(k,i) = B(i,k);
            B(i,l) = c*a_il + s*a_ik;
            B(l,i) = B(i,l);
          }
          r_ik = R(i,k);
          r_il = R(i,l);

          R(i,k) = c*r_ik - s*r_il;
          R(i,l) = c*r_il + s*r_ik;
        }


        if ((times_looped%10 == 0) && test_enabled){
            if(dot(B*test_vec1,B*test_vec2) != 0){
                cout << "Orthogonality was not conserved!" << endl;
                exit(1);
            }
        }


        times_looped++;
        done = (off_mat(B)<eps);
    }

    cout << times_looped << endl;

    eigenvalues = diagvec(B);
}

double eigensolver::off_mat(mat& M){
    double norm_of_mat = norm(M,"fro");
    double norm_of_diag = 0;

    for (int i = 0; i<col;i++){
        norm_of_diag += M(i,i)*M(i,i);
    }

    return norm_of_mat*norm_of_mat - norm_of_diag;
}

vec eigensolver::largest_non_diag(mat& M){

    vec max_el_index = zeros(2);
    double max_el = 0;


    for (int i = 0; i < col; i++){
        for (int j = 0; j < col; j++){
            if (j != i){
                if (abs(M(i,j)) > max_el){
                    max_el = abs(M(i,j));
                    max_el_index[0] = i;
                    max_el_index[1] = j;
                }
            }
        }
    }

    return max_el_index;
 }


vec eigensolver::get_eigenvalues(){
    return eigenvalues;
}

mat eigensolver::get_eingenvectors(){
    return R;
}


void eigensolver::test_max_value(double max_val, mat& M){
    vec calc_max_index = largest_non_diag(M);
    double calc_max = M(calc_max_index(0),calc_max_index(1));
    if (calc_max != max_val){
        cout << "The function did not calculate the correct largest element!" << endl;
        exit(1);
    }
}
