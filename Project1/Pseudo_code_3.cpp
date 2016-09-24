#include <armadillo>

for (i = 0; i < n)
    A(i,i) = 2;
    if (i-1 >= 0):
        A(i,i-1) = -1;
    if (i+1 <= n-1):
        A(i,i+1) = -1;

arma::lu(L,U,A);
w = solve(L,f);
x = solve(U,w);
