#include "schrodinger.h"
#include <armadillo>

using namespace std;
using namespace arma;

Schrodinger::Schrodinger(Potential *potential) {
    m_potential = potential;

}

 mat Schrodinger::setupMatrix(int n,int rho_max) {


     int rho_min = 0;
     double h = double(rho_max-rho_min)/double(n);


     mat A = zeros(n,n);

     for (int i = 0; i < n; i++){
         A(i,i) = (2.0/(h*h)) + m_potential->computePotential(rho_min+h*(i+1));
         if(i<n-1){
             A(i,i+1) = -1.0/(h*h);
         }
         if(i>0){
             A(i,i-1) = -1.0/(h*h);
         }
     }

     return A;
}
