#include "coloumb.h"

Coloumb::Coloumb(double omega_)
{
    omega = omega_;

}

double Coloumb::computePotential(double r){
    return omega*omega*r*r + 1.0/r;
}
