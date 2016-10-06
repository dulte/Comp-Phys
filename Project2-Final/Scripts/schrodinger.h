#pragma once
#include <potential.h>
#include <armadillo>

using namespace arma;

class Schrodinger {
public:
    Schrodinger(Potential* potential);//Takes a Potential class as input, so we can use this class to make matrices for different potenials
    double setupMatrix(int n);

    mat setupMatrix(int n, int rho_max);
private:
    Potential* m_potential = nullptr;
};

