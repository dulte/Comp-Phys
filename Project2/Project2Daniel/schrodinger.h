#pragma once
#include <potential.h>
#include <armadillo>

using namespace arma;

class Schrodinger {
public:
    Schrodinger(Potential* potential);
    double setupMatrix(int n);

    mat setupMatrix(int n, int rho_max);
private:
    Potential* m_potential = nullptr;
};

