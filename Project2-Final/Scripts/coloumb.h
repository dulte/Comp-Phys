#pragma once
#include <potential.h>
#include <armadillo>


class Coloumb : public Potential
{
public:
    Coloumb(double omega_);
    double computePotential(double r);

private:
    double omega;
};

