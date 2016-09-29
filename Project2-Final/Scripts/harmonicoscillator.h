#pragma once
#include <potential.h>

class HarmonicOscillator : public Potential {
public:
    HarmonicOscillator(double omega_);
    double computePotential(double r);

    double omega;
};
