#pragma once
#include <potential.h>

class HarmonicOscillator : public Potential {
public:
    HarmonicOscillator();
    double computePotential(double r);
};
