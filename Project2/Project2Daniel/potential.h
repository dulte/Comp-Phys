#pragma once

class Potential {
public:
    Potential();
    virtual double computePotential(double r) = 0;
};
