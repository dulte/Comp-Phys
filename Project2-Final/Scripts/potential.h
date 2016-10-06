#pragma once

class Potential {
public:
    Potential();
    virtual double computePotential(double r) = 0; //Sets as an abstract method, so only the daughter classes, defined for specific potentials, can return a potential
};
