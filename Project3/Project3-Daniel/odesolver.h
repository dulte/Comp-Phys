#pragma once
#include "particle.h"
#include "system.h"
#include <vector>

using namespace std;



class ODEsolver
{
public:
    ODEsolver(double dt_, System* sys);
    void velocityVerletOneStep(vector<Particle*> particles);
    void eulerOneStep();
    double dt;
    System *system;

};
