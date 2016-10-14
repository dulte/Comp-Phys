#pragma once
#include <vector>

using namespace std;



class ODEsolver
{
public:
    ODEsolver(double dt_, System* system);
    void velocityVerlet(vector<Particle*> particles);
    void eulerCromerOneStep(vector<Particle*> particles);
    double dt;

};
