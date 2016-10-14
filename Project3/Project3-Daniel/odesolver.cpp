#include "odesolver.h"
#include <vector>

using namespace std;

ODEsolver::ODEsolver(double dt_, System *system)
{
    dt = dt_;
}

void ODEsolver::eulerCromerOneStep(vector<Particle *> particles)
{
    system->compute_acceleration();
    for (Particle* particle : particles){
        particle->velocity += particle->acceleration*dt;
        particle->postions += particle->velocity*dt;
    }
}

void velocityVerletOneStep(vector<Particle*> particles){
    for (Particle* particle : particles){
        particle->velocity += 0.5*particle->acceleration*dt;
        particle->postions += particle->velocity*dt;
    }

    system->compute_acceleration();

    for (Particle* particle : particles){
        particle->velocity += 0.5*particle->acceleration*dt;
    }

}
