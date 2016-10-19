#include "odesolver.h"
#include "system.h"
#include <iostream>
#include <vector>

using namespace std;

ODEsolver::ODEsolver(double dt_, System *sys)
{
    dt = dt_;
    system = sys;
}

void ODEsolver::eulerCromerOneStep()
{
    system->compute_acceleration();
    //cout << system->number_of_bodies << endl;
    //system->list_of_particles[1].m_position.print();
    for (int i = 0; i < system->list_of_particles.size(); i++){
        //Particle particle = system->list_of_particles[i];
        //system->list_of_particles[i].m_position.print();
        system->list_of_particles[i].m_velocity += system->list_of_particles[i].acceleration*dt;
        system->list_of_particles[i].m_position += system->list_of_particles[i].m_velocity*dt;
    }
}

void ODEsolver::velocityVerletOneStep(vector<Particle*> particles){
    for (Particle* particle : particles){
        particle->m_velocity += 0.5*particle->acceleration*dt;
        particle->m_position += particle->m_velocity*dt;
    }

    system->compute_acceleration();

    for (Particle* particle : particles){
        particle->m_velocity += 0.5*particle->acceleration*dt;
    }

}
