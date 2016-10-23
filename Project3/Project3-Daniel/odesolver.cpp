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

void ODEsolver::eulerOneStep()
{
    system->compute_acceleration();
    for (int i = 0; i < system->list_of_particles.size(); i++){
        for(int a=0; a<3; a++) {
            system->list_of_particles[i].m_position[a] += system->list_of_particles[i].m_velocity[a]*dt;
            system->list_of_particles[i].m_velocity[a] += system->list_of_particles[i].acceleration[a]*dt;
        }
    }
}

void ODEsolver::velocityVerletOneStep(){

    double old_a[system->number_of_bodies][3];
    for (int i = 0; i < system->list_of_particles.size(); i++){


        old_a[i][0] = system->list_of_particles[i].acceleration[0];
        old_a[i][1] = system->list_of_particles[i].acceleration[1];
        old_a[i][2] = system->list_of_particles[i].acceleration[2];

    }
    for (int i = 0; i < system->list_of_particles.size(); i++){
        for(int a=0; a<3; a++) {
            system->list_of_particles[i].m_position[a] += system->list_of_particles[i].m_velocity[a]*dt + dt*dt/2.0 *old_a[i][a] ;

        }
    }
    system->compute_acceleration();
    //system->compute_acceleration_relativistic(0,1);

    for (int i = 0; i < system->list_of_particles.size(); i++){
        for(int a=0; a<3; a++) {

            system->list_of_particles[i].m_velocity[a] += (system->list_of_particles[i].acceleration[a] + old_a[i][a])*dt/2.0;


        }
    }

}
