#ifndef SYSTEM_H
#define SYSTEM_H
#include<armadillo>
#include "particle.h"
#include <vector>

using namespace::std;




class System
{
public:
    System();
    void create_Particle(vec3 initial_position, vec3 initial_velocity, double mass);
    int number_of_bodies;
    double Kinetic_energy(double t);
    double Potential_energy(double t);
    vec3 Angular_momentum(double t);
    void compute_acceleration();
    double G = 4*M_PI;
    vector<Particle> list_of_particles;

   };

#endif // SYSTEM_H
