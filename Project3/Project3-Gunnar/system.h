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
    double Kinetic_energy();
    vec3 Angular_momentum(char* one_or_many, int body);
    double compute_acceleration();
    void compute_acceleration_relativistic(int index_body, int index_star);
    double G = 4*M_PI;
    vector<Particle> list_of_particles;

   };

#endif // SYSTEM_H
