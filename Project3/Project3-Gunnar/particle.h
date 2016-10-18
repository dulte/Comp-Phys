#ifndef PARTICLE_H
#define PARTICLE_H
#include <armadillo>
#include "vec3.h"


class Particle
{
public:
    Particle(vec3 position, vec3 velocity, double mass);
    double m_mass;
    vec3 m_position;
    vec3 m_velocity;
    vec3 acceleration;
};

#endif // PARTICLE_H
