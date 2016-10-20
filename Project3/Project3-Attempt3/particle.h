#ifndef PARTICLES_H
#define PARTICLES_H
#include <armadillo>
#include "vec3.h"


class Particle
{
public:
    Particle(vec3 position, vec3 velocity, double mass);
    double mass;
    vec3 position;
    vec3 velocity;
    vec3 acceleration;
};

#endif // PARTICLES_H
