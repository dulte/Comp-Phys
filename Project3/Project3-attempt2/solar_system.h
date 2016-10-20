#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H
#include "celestial_body.h"


class Solar_System
{
public:
    Solar_System();
    vec3 compute_acceleration(vec3 position, vec3 velocity);
    double compute_energy(vec3 position, vec3 velocity);

};

#endif // SOLAR_SYSTEM_H
