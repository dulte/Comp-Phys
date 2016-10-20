#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H
#include "vec3.h"

class Celestial_Body
{
public:
    Celestial_Body();
    double mass_of_planet;
    vec3 position_of_planet;
    vec3 velocity_of_planet;

    double compute_energy(vec3 position_of_planet, vec3 velocity_of_planet, double mass_of_planet);
    Celestial_Body(vec3 position_of_planet, vec3 velocity_of_planet, double mass_of_planet);

};

#endif // CELESTIAL_BODY_H
