#include "celestial_body.h"
#include "vec3.h"


Celestial_Body::Celestial_Body()
{
}

double Celestial_Body::compute_energy(vec3 position_of_planet, vec3 velocity_of_planet, double mass_of_planet){

    double kinetic_energy;
    double potential_energy;
    double total_energy;

    kinetic_energy=0.5*mass_of_planet*velocity_of_planet.lengthSquared();
}
