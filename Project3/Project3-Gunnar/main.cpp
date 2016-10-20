#include <iostream>
#include "vec3.h"
#include "particle.h"
#include "system.h"

using namespace std;

int main()
{
    System solar_system;
    vec3 initial_position;
    vec3 initial_velocity;
    initial_position= vec3(1,0,0);
    initial_velocity= vec3(0,1,1);
    double mass=100;
    solar_system.create_Particle(initial_position, initial_velocity, mass);
    double ek=solar_system.Kinetic_energy();
    cout << ek << endl;

    return 0;
}
