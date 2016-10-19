#include <iostream>
#include "vec3.h"
#include "system.h"
#include "particle.h"
#include "odesolver.h"

using namespace std;

int main(int argc, char *argv[])
{
    int steps_per_year = 10000;
    int years = 20;
    double dt = 1.0/(steps_per_year);

    System* sysESEuler = new System("positionsEarthSunEuler.xyz");
    ODEsolver* solver = new ODEsolver(dt,sysESEuler);


    vec3 posEarth = vec3(1.0,0.0,0.0);
    vec3 velEarth = vec3(0.0,2*M_PI,0.0);

    vec3 posSun = vec3(0.0,0.0,0.0);
    vec3 velSun = vec3(0.0,0.0,0.0);

    vec3 posJup = vec3(4.0,0,0);
    vec3 velJup = vec3(0.0,M_PI,0);
    sysESEuler->create_Particle(posSun,velSun,1);
    sysESEuler->create_Particle(posEarth,velEarth,3e-6);
    sysESEuler->create_Particle(posJup,velJup,9.54e-4);
    sysESEuler->create_Particle(posJup+vec3(0.05,0,0),velJup + vec3(0,1,0),9.54e-7);

    for (int i = 0; i< years*steps_per_year; i++){
        solver->eulerCromerOneStep();
        sysESEuler->dumpPositionsToFile();
    }





    return 0;
}
