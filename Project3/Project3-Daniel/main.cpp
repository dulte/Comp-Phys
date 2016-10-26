#include <iostream>
#include "vec3.h"
#include "system.h"
#include "particle.h"
#include "odesolver.h"
#include <string>
#include <sstream>
#include <ctime>

using namespace std;

void makePlanets(string filename,System * sys);

int main(int argc, char *argv[])
{
    clock_t begin = clock();
    int steps_per_year = 10000;
    int years = 30;
    double dt = 1.0/(steps_per_year);

    System* sysESEuler = new System("positionsEarthSunEuler.xyz");
    ODEsolver* solver = new ODEsolver(dt,sysESEuler);


    string folder = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project3/";
    makePlanets(folder + "planetData.txt", sysESEuler);
    sysESEuler->setInitialEnergyAndMomentum();
    int printEvery = 100;

    for (int i = 0; i< years*steps_per_year; i++){
        solver->eulerOneStep();
        if(i % printEvery == 0) {
            sysESEuler->dumpPositionsToFile();
            sysESEuler->dumpErrorToFile();
        }
    }


    clock_t end = clock();
    double elapsed_secs = double(end - begin) / double(CLOCKS_PER_SEC);
    cout << "Finished simulation after " << elapsed_secs << " seconds." << endl;

    return 0;
}


void makePlanets(string filename,System * sys){

    ifstream inFile(filename);

    inFile.clear();
    inFile.seekg(0,ios::beg);

    for (string l; getline(inFile,l);)
    {
        stringstream ss(l);

        string name;
        double x,y,z,vx,vy,vz,mass;

        ss >> name >> x >> y >> z >> vx >> vy >> vz >> mass;

        sys->create_Particle(vec3(x,y,z),vec3(vx,vy,vz),mass);

    }


}
