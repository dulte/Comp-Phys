#include <iostream>
#include "vec3.h"
#include "system.h"
#include "particle.h"
#include "odesolver.h"
#include <string>
#include <sstream>
#include <ctime>
#include <fstream>
#include <iomanip>


using namespace std;

void makePlanets(string filename,System * sys);
void dumpParametersToFile(string filename, System* sys, double steps,double years);

int main(int argc, char *argv[])
{
    clock_t begin = clock();
    int steps_per_year = 1e7;//10000000;
    int years = 10;
    double dt = 1.0/(steps_per_year);

    System* sysESEuler = new System("positionsMercurySun2.xyz");
    ODEsolver* solver = new ODEsolver(dt,sysESEuler);


    string folder = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project3/";
    makePlanets(folder + "planetData_merc.txt", sysESEuler);
    //sysESEuler->setInitialEnergyAndMomentum();
    int printEvery = 100;

    for (int i = 0; i< years*steps_per_year; i++){
        if(i % printEvery == 0) {
            sysESEuler->dumpPositionsToFile();
            //sysESEuler->dumpErrorToFile();
            //cout << setprecision(15) << sysESEuler->total_potential_energy + sysESEuler->Kinetic_energy() << endl;
        }
        //solver->eulerOneStep();
        solver->velocityVerletOneStep();
        sysESEuler->dumpThetaToFile();
    }

    dumpParametersToFile("parameters.txt",sysESEuler,steps_per_year*years,years);


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

        if (name != "Sun"){
            sys->create_Particle(vec3(x,y,z),vec3(vx,vy,vz),mass);

        }

    }

    vec3 totalMomentum = vec3();
    vec3 sumPositions = vec3();

    for (Particle & part: sys->list_of_particles){
        totalMomentum += part.m_velocity*part.m_mass;
        sumPositions += part.m_position*part.m_mass;
    }

    //sys->create_Particle(-1*sumPositions,-1*totalMomentum,1);
    sys->create_Particle(vec3(),vec3(),1);

}

void dumpParametersToFile(string filename, System* sys, double steps,double years){
    ofstream outFile(filename);
    outFile << "NumberOfPlanets" << " " << sys->number_of_bodies << "\n";
    outFile << "Years" << " " << years << "\n";
    outFile << "Timesteps" << " " << steps << "\n";

    outFile.close();

}
