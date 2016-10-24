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
    int steps_per_year = 1e6;
    int years = 10;
    double dt = 1.0/(steps_per_year);

    System* sys = new System("positions.xyz");
    ODEsolver* solver = new ODEsolver(dt,sys);


    string folder = "/home/daniel/Dokumenter/Skole/Comp-Phys/Project3/";
    makePlanets(folder + "planetData.txt", sys);
    sys->setInitialEnergyAndMomentum();
    int printEvery = 100;

    for (int i = 0; i< years*steps_per_year; i++){
        //solver->eulerOneStep();
        solver->velocityVerletOneStep();
        sys->dumpThetaToFile();
        if(i % printEvery == 0) {
            sys->dumpPositionsToFile();
            sys->dumpErrorToFile();
        }
    }

    dumpParametersToFile("parameters.txt",sys,steps_per_year*years,years);


    clock_t end = clock();
    double elapsed_secs = double(end - begin) / double(CLOCKS_PER_SEC);
    cout << "Finished simulation after " << elapsed_secs << " seconds." << endl;

    return 0;
}


//Gets the initial data of the planets from the file made by planetInfo.py
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

        if (name != "Sun" && name != ""){
            sys->create_Particle(vec3(x,y,z),vec3(vx,vy,vz),mass);

        }

    }

    //The code below gives the sun an initial position and velocity such that the center of mass is always at (0,0,0)
    vec3 totalMomentum = vec3();
    vec3 sumPositions = vec3();

    for (Particle & part: sys->list_of_particles){
        totalMomentum += part.m_velocity*part.m_mass;
        sumPositions += part.m_position*part.m_mass;
    }

    sys->create_Particle(-1*sumPositions,-1*totalMomentum,1);
    //sys->create_Particle(vec3(),vec3(),1);

}

void dumpParametersToFile(string filename, System* sys, double steps,double years){
    ofstream outFile(filename);
    outFile << "NumberOfPlanets" << " " << sys->number_of_bodies << "\n";
    outFile << "Years" << " " << years << "\n";
    outFile << "Timesteps" << " " << steps << "\n";


    outFile.close();

}
