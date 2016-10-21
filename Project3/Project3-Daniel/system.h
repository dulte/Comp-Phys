#pragma once
#include<armadillo>
#include "particle.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace::std;




class System
{
public:
    System(string filename);
    ~System();

    void create_Particle(vec3 initial_position, vec3 initial_velocity, double mass);
    int number_of_bodies;
    double Kinetic_energy();
    vec3 Angular_momentum(char* one_or_many, int body);
    void compute_acceleration();
    void compute_acceleration_relativistic(int index_body, int index_star);

    void setInitialEnergyAndMomentum();
    void dumpThetaToFile();

    double G = 4*M_PI*M_PI;
    vector<Particle> list_of_particles;
    ofstream outFile;
    ofstream outFileEnergyMomentum;
    ofstream outFileTheta;

    double total_potential_energy;
    double total_energy;
    double initialEnergy;
    double initialAngularMomentum;

    double theta = 0;
    double thetaPrev = 0;
    double prevTheta = 0;
    double rPrev = 0;
    double rPrevPrev = 0;



    void dumpPositionsToFile();
    void dumpErrorToFile();

};

