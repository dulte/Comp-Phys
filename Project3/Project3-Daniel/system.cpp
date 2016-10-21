#include "system.h"
#include "particle.h"
#include "vec3.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

System::System(string filename)
{
    number_of_bodies = 0;
    outFile.open(filename);
    outFileEnergyMomentum.open("errorEnergyMomentum.txt");
}

System::~System(){
    outFile.close();
    outFileEnergyMomentum.close();
}

void System::create_Particle(vec3 initial_position, vec3 initial_velocity, double mass){
    list_of_particles.push_back(Particle(initial_position, initial_velocity, mass));
    number_of_bodies=list_of_particles.size();

}

void System::compute_acceleration(){
    //double total_potential_energy=0;
    total_potential_energy = 0.0;
    for(Particle &body : list_of_particles) {
        body.acceleration=0.0;
    }

    for(int i=0; i < number_of_bodies; i++){
        Particle &body_i = list_of_particles[i];

        for(int j=i+1; j<number_of_bodies; j++){
            Particle &body_j = list_of_particles[j];

            double dx = body_j.m_position[0] - body_i.m_position[0];
            double dy = body_j.m_position[1] - body_i.m_position[1];
            double dz = body_j.m_position[2] - body_i.m_position[2];

            double r2 = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r2);
            double prefac = G*body_i.m_mass*body_j.m_mass/(r*r2);

            double Fx = dx*prefac;
            double Fy = dy*prefac;
            double Fz = dz*prefac;

            body_i.acceleration[0] += Fx * body_i.m_massInverse;
            body_i.acceleration[1] += Fy * body_i.m_massInverse;
            body_i.acceleration[2] += Fz * body_i.m_massInverse;

            body_j.acceleration[0] -= Fx * body_j.m_massInverse;
            body_j.acceleration[1] -= Fy * body_j.m_massInverse;
            body_j.acceleration[2] -= Fz * body_j.m_massInverse;
            total_potential_energy -= prefac*r2;

            //cout << prefac << endl;

        }
        //current_body.acceleration.print();
    }
    // total_potential_energy=total_potential_energy/2.0;
    //return total_potential_energy;
}

void System::compute_acceleration_relativistic(int index_body, int index_star){
    double c = 63198;
    vec3 angular_moment;
    angular_moment= System::Angular_momentum("one", index_body);
    Particle current_body=list_of_particles[index_body];
    current_body.acceleration=0.0;
    vec3 r = (list_of_particles[index_star].m_position-current_body.m_position);
    double relativistic_factor=1+(3.0*angular_moment.lengthSquared())/(r.lengthSquared()*c*c);
    double prefac=G*list_of_particles[index_star].m_mass/(r.lengthSquared()*r.length());
    current_body.acceleration=r*relativistic_factor*prefac;
}


void System::setInitialEnergyAndMomentum(){

    compute_acceleration();
    initialAngularMomentum = (Angular_momentum("many",0)).length();
    initialEnergy = total_potential_energy;
}



void System::dumpPositionsToFile()
{
    outFile << number_of_bodies << endl;
    outFile << "New Time step" << endl;
    for (Particle particle:list_of_particles){
        outFile << particle.m_position[0] << " " << particle.m_position[1] << " " << particle.m_position[2] << "\n";
    }

}

void System::dumpErrorToFile()
{
    compute_acceleration();
    double errorMomentum = (abs(((Angular_momentum("many",0)).length() - initialAngularMomentum))/(initialAngularMomentum));
    double errorEnergy = (abs((total_potential_energy - initialEnergy))/initialEnergy);

    outFileEnergyMomentum << errorMomentum << " " << errorEnergy << "\n";
}

vec3 System::Angular_momentum(char* one_or_many, int body){
    vec3 total_angular_momentum;
    if (strcmp(one_or_many, "many") == 0){
        for(int j=0; j<number_of_bodies; j++){
            vec3 linear_momentum=(list_of_particles[j].m_velocity)*(list_of_particles[j].m_mass);
            vec3 angular_momentum=(list_of_particles[j].m_position).cross(linear_momentum);
            total_angular_momentum += angular_momentum;
        }
    }

    else if (strcmp(one_or_many, "one") == 0){
        vec3 linear_momentum=(list_of_particles[body].m_velocity)*(list_of_particles[body].m_mass);
        vec3 angular_momentum=(list_of_particles[body].m_position).cross(linear_momentum);
        total_angular_momentum = angular_momentum;
    }

    else{
        cout << "Need to specify if computing angular momentum for one body or many" << endl;
        exit(0);
    }
    return total_angular_momentum;
}

double System::Kinetic_energy(){
    double total_kinetic_energy=0;
    for(int j=0; j<number_of_bodies; j++){
        double speed_squared=(list_of_particles[j].m_velocity).lengthSquared();
        total_kinetic_energy += 0.5*(list_of_particles[j].m_mass)*speed_squared;
    }
    return total_kinetic_energy;
}