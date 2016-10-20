#include "system.h"
#include "particle.h"
#include "vec3.h"

System::System()
{

}


void System::create_Particle(vec3 initial_position, vec3 initial_velocity, double mass){
    list_of_particles.push_back(Particle(initial_position, initial_velocity, mass));
    int number_of_bodies=list_of_particles.size();

}

double System::compute_acceleration(){
    total_potential_energy=0;
    for(int i=0; i < number_of_bodies; i++){
        Particle current_body=list_of_particles[i];
        current_body.acceleration=0.0;
        for(int j=0; j<number_of_bodies; j++){
            if(i != j){
                vec3 r = (list_of_particles[j].position-current_body.position);
                double prefac=G*list_of_particles[j].mass/(r.length());
                current_body.acceleration+= prefac*(r/(r.lengthSquared()));
                total_potential_energy -= prefac*list_of_particles[i].mass;
            }
        }
    }
    total_potential_energy=total_potential_energy/2.0;
    return total_potential_energy;
}



vec3 System::Angular_momentum(){
    vec3 total_angular_momentum;
    for(int j=0; j<number_of_bodies; j++){
        vec3 linear_momentum=list_of_particles[j].mass*list_of_particles[j].velocity;
        vec3 angular_momentum=(list_of_particles[j].position).cross(linear_momentum);
        total_angular_momentum += angular_momentum;
    }
    return total_angular_momentum;

}

double System::Kinetic_energy(){
    double total_kinetic_energy=0;
    for(int j=0; j<number_of_bodies; j++){
        double speed_squared=(list_of_particles[j].velocity).lengthSquared();
        total_kinetic_energy += 0.5*(list_of_particles[j].mass)*speed_squared;
    }
    return total_kinetic_energy;
}

