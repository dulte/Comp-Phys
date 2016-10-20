#include "system.h"
#include "particle.h"
#include "vec3.h"

System::System()
{

}


void System::create_Particle(vec3 initial_position, vec3 initial_velocity, double mass){
    list_of_particles.push_back(Particle(initial_position, initial_velocity, mass));
}

void System::compute_acceleration(){
    int number_of_bodies=list_of_particles.size();
    for(int i=0; i < number_of_bodies; i++){
        Particle current_body=list_of_particles[i];
        current_body.acceleration=0.0;
        for(int j=0; j<number_of_bodies; j++){
            if(i != j){
                vec3 r = (list_of_particles[j].position-current_body.position);
                current_body.acceleration+= r*G*list_of_particles[j].mass/((r.lengthSquared()*r.length()));
            }
        }
    }
}
