#include "particle.h"
#include "system.h"

Particle::Particle(vec3 position, vec3 velocity, double mass)
{
    m_mass=mass;
    if (m_mass != 0){
        m_massInverse = 1.0 / mass;
    }else{
        m_massInverse = 0;
    }
    m_velocity=velocity;
    m_position=position;


}
