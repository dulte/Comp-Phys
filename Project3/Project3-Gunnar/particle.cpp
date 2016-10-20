#include "particle.h"
#include "system.h"

Particle::Particle(vec3 position, vec3 velocity, double mass)
{
    m_mass=mass;
    m_velocity=velocity;
    m_position=position;


}
