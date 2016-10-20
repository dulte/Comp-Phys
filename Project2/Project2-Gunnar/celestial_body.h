#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H

class Celestial_Body
{
public:
    Celestial_Body();
    double mass_of_planet;
    vec initial_pos;
    vec initial_vel;
    double Celestial_Body::acceleration(vec position);

};

#endif // CELESTIAL_BODY_H
