#include "vec3.h"
#include <cmath>

vec3::vec3()
{
 components[0]=0;
 components[1]=0;
 components[2]=0;
}

double vec3::lengthSquared()
{
        return components[0]*components[0]+components[1]*components[1]+components[2]*components[2];
}

double vec3::length(){

    return sqrt(lengthSquared());
}
vec3 &vec3::operator-=(vec3 rhs)
{
    components[0]-=rhs[0];
    components[1]-=rhs[1];
    components[2]-=rhs[2];
    return *this;
}

vec3 &vec3::operator-=(double rhs)
{
    components[0]-=rhs;
    components[1]-=rhs;
    components[2]-=rhs;
    return *this;
}


