#include "vec3.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>

using namespace std;

vec3::vec3()
{
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
}

vec3::vec3(double x, double y, double z)
{
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

vec3 &vec3::operator-=(vec3 rhs)
{
    components[0] -= rhs[0];
    components[1] -= rhs[1];
    components[2] -= rhs[2];
    return *this;
}

vec3 &vec3::operator+=(vec3 rhs)
{
    components[0] += rhs[0];
    components[1] += rhs[1];
    components[2] += rhs[2];
    return *this;
}

vec3 &vec3::operator *=(double c){
    components[0] *= c;
    components[1] *= c;
    components[2] *= c;
    return *this;
}


vec3 &vec3::operator /=(double c){
    components[0] /= c;
    components[1] /= c;
    components[2] /= c;
    return *this;
}

vec3 &vec3::operator /=(vec3 rhs){
    components[0] /= rhs[0];
    components[1] /= rhs[1];
    components[2] /= rhs[2];
    return *this;
}

vec3 &vec3::operator=(vec3 rhs)
{
    components[0] = rhs[0];
    components[1] = rhs[1];
    components[2] = rhs[2];
    return *this;
}

vec3 &vec3::operator=(double s)
{
    components[0] = s;
    components[1] = s;
    components[2] = s;
    return *this;
}



vec3 vec3::randint(int min, int max)
{
    components[0] = rand()%(max-min + 1) + min;
    components[1] = rand()%(max-min + 1) + min;
    components[2] = rand()%(max-min + 1) + min;
    return *this;
}

void vec3::normalize(){
    *this = *this/length();
}

void vec3::print(){
    cout << "[" << components[0] << "," << components[1] << "," << components[2] << "]" << endl;
}

vec3 vec3::cross(vec3 otherVec){
    return vec3(y()*otherVec.z()-z()*otherVec.y(), (z()*otherVec.x()-x()*otherVec.z()), x()*otherVec.y()-y()*otherVec.x());
}

double vec3::lengthSquared()
{
    return components[0]*components[0] + components[1]*components[1] + components[2]*components[2];
}

double vec3::length()
{
   return sqrt(lengthSquared());
}
