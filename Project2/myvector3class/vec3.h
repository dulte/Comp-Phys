#ifndef VEC3_H
#define VEC3_H

class vec3
{
public:
    double components[3];
    vec3();

    double lengthSquared();
    double length();
    double &operator()(int index) {return components[index];}
    double &operator[](int index) {return components[index];}
    vec3 &operator-=(vec3 rhs);
    vec3 &operator -=(double rhs);
};

inline vec3 operator-(vec3 lhs, vec3 rhs){
    lhs -= rhs;
    return lhs;
}

inline vec3 operator-(double lhs, vec3 rhs){
    rhs -= lhs;
    return rhs;
}

inline vec3 operator-(vec3 lhs, double rhs){
    lhs -= rhs;
    return lhs;
}

#endif // VEC3_H
