#include <iostream>
#include "vec3.h"

using namespace std;

int main()
{
    vec3 A;
    vec3 B;
    vec3 C;
    A=vec3 (1, 0, 0);
    B=vec3 (0, 1, 0);
    C=A.cross(B);
    C.print();
    return 0;
}

