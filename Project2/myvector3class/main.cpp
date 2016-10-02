#include <stdio.h>
#include <iostream>
#include "vec3.h"

using namespace std;

int main(int argc, char *argv[])
{
    vec3 a;
    vec3 b;
    a(0)+=2;
    b(1)=5;
    a-b;
    b -= 5;


    return 0;
}

