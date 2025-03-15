#include <stdlib.h>
#include <stdio.h>
const float precision = 1e-7;

float a11 (float *x){
    return 1;
}

float a22 (float *x){
    return 1;
}

float a00 (float *x){
    return 1;
}

float a12 (float *x){
    return 0;
}

float bN (float *x){
    return 0;
}

float FOMEGA (float *x){
    return 1;
}

float FN (float *x){
    return 1;
}

float UD (float *x){
    return 100* x[0] + x[1];
}
