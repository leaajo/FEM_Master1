#include <stdlib.h>
#include <stdio.h>
const float precision = 1e-7;

void WW (int nbneel, float *fctbas, float eltdif, float cofvar, float **matelem){
    int i, j;
    float coeff;

    for (i=0; i<nbneel; i++){
        coeff = eltdif* cofvar* fctbas[i];
        for(j=0; j<nbneel; j++){
            matelem[i][j]= matelem[i][j]+ coeff*fctbas[j];
        }
    }
}

void W (int nbneel, float *fctbas, float eltdif, float cofvar, float *matelem){
    int i;
    float coeff;

    for (i=0; i<nbneel; i++){
        coeff = eltdif* cofvar* fctbas[i];
        matelem[i]= matelem[i]+coeff;      
    }
}

void ADWDW(float **derfctbas, float **cofvar, float eltdif, int nbneel, float **matelem){
    int alpha, beta;
    float coeff; 

    for (alpha =0; alpha <2; alpha++){
        for (beta = 0; beta <2; beta++){
            for (int i =0; i<nbneel-1; i++){
                 coeff = eltdif*cofvar[alpha][beta]*derfctbas[i][alpha];
                 for (int j =0; j<nbneel-1; j++){
                    matelem[i][j]= matelem[i][j]+coeff*derfctbas[j][beta];
                 }
            }
        }
    }
}
