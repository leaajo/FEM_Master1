#include <stdlib.h>
#include <stdio.h>
const float precision = 1e-7;

void cal1Elem(int nRefDom,
                int nbRefD0,
                float* numRefD0,
                int nbRefD1,
                float *numRefD1,
                int nbRefF1,
                float *numRefF1,
                int typeEl,
                int nbneel,
                float **coorEl,
                int nbaret,
                float nRefArEl,
                float **MatElem, float *SMbrElem, float *NuDElem, float *uDElem){

    // Calcul de q
    int q;
    q = nppquad(typeEl);

    // calcul de xhat et omeka k
    float **xhat;  xhat = alloctab(q,2);
    float *omegak; 
    ppquad(typeEl ,q, omegak, xhat);

    // à l'intérieur du domaine :
    intElem(typeEl, q, nbneel, coorEl, omegak, xhat, MatElem, SMbrElem);

}