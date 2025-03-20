// test de intElem
#include <stdio.h>
#include <stdlib.h>
#include "decl_tp2.h"

int main(){

    int nRefDom = 0;
    int nbRefD0 = 0;
    int nbRefD1 = 0;
    int nbRefF1 = 0;

    float *numRefD0, *numRefD1, *numRefF1;
    numRefD0 = malloc(2*sizeof(float));
    numRefD1 = malloc(3*sizeof(float));
    numRefF1 = malloc(2*sizeof(float));


    //char* ficmai = "car3x3t_3";
    char ficmai[] = "car1x1t_1.txt"; 
    int ptypeEl, pnbtng, pnbtel, pnbeel, pnbaret;
    int **pngnel;
    int **pnRefArEl;
    float **coorAll;

    lecfima(&(*ficmai), &ptypeEl, &pnbtng, &coorAll, &pnbtel, &pngnel, &pnbeel, &pnbaret, &pnRefArEl);

    float **coorEl = alloctab(pnbeel, 2);
    float **MatElem; MatElem = alloctab(pnbeel, pnbeel);
    float SMbrElem[pnbeel];
    int NuDElem[pnbeel];
    float uDElem[pnbeel];

    int* NB = malloc(pnbeel * sizeof(int));

    for (int K = 0; K<pnbtel; K++){

        for( int i = 0; i<pnbeel; i++){
            NB[i]= pngnel[K][i];
        }

        selectPts(pnbeel, NB, coorAll, coorEl);

        cal1Elem(nRefDom, nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1, 
                ptypeEl, pnbeel, coorEl, pnbaret, &pnRefArEl,
                MatElem, SMbrElem, NuDElem, uDElem);
        
        impCalEl(K, ptypeEl, pnbeel, MatElem, SMbrElem, NuDElem, uDElem);
    }
    return 0;
}