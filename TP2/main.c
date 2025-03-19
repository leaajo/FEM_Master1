// test de intElem
#include <stdio.h>
#include <stdlib.h>
#include "decl_tp2.h"
const float precision = 1e-7;

int main(){
    int nRefDom = 0;
    int nbRefD0 = 0;
    float* numRefD0;
    int nbRefD1 = 0;
    float *numRefD1;
    int nbRefF1 = 0;
    float *numRefF1; 

    char* ficmai = "car1x1t_1";
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

    int NB[pnbeel];
    for (int K = 0; K<pnbtel; K++){
        for( int i = 0; i<pnbeel; i++){
            NB[i]= pngnel[K][i];
        }
        selectPts(pnbeel, NB, coorAll, coorEl);
        cal1Elem(nRefDom, nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1, ptypeEl, pnbeel, coorEl, pnbaret, pnRefArEl,
                 MatElem, SMbrElem, NuDElem, uDElem);
        impCalEl(K, ptypeEl, pnbeel, MatElem, SMbrElem, NuDElem, uDElem);
    }
    return 0;
}