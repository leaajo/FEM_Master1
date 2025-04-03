#include <stdio.h>
#include <stdlib.h>

void main(){
    int Nblign;
    int NbCoef;
    int NbCoef;
    float *SecMembre   = malloc(Nblign*sizeof(float));
    float *ValDlDir    = malloc(Nblign*sizeof(float));
    float *AdPrCoefLi  = malloc(Nblign*sizeof(float));

    // Ils dépendent de AdPrCoefLi 
                        NbCoef = AdPrCoefLi[Nblign]-1;
    float *Matrice     = malloc((Nblign+NbCoef)*sizeof(float));
    float *NumCol      = malloc(NbCoef*sizeof(float));
    float *AdSuccLi    = malloc(NbCoef*sizeof(float));

    //char* ficmai = "car3x3t_3";
    char ficmai[] = "car3x3t_3.txt"; 
    int ptypeEl, pnbtng, pnbtel, pnbeel, pnbaret;
    int **pngnel;
    int **pnRefAr;
    float **coorAll;

    lecfima(&(*ficmai), &ptypeEl, &pnbtng, &coorAll, &pnbtel, &pngnel, &pnbeel, &pnbaret, &pnRefAr);


    int *numRefD0; // = (int *)malloc(sizeof(int) * 1); numRefD0[0] = 1;
    int *numRefD1; // = (int *)malloc(sizeof(int) * 1); numRefD1[0] = 4;
    int *numRefF1; // = (int *)malloc(sizeof(int) * 2); numRefF1[0] = 2; numRefF1[1] = 3;
    int nRefDom; //= 0;
    int nbRefD0; // = 1;
    int nbRefD1; // = 1;
    int nbRefF1; // = 2;

    char* fic_ref = "NUMREF.Test";
    lecture_conditions_bords(fic_ref, &numRefD0, &numRefD1, &numRefF1, &nbRefD0, &nbRefD1, &nbRefF1, &nRefDom);


    float **coorEl = alloctab(pnbeel, 2);
    float **MatElem; MatElem = alloctab(pnbeel, pnbeel);
    float SMbrElem[pnbeel];
    int NuDElem[pnbeel];
    float uDElem[pnbeel];
    int pnRefArEl[pnbaret];

    int* NB = malloc(pnbeel * sizeof(int));

    for (int K = 0; K<pnbtel; K++){

        for(int i = 0; i<pnbeel; i++){
            NB[i]= pngnel[K][i];
        }
        selectPts(pnbeel, NB, coorAll, coorEl);
        for (int i =0; i <pnbaret; i++){
            pnRefArEl[i]=pnRefAr[K][i];
        }

        cal1Elem(nRefDom, nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1, 
                ptypeEl, pnbeel, coorEl, pnbaret, pnRefArEl,
                MatElem, SMbrElem, NuDElem, uDElem);
        
        impCalEl(K+1, ptypeEl, pnbeel, MatElem, SMbrElem, NuDElem, uDElem);
    }
    return 0;

}