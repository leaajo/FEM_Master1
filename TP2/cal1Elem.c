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
                int *nRefArEl, // à construire
                float **MatElem, float *SMbrElem, int *NuDElem, float *uDElem){

    for (int i = 0; i<nbneel; i++){
        for (int j = 0; j<nbneel; j++){
            MatElem[i][j]=0;
        } 
        SMbrElem[i]=0;
        NuDElem[i]=1;
        uDElem[i]=0;   
    }
    // Calcul de q
    int q;
    q = nppquad(typeEl);

    // calcul de xhat et omega k
    float **xhat;  xhat = alloctab(q,2);
    float omegak[q];
    ppquad(typeEl ,q, omegak, xhat);

    // à l'intérieur du domaine :
    intElem(typeEl, q, nbneel, coorEl, omegak, xhat, MatElem, SMbrElem);

    float **coorAr= alloctab(2,2);

    for(int i=0; i<nbaret; i++){
        if (nRefArEl[i]!=nRefDom){

        // On récupère les coordonnées des noeuds liés à l'arête i :
         int *numNeAret = (int*)malloc(2*sizeof(int));
         numNaret(typeEl, i, numNeAret);
         selectPts(2,numNeAret, coorEl, coorAr);

         for (int l=0; l<2; l++){

            // Condition de Dirichlet homogène
            for (int j= 0; j<nbRefD0; j++){
                if (nRefArEl[i]==numRefD0[j]){
                    NuDElem[ numNeAret[l]-1 ]=0;
                }
            }

            // Condition de Dirichlet non-homogène
            for (int j= 0; j<nbRefD1; j++){
                if (nRefArEl[i]==numRefD1[j]){
                    NuDElem[ numNeAret[l]-1 ]=-1;
                    uDElem[ numNeAret[l]-1 ] = UD(coorAr);
                }
            }
         }

        // Condition de Neumann ou Fourier 

        }
    }
    

}