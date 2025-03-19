#include <stdio.h>
#include <stdlib.h>

// Paramètres :
// @t         : type de l'élément (quadrangle, triangle, rectangle)
// @q         : nombre de points de quadrature
// @nbneel    : nombre de noeuds par élément 
// @coorEl    : coordonnées des noeuds sur l'élément
// @omegak    : poids des points de quadrature 
// @xhat      : coordonnées des points de quadrature
// @vecElem   : calcul et stock de l'intégrale issue de W
// @matElem   : calcul et stock de l'intégrale issu de WW et ADWDW

void intAret( int t,
            int q,
            int nbneel,
            float **coorEl,
            float *omegak,
            float **xhat,
            float **matElem,
            float **vecElem){

    //------- Initialisation des variables -------------------------------------------------------------------
    int d =2; // dimension
    if (t==1) d=1;

    float wk[nbneel];                               // calcul de Fkx^
    float Fkxhat[d];                                // vecteur qui contiendra 

    float **dwk; dwk = alloctab(nbneel, d);         // 
    float **dwkEl; dwkEl = alloctab(nbneel,d);      // pour ADWDW

    float **dFkxhat; dFkxhat = alloctab(2,2);         // Jacobienne de Fk(x^)

    float **dFkxhatinv; dFkxhatinv = alloctab(2,2); // inverse de la matrice D(Fk(x^))
    float det;                                    // déterminant de celle-ci

    // Initialisation des coefficients qui nous serviront à remplir matElem et vecElem
    float cofvar_W, cofvar_WW;
    float eltdif; 

    //----------- Boucle pour le remplissage de matElem et vecElem----------------------------------------------
    for(int k =0; k<q; k++){

        // Calcul de eltdif
        calDerFbase(t, xhat[k], dwk);
        matJacob(q,d, dwk, coorEl[k], dFkxhat);
        det = invertM2x2(dFkxhat, dFkxhatinv);
        eltdif = omegak[k]*absf(det);

        // Calcul de cofvar_w
        calFbase(t, xhat[k], wk);
        transFk(q, wk, coorEl, Fkxhat);
        cofvar_W = FN(Fkxhat);

        // contribution avec la procédure W
        W(nbneel, wk, eltdif, cofvar_W, vecElem);

        // calcul de cofvar_ww
        cofvar_WW = BN(Fkxhat);

        // contribution avec la procédure WW
        WW(nbneel, wk, eltdif, cofvar_WW, matElem); 

        // libération de la mémoire 
        free(dwk); free(dwkEl); free(wk); free(Fkxhat); 
    }
}
