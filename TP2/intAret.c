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
            int nbnAr,
            float **coorAr,
            float *omegak,
            float **xhat,
            float **matAret,
            float **vecAret){

    //------- Initialisation des variables -------------------------------------------------------------------
    int d =1;

    float wk[nbnAr];                               // calcul de Fkx^
    float Fkxhat[d];                                // vecteur qui contiendra 

    float **dwk; dwk = alloctab(q, 2);         // 
    float **dwkEl; dwkEl = alloctab(q,2);      // pour ADWDW

    float **dFkxhat; dFkxhat = alloctab(2,2);         // Jacobienne de Fk(x^)

    // Initialisation des coefficients qui nous serviront à remplir matElem et vecElem
    float cofvar_W, cofvar_WW;
    float eltdif; 

    //----------- Boucle pour le remplissage de matElem et vecElem----------------------------------------------
    for(int k =0; k<q; k++){

        // Calcul de eltdif
        calDerFbase(t, xhat[k], dwk);
        matJacob(q,d, dwk, coorAr[k], dFkxhat);
        eltdif = omegak[k]*sqrtf(dFkxhat[0][0] * dFkxhat[0][0] + dFkxhat[1][0] * dFkxhat[1][0]);

        // Calcul de cofvar_w
        calFbase(t, xhat[k], wk);
        transFk(q, wk, coorAr, Fkxhat);
        cofvar_W = FN(Fkxhat);

        // contribution avec la procédure W
        W(nbnAr, wk, eltdif, cofvar_W, vecAret);

        // calcul de cofvar_ww
        cofvar_WW = BN(Fkxhat);

        // contribution avec la procédure WW
        WW(nbnAr, wk, eltdif, cofvar_WW, matAret); 
    }
}
