void intElem( int t,
            int q,
            int nbneel,
            float **coorEl,
            float *omegak,
            float **xhat,
            float **matElem,
            float *vecElem){

    for (int i=0; i<nbneel; i++){
        printf("vecElem[%d]= %f \n", i, vecElem[i]);
    }
    //------- Initialisation des variables -------------------------------------------------------------------
    int d =2; // dimension
    if (t==3) d=1;

    float *wk; wk = malloc(nbneel * sizeof(float));                        // calcul de Fkx^
    float *Fkxhat = malloc(d * sizeof(float));                             // vecteur qui contiendra les résultats de transFk

    float **valDerFbx = alloctab(nbneel,2); 
    
    float **dwkEl = alloctab(nbneel,2);                // pour ADWDW

    float **dFkxhat = alloctab(2,2);                   // Jacobienne de Fk(x^)

    float **dFkxhatinv = alloctab(2,2);                // inverse de la matrice D(Fk(x^))
    float det;                                         // déterminant de celle-ci

    // Initialisation des coefficients qui nous serviront à remplir matElem et vecElem
    float cofvar_W, cofvar_WW;
    float **cofvar_ADWDW = alloctab(2,2);
    float eltdif; 

    //----------- Boucle pour le remplissage de matElem et vecElem----------------------------------------------
    for(int k =0; k<q; k++){

        // Calcul de eltdif 
        calDerFbase(t, xhat[k], valDerFbx);
        matJacob(valDerFbx, d, q, coorEl, dFkxhat);
        det = invertM2x2(dFkxhat, dFkxhatinv);
        eltdif = omegak[k]*fabsf(det);


        // Calcul de cofvar_w
        calFbase(t, xhat[k], wk);
        transFk(wk, coorEl, q, Fkxhat); 
        cofvar_W = FOMEGA(Fkxhat);
        
        // contribution avec la procédure W
        W(nbneel, wk, eltdif, cofvar_W, vecElem);

        // calcul de cofvar_ww
        cofvar_WW = a00(Fkxhat);

        // contribution avec la procédure WW
        WW(nbneel, wk, eltdif, cofvar_WW, matElem); 

        // calcul de cofvar_ADWDW
        cofvar_ADWDW[0][0]= a11(Fkxhat);
        cofvar_ADWDW[0][1]= a12(Fkxhat);
        cofvar_ADWDW[1][0]= a12(Fkxhat); 
        cofvar_ADWDW[1][1]= a22(Fkxhat);

        for (int i=0; i<nbneel; i++){
            dwkEl[i][0] = valDerFbx[i][0]*dFkxhatinv[0][0] + valDerFbx[i][1]*dFkxhatinv[1][0];
            dwkEl[i][1] = valDerFbx[i][0]*dFkxhatinv[0][1] + valDerFbx[i][1]*dFkxhatinv[1][1];
        }

        // contribution avec la procédure ADWDW 
        ADWDW(dwkEl, cofvar_ADWDW, eltdif, nbneel, matElem);

    }
    // libération de la mémoire 
        matrix_free(valDerFbx);
        matrix_free(dwkEl); 
        free(wk);
        free(Fkxhat);
        matrix_free(dFkxhat);
        matrix_free(dFkxhatinv);
        matrix_free(cofvar_ADWDW);
}