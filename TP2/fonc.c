#include <stdlib.h>
#include <stdio.h>
#include <math.h>
const float precision = 1e-7;

// FONCTIONS DU TP1
int **alloctabi(int dim1, int dim2) {
  int **ptr;

  ptr = malloc(dim1 * sizeof(int *));
  if (ptr != NULL) {
    int i;
    int *tmp = malloc(dim1 * dim2 * sizeof(int));  // Allocation de la mémoire pour toutes les lignes
    if (tmp != NULL) {
      for (i = 0; i < dim1; i++) {
        ptr[i] = tmp;    // Assigner le début de chaque ligne
        tmp += dim2;     // Avancer de dim2 éléments pour la ligne suivante
      }
    }
    else {
      free(ptr);
      ptr = NULL;
    }
  }
  return ptr;
}

float **alloctab(int dim1, int dim2) {
  float **ptr;

  ptr = malloc(dim1 * sizeof(float *));
  if (ptr != NULL) {
    int i;
    float *tmp = malloc(dim1 * dim2 * sizeof(float));  // Allocation de la mémoire pour toutes les lignes
    if (tmp != NULL) {
      for (i = 0; i < dim1; i++) {
        ptr[i] = tmp;     // Assigner le début de chaque ligne
        tmp += dim2;      // Avancer de dim2 éléments pour la ligne suivante
      }
    }
    else {
      free(ptr);
      ptr = NULL;
    }
  }
  return ptr;
}


void matrix_free_int(int** ppMatrix) {
    if (ppMatrix) {
        free(ppMatrix[0]); 
        free(ppMatrix);    
    }
}


void matrix_free(float** ppMatrix) {
    if (ppMatrix) {
        free(ppMatrix[0]); 
        free(ppMatrix);    
    }
}


void etiqAr(int t, int n1_1, int n2_1, int nrefdom, int *nrefcot, int nbtel, int nbaret, int **nRefAr){

    // On initialise le nRefAr à nrefdom :
    for (int i =0; i< nbtel; i++){
        for (int j = 0; j< nbaret; j++) {
            nRefAr[i][j]= nrefdom;
        }
    }

  if (t==1){ // Quadrangle
    for (int i = 0; i<n1_1; i++){
      nRefAr[i][3] = nrefcot[0];
    }
    for (int i = n1_1*n2_1-1; i>(n1_1)*(n2_1)-n1_1-1; i--){
      nRefAr[i][1] = nrefcot[2];
    }
    for (int i = n1_1 -1; i<(n1_1)*(n2_1); i+=n1_1){
      nRefAr[i][0] = nrefcot[1];
    }
    for (int i = 0; i<(n1_1)*(n2_1); i+=n1_1){
      nRefAr[i][2] = nrefcot[3];
    }
  }
  if (t==2){ // Triangle
    for (int i = 0; i<2*(n1_1); i+=2){
      nRefAr[i][2] = nrefcot[0];
    }
    for (int i = 2*n1_1*n2_1-1; i > 2*n1_1*n2_1-2*n1_1; i = i -2){
      nRefAr[i][2] = nrefcot[2];
    }
    for (int i = 2*n1_1-1 ; i<2*(n1_1)*(n2_1); i+=2*(n1_1)){
      nRefAr[i][1] = nrefcot[1];
    }
    for (int i = 0; i<2*(n1_1)*(n2_1); i+=2*(n1_1)){
      nRefAr[i][1] = nrefcot[3];
    }
  }
}

void creation_fichier_maillage(float a, float b, float c, float d, 
                               int n1, int n2, int t, 
                              int nrefdom, int *nrefcot){
  // Création du fichier en mode écriture
  FILE *meshfile = fopen("meshfile.txt", "w");
  if (meshfile == NULL) {
    perror("Erreur lors de l'ouverture du fichier");
    return;
  }

  // On écrit le nombre de noeuds totaux : n = n1*n2
  fprintf(meshfile, "%d\n", n1*n2);

  // On stocke les variables n1-1, n2-1 etc pour ne pas devoir les calculer à chaque fois
  int n1_1 = n1 - 1;
  int n2_1 = n2 - 1;

  // On écrit les coordonnées des noeuds géométriques
  for (int j=0; j<n2; j++){
    for (int i=0; i<n1; i++){
      float absc = a + i * ((b - a) / n1_1);
      float ordo = c + j * ((d - c) / n2_1);
      fprintf(meshfile, "%f %f\n", absc, ordo); 
    }
  }

  if (t == 1) { // Cas des quadrangles
    int m = (n1_1)*(n2_1); // nombre d'éléments
    int p = 4;             // nombre de noeuds par élément
    int q = 4;             //  nombre d'arêtes par élément

    // On crée le tableau nRefAr
    int** nRefAr = alloctabi(m,q);
    etiqAr(t ,n1_1, n2_1, nrefdom, nrefcot, m, q, nRefAr);

    
    fprintf(meshfile, "%d %d %d %d\n", m, t, p, q);
    int k =0 ;
      for (int j = 0; j<n2-1; j++){
	for (int i = 0; i<n1-1;  i++) {
	  fprintf(meshfile, "%d %d %d %d ",
		  n1*j + i + 2,
		  n1*(j+1) + i + 2,
		  n1*(j+1) + i +1,
		  n1*j + i + 1);
	  fprintf(meshfile, "%d %d %d %d\n", nRefAr[k][0], nRefAr[k][1],  nRefAr[k][2], nRefAr[k][3]);
	  k++;
	}
      }
      matrix_free_int(nRefAr);
    }

  if (t == 2) { // Dans le cas des triangles, on les traite deux par deux
     int m = 2*(n1_1)*(n2_1); // nombre d'éléments
     int p = 3;             // nombre de noeuds par élément
     int q = 3;             //  nombre d'arêtes par élément

     // On crée le tableau nRefAr
    int** nRefAr = alloctabi(2*(n1_1)*(n2_1),3);
    etiqAr(t ,n1_1, n2_1, nrefdom, nrefcot, m, q,  nRefAr);
    
    
      fprintf(meshfile, "%d %d %d %d\n", m, t, p, q);
      int k=0;
      for (int j = 0;  j<n2-1; j++){
	for (int i = 0;  i<n1-1; i++) {
	  fprintf(meshfile, "%d %d %d ",
		  n1*j + i + 2,
		  n1*(j+1) + i + 1,
		  n1*j + i +1);
	  fprintf(meshfile, "%d %d %d\n", nRefAr[k][0], nRefAr[k][1],  nRefAr[k][2]);
	  fprintf(meshfile, "%d %d %d ",
		  n1*(j+1) + i + 1,
		  n1*j + i + 2,
		  n1*(j+1) + i +2);
	  fprintf(meshfile, "%d %d %d\n", nRefAr[k+1][0], nRefAr[k+1][1],  nRefAr[k+1][2]);
	  k+=2;
	}
      }
      matrix_free_int(nRefAr);
    }
  
  fclose(meshfile);
}

int lecfima(char *ficmai, int *ptypel, int *pnbtng, float ***pcoord, int *pnbtel, int ***pngnel, int *pnbeel, int *pnbaret, float **pnRefAr) {
    // Ouverture du fichier en mode lecture
    FILE* meshfile = fopen(ficmai, "r");
    if (meshfile == NULL) {
        printf("Erreur : Impossible d'ouvrir le fichier %s\n", ficmai);
        return -1;
    }

    // Lecture du nombre total de noeuds
    fscanf(meshfile, "%d", pnbtng); 

    // Allocation d'un tableau pcoord pour obtenir les coordonnées des points du fichier 
    *pcoord = alloctab(*pnbtng, 2);
    // Remplissage du tableau pcoord
    for (int i = 0; i< *pnbtng; i++) fscanf(meshfile,"%f %f", &(*pcoord)[i][0], &(*pcoord)[i][1]);

    // Lecture de m, t, p, q
    fscanf(meshfile, "%d %d %d %d", pnbtel,ptypel, pnbeel, pnbaret);

    // Allocation mémoire pour le numéro global des noeuds de chaque élément + numéro de référence de chaque arête de chaque élément
    *pngnel = alloctabi(*pnbtel, *pnbeel);
    *pnRefAr = alloctab(*pnbtel, *pnbaret);

    // Remplissage des deux tableaux (en "même" temps):
    for (int i =0; i<*pnbtel; i++){
        for (int j = 0; j<*pnbeel; j++){
            fscanf(meshfile, "%d", &(*pngnel)[i][j]);
        }
        for (int j =0; j< *pnbaret; j++){
            fscanf(meshfile, "%d", pnRefAr[i][j]);
        }
    }

    fclose(meshfile);
    return 0;
}

// FONCTIONS DU TP2 A
int nppquad(int t){
    if (t==1) return 9;
    return 3;
}

void ppquad(int t, int q, float *omegak, float **xk){
    float a = 1./36; float b = 1./9; float c = 1./2; float d = 1./6;
    switch(t){
        case 1: // Quadrangles
            for (int k=0; k<4; k++) omegak[k]=a;
            for (int k=4; k<8; k++) omegak[k]=b;
            omegak[8]= 4*b;
            xk[0][0]=1; xk[0][1]=0;
            xk[1][0]=1; xk[1][1]=1;
            xk[2][0]=0; xk[2][1]=1;
            xk[3][0]=0; xk[3][1]=0;
            xk[4][0]=1; xk[4][1]=c;
            xk[5][0]=c; xk[5][1]=1;
            xk[6][0]=0; xk[6][1]=c;
            xk[7][0]=c; xk[7][1]=0;
            xk[8][0]=c; xk[8][1]=c;
            break;
        case 2: // Triangles
            for (int k= 0; k<3; k++) omegak[k]= d;
            xk[0][0]=c; xk[0][1]=c;
            xk[1][0]=0; xk[1][1]=c;
            xk[2][0]=c; xk[2][1]=0;
            break;
        case 3: // Segment
            omegak[0]= d; omegak[1]= d; omegak[2]= 2*d;
            xk[0][0]=1;
            xk[1][0]=0;
            xk[2][0]=c;
            break;
        }
}

void calFbase(int t, float *x_hat, float *w_hat){
    switch(t){
        case 1: //Quadrangle
            w_hat[0] = x_hat[0]*(1-x_hat[1]);
            w_hat[1] = x_hat[0]*x_hat[1];
            w_hat[2] = (1-x_hat[0])*x_hat[1];
            w_hat[3] = (1-x_hat[0])*(1-x_hat[1]);
            break;
        case 2: // Triangle   
            w_hat[0] = x_hat[0]*(1-x_hat[1]);
            w_hat[1] = (1-x_hat[0])*x_hat[1];
            w_hat[2] = 1 - x_hat[1]-x_hat[2];
            break;
        case 3: // Segment
            w_hat[0] = x_hat[1];
            w_hat[1] = 1-x_hat[1];
            break;
        }
}

void calDerFbase (int t, float *pt, float **valDerFbx){
    switch(t){
        case 1: // Quadrangle
            valDerFbx[0][0] = 1-pt[1]  ;valDerFbx[0][1] = -pt[0];
            valDerFbx[1][0] = pt[1]    ;valDerFbx[1][1] = pt[0];
            valDerFbx[2][0] = -pt[1]   ;valDerFbx[2][1] = 1-pt[0];
            valDerFbx[3][0] = -1+pt[1] ;valDerFbx[3][1] = pt[0]-1;
            break;
        case 2: // Triangle
            valDerFbx[0][0] = 1  ;valDerFbx[0][1] = 0;
            valDerFbx[1][0] = 0  ;valDerFbx[1][1] = 1;
            valDerFbx[2][0] = -1 ;valDerFbx[2][1] = -1;
            break;
        case 3: // Segment
            valDerFbx[0][0] = 0   ;valDerFbx[0][1] = 1;
            valDerFbx[1][0] = -1  ;valDerFbx[1][1] = 0;
            break;
    }
}

void transFk(float *Fbasexhat, float **ak, int p, float *Fkxhat){
    Fkxhat[0]=0; Fkxhat[1]=0;  
    for (int i=0; i<p; i++){
        Fkxhat[0]+= ak[i][0]*Fbasexhat[i];
        Fkxhat[1]+= ak[i][1]*Fbasexhat[i];
    }
}

void matJacob(float **ak, int d, int p, float **Derfbasexhat, float **JacobFk){
    // On initialise tout à 0
    JacobFk[0][0] = 0; JacobFk[0][1] = 0;
    JacobFk[1][0] = 0; JacobFk[1][1] = 0;
    if (d==1){
        for(int i=0; i<p; i++){
            JacobFk[0][0] += ak[i][0]* Derfbasexhat[i][0];
            JacobFk[1][0] += ak[i][1]* Derfbasexhat[i][0];
        }
    }
    if (d==2){
        for(int i=0; i<p; i++){
            JacobFk[0][0] += ak[i][0]* Derfbasexhat[i][0];
            JacobFk[0][1] += ak[i][0]* Derfbasexhat[i][1];
            JacobFk[1][0] += ak[i][1]* Derfbasexhat[i][0];
            JacobFk[1][1] += ak[i][1]* Derfbasexhat[i][1];
        }
    }
}

float invertM2x2 (float **M, float **Minv){
    float det = M[0][0]*M[1][1]- M[0][1]*M[1][0];
    if(fabsf(det)<precision){
        printf("matrice non inversible");
        exit(1);
    }
    float inv_det = 1./det;
    Minv[0][0] = M[1][1]*inv_det;
    Minv[1][1] = M[0][0]*inv_det;
    Minv[1][0] = -M[1][0]*inv_det;
    Minv[0][1] = -M[0][1]*inv_det;
    return det; 
}

void numNaret(int typel, int numaret, int *numNeAret){
    switch(typel){
        case 1: // Quadrangles
            numNeAret[0] = numaret;
            numNeAret[1] = (numaret+1)%4;
            break;
        case 2: // Triangle
            numNeAret[0] = numaret;
            numNeAret[1] = (numaret+1)%3;
            break;
    }
}

void selectPts(int nb, const int pts_nb[], float *coordSet[], float *select_coord[]){
    for(int i=0; i<nb; i++) {
        select_coord[i] = coordSet[pts_nb[i]-1]; 
    }
}

// FONCTION TP 2B
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


// Paramètres :
// @t         : type de l'élément (quadrangle, triangle, rectangle)
// @q         : nombre de points de quadrature
// @nbneel    : nombre de noeuds par élément 
// @coorEl    : coordonnées des noeuds sur l'élément
// @omegak    : poids des points de quadrature 
// @xhat      : coordonnées des points de quadrature
// @vecElem   : calcul et stock de l'intégrale issue de W
// @matElem   : calcul et stock de l'intégrale issu de WW et ADWDW

void intElem( int t,
            int q,
            int nbneel,
            float **coorEl,
            float *omegak,
            float **xhat,
            float **matElem,
            float *vecElem){

    //------- Initialisation des variables -------------------------------------------------------------------
    int d =2; // dimension
    if (t==1) d=1;

    float wk[nbneel];                               // calcul de Fkx^
    float *Fkxhat = malloc(d * sizeof(float));                               // vecteur qui contiendra 

    float **dwk; dwk = alloctab(nbneel, d);         // 
    float **dwkEl; dwkEl = alloctab(nbneel,d);      // pour ADWDW

    float **dFkxhat; dFkxhat = alloctab(2,2);         // Jacobienne de Fk(x^)

    float **dFkxhatinv; dFkxhatinv = alloctab(2,2); // inverse de la matrice D(Fk(x^))
    float det;                                    // déterminant de celle-ci

    // Initialisation des coefficients qui nous serviront à remplir matElem et vecElem
    float cofvar_W, cofvar_WW;
    float **cofvar_ADWDW;  cofvar_ADWDW = alloctab(2,2);
    float eltdif; 

    //----------- Boucle pour le remplissage de matElem et vecElem----------------------------------------------
    for(int k =0; k<q; k++){

        // Calcul de eltdif
        calDerFbase(t, xhat[k], dwk);
        matJacob(dwk, d, q, coorEl[k], dFkxhat);
        det = invertM2x2(dFkxhat, dFkxhatinv);
        eltdif = omegak[k]*fabsf(det);

        // Calcul de cofvar_w
        calFbase(t, xhat[k], wk);
        transFk(dwk, wk, q, Fkxhat); // pas sûre de celle là, il faut les ak 
        cofvar_W = FOMEGA(Fkxhat);

        // contribution avec la procédure W
        W(nbneel, wk, eltdif, cofvar_W, vecElem);

        // calcul de cofvar_ww
        cofvar_WW = a00(Fkxhat);

        // contribution avec la procédure WW
        WW(nbneel, wk, eltdif, cofvar_WW, matElem); 

        // calcul de cofvar_ADWDW
        cofvar_ADWDW[0][0]= a00(Fkxhat);
        cofvar_ADWDW[0][1]= a12(Fkxhat);
        cofvar_ADWDW[1][0]= a12(Fkxhat); 
        cofvar_ADWDW[1][1]= a22(Fkxhat);

        for (int i=0; i<nbneel; i++){
            dwkEl[i][0] = dwk[i][0]*dFkxhatinv[0][0] + dwk[i][1]*dFkxhatinv[1][0];
            dwkEl[i][1] = dwk[i][1]*dFkxhatinv[0][1] + dwk[i][1]*dFkxhatinv[1][1];
        }

        // contribution avec la procédure ADWDW 
        ADWDW(dwkEl, cofvar_ADWDW, eltdif, nbneel, matElem);

        // libération de la mémoire 
        matrix_free(dwk); matrix_free(dwkEl); free(wk); free(Fkxhat); 
    }
}

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
                int **nRefArEl,
                float **MatElem, float *SMbrElem, int *NuDElem, float *uDElem){
    // Calcul de q
    int q;
    q = nppquad(typeEl);

    // calcul de xhat et omega k
    float **xhat;  xhat = alloctab(q,2);
    float omegak[q];

    ppquad(typeEl ,q, omegak, xhat);

    // à l'intérieur du domaine :
    intElem(typeEl, q, nbneel, coorEl, omegak, xhat, MatElem, SMbrElem);

}

void impCalEl(int K, int typEl, int nbneel, float **MatElem, float *SMbrElem,
              int *NuDElem, float *uDElem) {
/************************************************************************
  Imprime les resultats de la matrice et du second membre elementaires
  ainsi que les conditions Dirichlet en chaque noeud
  et les valeurs des conditions Dirichlet non homogene
 
*** Arguments *** 
   K        : Numero de l'element
   typEl    : Numero de type de l'element
   nbneel   : Nombre de noeuds de l'element
   MatElem  : Matrice elementaire de dimensions (nbneel,nbneel)
   SMbrElem : Second membre elementaire de dimension nbneel
   NuDElem  : Tableau de reperage des noeuds porteurs de conditions de Dirichlet
   uDElem   : Tableau des valeurs de blocage
              pour les noeuds Dirichlet non homogene
************************************************************************/
  int i, j;
  printf("\n");
  printf(" ELEMENT=%3d    DE TYPE=%5d    NB NOEUDS=%2d\n", K,typEl,nbneel);
  printf(" NuDElem   uDElem    SMbrElem    MatElem\n");
  for (i=0; i < nbneel; i++) {
    printf(" %6d %10.4e %10.4e", NuDElem[i],uDElem[i],SMbrElem[i]);
    for (j=0; j <= i; j++) { printf(" %10.4e", MatElem[i][j]); }
    printf("\n");
  }
}
