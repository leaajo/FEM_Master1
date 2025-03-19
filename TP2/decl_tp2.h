#include <stdio.h>
#include <stdlib.h>

extern const float precision;


int  **alloctabi(int dim1, int dim2);

float  **alloctab(int dim1, int dim2);

void matrix_free_int(int** ppMatrix);

void matrix_free(float** ppMatrix);

void etiqAr(int t, int n1_1, int n2_1, int nrefdom, int *nrefcot, int nbtel, int nbaret, int **nRefAr);

void creation_fichier_maillage(float a, float b, float c, float d, 
                               int n1, int n2, int t, 
                              int nrefdom, int *nrefcot);
                              
int lecfima(char *ficmai, int *ptypel, int *pnbtng, float ***pcoord, 
            int *pnbtel, int ***pngnel, int *pnbeel, int *pnbaret, 
            int ***pnRefAr);

// Fonctions du TP2a
int nppquad(int t);
void ppquad(int t, int q, float *omegak, float **xk);
void calFbase(int t, float *x_hat, float *w_hat);
void calDerFbase (int t, float *pt, float **valDerFbx);
void transFk(float *Fbasexhat, float **ak, int p, float *Fkxhat);
void matJacob(float **ak, int d, int p, float **Derfbasexhat, float **JacobFk);
float invertM2x2 (float **M, float **Minv);
void numNaret(int typel, int numaret, int *numNeAret);
void selectPts(int nb, const int pts_nb[], float *coordSet[], float *select_coord[]);

// Fonction du TP2b
float a00 (float *x);
float a11 (float *x);
float a12 (float *x);
float a22 (float *x);
float bN (float *x);
float FOMEGA (float *x);
float FN (float *x);
float UD (float*x);

void WW (int nbneel, float *fctbas, float eltdif, float cofvar, float **matelem);
void W (int nbneel, float *fctbas, float eltdif, float cofvar, float *matelem);
void ADWDW(float **derfctbas, float **cofvar, float eltdif, int nbneel, float **matelem);

void intElem( int t,
            int q,
            int nbneel,
            float **coorEl,
            float *omegak,
            float **xhat,
            float **matElem,
            float *vecElem);
            
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
                int ***nRefArEl,
                float **MatElem, float *SMbrElem, int *NuDElem, float *uDElem);

void impCalEl(int K, int typEl, int nbneel, float **MatElem, float *SMbrElem,
              int *NuDElem, float *uDElem);