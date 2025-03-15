#include <stdio.h>
#include <stdlib.h>

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

