#include <stdlib.h>
#include <stdio.h>
const float precision = 1e-7;

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
            for (int k= 0; k<3; k++) omegakk[k]= d;
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

void calDerFbase (int t, float *pt, float **valDerFbx,){
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

void selectPts(int nb, const int pts_nb[], float *coordSet[], float *select_coord[])
{
    for(int i=0; i<nb; i++) {
        select_coord[i] = coordSet[pts_nb[i]-1]; 
    }
}