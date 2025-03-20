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
                int ***nRefArEl,
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

    for(int i=0; i<nbneel; i++){
        nRefArEl
    }

}

        if (ref_edg_K[i] == ref_interior) { flag_categorized = 1; } // Most current case : edge is inside the domain
        else {
            flag_categorized = 0;                   // edge isn't categorized yet
            vertices_Edge(i + 1, t, edg_nod_ind);
            for (j = 0; j < N_NOD_EDG; j++) edg_nod_ind[j]--;
            // We retrieve indices (2) ? Shouldn't do it directly in vertices edges and change the function's name

            for (j = 0; j < n_Dh & flag_categorized == 0; j++) {      // Dirichlet homogeneous edge
                if (ref_edg_K[i] == ref_Dh[j]) {
                    for (k = 0; k < N_NOD_EDG; k++) nodes_D[edg_nod_ind[k]] = 0;
                    flag_categorized = 1;
                }
            }
            for (j = 0; j < n_Dnh & flag_categorized == 0; j++) {     // Dirichlet non-homogeneous edge
                if (ref_edg_K[i] == ref_Dnh[j]) {
                    for (k = 0; k < N_NOD_EDG; k++) {
                        nodes_D[edg_nod_ind[k]] = -1;
                        uD_aK[edg_nod_ind[k]] = uD(a_K[edg_nod_ind[k]]);
                    }
                    flag_categorized = 1;
                }
            }
            for (j = 0; j < n_NF & flag_categorized == 0; j++) {      // Neumann or Fourier edge
                if (ref_edg_K[i] == ref_NF[j]) {
                    selectPts(N_NOD_EDG, edg_nod_ind, (float **) a_K, edg_nodes_coords);/// @warning lost "const" on a_K
                    wp_quad(SEG_TYPE, x_quad_hat_edg, weights_edg);

                    // Re-intialize A_K_edg and l_K_edg to 0;
                    for (k = 0; k < N_NOD_EDG; k++) {
                        l_K_edg[k] = 0;
                        for (l = 0; l < N_NOD_EDG; l++) A_K_edg[k][l] = 0;
                    }

                    // Update results (we have the contribution for each node at the end of the edge (2 nodes)
                    intEdge(n_quad_pts_edg, (const float **) edg_nodes_coords, (const float **) x_quad_hat_edg,
                            weights_edg, A_K_edg, l_K_edg);
                    for (k = 0; k < N_NOD_EDG; k++) {
                        l_K[edg_nod_ind[k]] += l_K_edg[k];
                        for (l = 0; l < N_NOD_EDG; l++) A_K[edg_nod_ind[k]][edg_nod_ind[l]] += A_K_edg[k][l];
                    }
                    flag_categorized = 1;
                }
            }
        }