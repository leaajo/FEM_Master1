#include <stdio.h>
#include <stdlib.h>

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

int lecfima(char *ficmai, int *ptypel, int *pnbtng, float ***pcoord, int *pnbtel, int ***pngnel, int *pnbeel, int *pnbaret, int ***pnRefAr) {
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
    *pnRefAr = alloctabi(*pnbtel, *pnbaret);

    // Remplissage des deux tableaux (en "même" temps):
    for (int i =0; i<*pnbtel; i++){
        for (int j = 0; j<*pnbeel; j++){
            fscanf(meshfile, "%d", &(*pngnel)[i][j]);
        }
        for (int j =0; j< *pnbaret; j++){
            fscanf(meshfile, "%d", &(*pnRefAr)[i][j]);
        }
    }

    fclose(meshfile);
    return 0;

}
