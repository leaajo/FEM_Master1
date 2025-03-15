#include <stdio.h>
#include <stdlib.h>
#include "decl.h"


int main() {
    // Paramètres pour tester la fonction
    float a = 0.0, b = 1.0, c = 0.0, d = 1.0;
    int n1 = 5, n2 = 5;
    int t = 1; // Quadrangles
    int nrefdom = 0;
    int nrefcot[] = {1, 2, 3, 4};

    printf("Avant le code.\n");
    // Appel de la fonction
    creation_fichier_maillage(a, b, c, d, n1, n2, t, nrefdom, nrefcot);

    printf("Fichier meshfile.txt généré avec succès.\n");

    // Déclarations des variables nécessaires
    char ficmai[] = "car3x3t_3.txt";  
    int ptypel = 0, pnbtng = 0, pnbtel = 0, pnbeel = 0, pnbaret = 0;
    float **pcoord = NULL;
    int **pngnel = NULL, **pnRefAr = NULL;

    // Appel de la fonction
    int result = lecfima(ficmai, &ptypel, &pnbtng, &pcoord, &pnbtel, &pngnel, &pnbeel, &pnbaret, &pnRefAr);

    // Vérification du résultat
    if (result == 0) {
        printf("Lecture réussie\n");
        printf("pnbtng = %d\n", pnbtng);
        printf("Coordonnées des noeuds:\n");
        for (int i = 0; i < pnbtng; i++) {
            printf("Noeud %d: (%f, %f)\n", i, pcoord[i][0], pcoord[i][1]);
        }
    } else {
        printf("Erreur lors de la lecture du fichier\n");
    }
    return 0;
}
