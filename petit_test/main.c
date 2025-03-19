#include <stdio.h>
#include <stdlib.h>

int main() {
    int pnbeel = 2;
    int pnbtel = 3;  // Assurez-vous que pnbtel ne dépasse pas 2
    int pngnel[3][2] = {
        {2, 3},
        {3, 2},
        {1, 4}
    };
    
    // Allocation mémoire pour le tableau NB
    int* NB = malloc(pnbeel * sizeof(int));

    if (NB == NULL) {
        // Si malloc échoue, afficher un message d'erreur et quitter
        printf("Memory allocation failed\n");
        return 1;
    }

    // Boucle à travers les indices de pngnel
    for (int K = 0; K < pnbtel; K++) {  // Limité à pnbtel qui est 2
        // On affiche les valeurs de pngnel pour vérifier que l'on accède correctement
        for (int i = 0; i < pnbeel; i++) {
            // Affichage de l'élément actuel
            printf("K = %d, i = %d, pngnel[K][i] = %d\n", K, i, pngnel[K][i]);
            // Affectation de la valeur à NB
            NB[i] = pngnel[K][i];
            // Affichage de NB[i]
            printf("%d ", NB[i]);
        }
        printf("\n"); // Ajouter un saut de ligne après chaque itération de K
    }

    // Libération de la mémoire allouée
    free(NB);

    return 0;
}
