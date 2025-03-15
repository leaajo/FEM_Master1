#include <stdio.h>
#include <stdlib.h>

// Read a mesh file and assign parameters values, coordinates, node global numbers, and edges references numbers
int read_mesh(char* input_mesh, int *type, int* n_nodes, float*** p_coords, int* n_elem, int*** p_gNb_node, int *n_nod_elem, int* n_edges, int*** p_refEdg)
{
    int i,j;

    // Opening meshfile in readmode
    FILE *mesh_f = fopen(input_mesh, "r");	// mode 'w' for reading
    if (mesh_f == (FILE*)NULL) {
        printf("Error when trying to open the input (mesh) file : '%s' in read_mesh function", input_mesh);
        perror("^^ IO - Error ^^");
        return -1;
    }

    fscanf(mesh_f,"%d",n_nodes);      // number of nodes (n)
    *p_coords = alloctab(*n_nodes,2);   // allocate a n x 2 tab for coordinates

    // fill coordinates x_i y_i
    for (i=0; i<*n_nodes; i++) fscanf(mesh_f,"%f %f",&(*p_coords)[i][0], &(*p_coords)[i][1]);

    // assign nb element (m), type (t), nb of geometrical node per element (p), nb of egdes per element (q);
    fscanf(mesh_f,"%d %d %d %d",n_elem,type,n_nod_elem,n_edges);

    // fill global node number reference for each element and edges reference number
    *p_gNb_node = alloctabI(*n_elem,*n_nod_elem);
    *p_refEdg = alloctabI(*n_elem,*n_edges);
    for (i=0; i<*n_elem; i++){
        for (j=0; j<*n_nod_elem; j++) fscanf(mesh_f,"%d ",&(*p_gNb_node)[i][j]);
        for (j=0; j<*n_edges; j++) fscanf(mesh_f,"%d ",&(*p_refEdg)[i][j]);
    }
    fclose(mesh_f); // Close mesh output file
    return 0;
}