#include "hpc.h"
#include "hpc_linked_list.h"

/* perform the ICNE(0)-decomposition on an matrix stored
   in the sed format, the original matrix gets overwritten */

index sed_icne0(sed* A, double alpha, sed* L){

    double* x;
    index* ind;
    index n;
    index nzmax;

    x = A->x;
    ind = A->i;
    n = A->n;

    /* helper variables for storing intermediate results */
    double l_ii;
    nzmax = n; // number of nenzro elements of L + d, is set to n initially
    double* a_j = malloc(n*sizeof(double)); // stores rows/cols of A intermediately
    double* l_j = malloc((n-1)*sizeof(double));
    node* u_k;

    /* array of linked lists for storing columns of L */
    node** cols_head = malloc(n*sizeof(node*)); // first nodes of rows
    node** cols_tail = malloc(n*sizeof(node*)); // last nodes of rows
    for(index i=0; i<n; i++) { // initialize linked lists
        cols_head[i] = create_slist();
        cols_tail[i] = cols_head[i]; 
    }

    double* d = malloc(n*sizeof(double));
    d[0] = x[0];

    for(index i = 1; i < n; i++){

        for(index j =0; j<n; j++){ a_j[j] = 0;} // set a_j to zero
        l_ii = 0; //set l_ii to zero

        for(index p = ind[i]; p < ind[i+1]; p++){ //iterate through column j
                 
            l_ii += pow(x[p], 2);  // calucualte l_ii = (a_i, a_i) + alpha
            
            a_j[ind[p]] = x[p];    // extract collumn_j
        }
        a_j[i] = x[i]; // extract collumn_j

        l_ii += pow(x[i], 2);  // calucualte l_ii = (a_i, a_i) + alpha for the diagonal element
        l_ii += alpha;         // calucualte l_ii = (a_i, a_i) + alpha 
        
        /* calculate all nonzero inner products: l_ij = (a_i,a_j) */
        for(index j=0; j<n; j++){ l_j[j] = 0;} // set l_j to zero

        for(index j = 0; j < i; j++){
            
            for(index p = ind[j]; p<ind[j+1]; p++ ){ // iterate through col j

                l_j[j] += a_j[ind[p]] * x[p];
            }
            l_j[j] += x[j] * a_j[j];

            if(l_j[j] == 0) { continue; } // storing only nonzero inner-products
           
            //store l_j in a linked list datastructure
            cols_tail[j] = append_node(cols_tail[j], l_j[j] /*data */, i /* row index */);
            nzmax++; // increment number of nonzero elements
        }
        l_j[i] = l_ii;
        //printf("ld_%ld = ", i); print_buffer_double(l_j, n);

        for(index k = 0; k < i; k++){

            if(l_j[k] == 0) { continue; } //proceed only for nonzero inner products

            u_k = cols_head[k]; // u_k = k-th column of L
            //printf("u_%ld = ", k);print_list_data(u_k);

            // -2
            l_j[k] = l_j[k] / d[k];
            //printf("l_j = "); print_buffer_double(l_j, n);

            while(u_k->next != NULL){

                u_k = u_k->next;
                l_j[k] = l_j[k] - (u_k->data * l_j[u_k->ind]);
            }
            u_k->data = l_j[k];
        }
        d[i] = l_ii;
        //printf("l_j = "); print_buffer_double(l_j, n);
    }

    /* build sed_matrix containing L + d */
    sed_realloc(L, nzmax);
    // store d on the diagonal of L
    for(index i = 0; i < n; i++){
        L->x[i] = d[i];
    }

    index count = n+1;
    for(index i=0; i<n; i++){
        
        node* current_node = cols_head[i];

         // store column pointers
        L->i[i] = count;

        while(current_node->next != NULL){ // add colum i to L->x

            current_node = current_node->next;
            L->x[count] = current_node->data;
            L->i[count] = current_node->ind;

            count++;
        }
    }
    L->i[n] = count;

    return (1);
}

