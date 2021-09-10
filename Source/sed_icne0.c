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
    nzmax = A->nzmax;

    /* helper variables for storing intermediate results */
    double l_ii;
    index nzmax_l = n+1; // number of nenzro elements of L + d, is set to n+1 initially
    double* a_j = malloc(n*sizeof(double)); // stores rows/cols of A intermediately
    double* l_j = malloc((n-1)*sizeof(double)); // NOTE evl linked list
    node* u_k;

    /* array of linked lists for storing columns of L */
    node** cols_head = malloc((n-1)*sizeof(node*)); // first nodes of rows
    node** cols_tail = malloc((n-1)*sizeof(node*)); // last nodes of rows
    for(index i=0; i<n-1; i++) { // initialize linked lists
        cols_head[i] = create_slist();
        cols_tail[i] = cols_head[i]; 
    }

    double* d = malloc(n*sizeof(double));
    d[0] = x[0];

    for(index j = 1; j < n; j++){

        for(index i =0; i<n; i++){ a_j[i] = 0;} // set a_j to zero
        
        for(index i = ind[j]; i < ind[j+1]; i++){ //iterate through column j
            
            /* calucualte l_ii = (a_i, a_i) + alpha */
            l_ii += pow(x[i],2);
            l_ii += alpha;

            a_j[ind[i]] = x[i]; // extract collumn_j
        }
        a_j[j] = x[j];

        /* calculate all nonzero inner products: l_ij = (a_i,a_j) */
        for(index i = 0; i < j; i++){
            for(index jj = ind[i]; jj<ind[i+1]; jj++ ){ // iterate through col i

                l_j[i] += a_j[ind[jj]] * x[jj];
            }
            l_j[i] += x[i] * a_j[i];

            if(l_j[i] == 0) { continue; } // storing only nonzero inner-products

            //store l_j in a linked list datastructure
            cols_tail[i] = append_node(cols_tail[i], l_j[i] /*data */, j /* row index */);
            nzmax_l++; // increment number of nonzero elements
        }

        for(index k = 0; k < j; k++){

            if(l_j[k] == 0) { continue; } //proceed only for nonzero inner products

            u_k = cols_head[k]; // u_k = k-th column of L
            l_j[k] = l_j[k] / d[k];
            
            while(u_k->next != NULL){
                l_j[k] = l_j[k] - u_k->data * l_j[u_k->ind];
            }
        }
        d[j] = l_ii;
    }

    /* store d and  L in a SED-Matrix  */
    L = sed_alloc(n, nzmax_l, 1);
    
    // store d on the diagonal of L
    for(index i = 0; i < n; i++){
        L->x[i] = d[i];
    }

    index count = n+1;
    for(index i=0; i<n-1; i++){
        
        node* current = cols_head[i];

        while(current->next != NULL){

        }
    }
    return (1);
}

