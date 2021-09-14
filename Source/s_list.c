#include "hpc.h"
#include "hpc_linked_list.h"
/* Basic implementation of a singly linked list datastructure, which is used to
   store colums/rows of sparse matrices. */


/* creates a singly linked list
    returns the head-node of the list */
node* create_slist(){
    node* head = NULL;
    return (node *) malloc(sizeof(node));
}

/* adds a node after current node,
   retruns the new current node */
node* append_node(node* tail, double data, index ind){
    
    node* new = malloc(sizeof(node));
    new->data = data;
    new->ind = ind;
    new->next = NULL;
    tail->next = new;
    return new;
}

/* given the head node, the merhode prints
   all data mambers of the linked list */
void print_list_data(node* head){

    node* node = head;
    while(node->next != NULL){
        node = node->next;
        printf("%2.2f ", node->data);
    }
    printf("\n");
}

/* given the head node, the merhode prints
   all index mambers of the linked list */
void print_list_ind(node*head){

    node* node = head;
    while(node->next != NULL){
        node = node->next;
        printf("%d ", (int)node->ind);
    }
    printf("\n");
}