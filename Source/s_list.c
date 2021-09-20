#include "hpc.h"
#include "hpc_linked_list.h"
/* Basic implementation of a singly linked list datastructure, which is used to
   store colums/rows of sparse matrices. */


/* creates a singly linked list
    returns the head-node of the list */
node* create_slist(){
    node* head = malloc(sizeof(node));
    head->data = 0.0;
    head->ind = 0;
    head->next = NULL;
    return (node *) malloc(sizeof(node));
}

/* adds a node after tail node,
   retruns the new tail node */
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

    if (head->next == NULL)
    {
        printf("Linked List is empty\n");
        return;
    }

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

    if (head->next == NULL)
    {
        printf("Linked List is empty\n");
        return;
    }

    node* node = head;
    while(node->next != NULL){
        node = node->next;
        printf("%d ", (int)node->ind);
    }
    printf("\n");
}

void freeList(struct node* head)
{
   struct node* tmp;

   while (head != NULL)
    {
       tmp = head;
       head = head->next;
       free(tmp);
    }

}