/* --------- linked list datastructure ------------- */
/* 
Author: Sebastian Acerbi 
*/

// Structure for nodes in a singly linked list
typedef struct node
{
    double data; 
    index ind; //index associate with the data
    struct node *next;
}node;

node* create_slist();
node* append_node(node* tail, double data, index ind);
void print_list_data(node* head);
void print_list_ind(node* head);
void freeList(struct node* head);