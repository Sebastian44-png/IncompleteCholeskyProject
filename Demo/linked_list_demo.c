#include "hpc.h"
#include "hpc_linked_list.h"

int main(){

    node* head_list1 = create_slist();

    node* current = head_list1;

    for (index i=0; i<10; i++){
        current = append_node(current, i, (double)i+1);
    }
    
    node* head_list2 = create_slist();

    current = head_list2;

    for (index i=0; i<10; i++){
        current = append_node(current, i*i, (double)i+1);
    }

    node** cols = malloc(sizeof(node*)*2);
    cols[0] = head_list1;
    cols[1] = head_list2;

    for(int i =0; i< 2; i++){
        print_list_data(cols[i]);
    }

    free(head_list1);
    free(head_list2);
}