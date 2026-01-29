
































//#include<iostream>
#include<stdlib.h>

struct nb_lst{
     unsigned int id;
     struct nb_lst* next;


}
struct nd_graph{
     unsigned int id;
     unsigned int nb_sz;
     struct nb_lst head;
     struct nb_lst tail;
     struct nb_lst* nbrhd;

};

void nd_schdlr(unsigned short int ** adj_bln, unsigned int n_nodes){
     static struct nd_graph*=(struct nd_graph*) malloc(sizeof(nd_grph)*n_nodes);
     if(adj_bln!=NULL){
          for(unsigned int i=0;i<n_nodes;i++){//TODO keep updated when node is deleted;
               nd_graph[i].id=i+1;//Matlab indexing; identificator of the node is its number;
               
          
          }
          
     
     }
     else{
     
     }
     
     

     

}



