
































#include<stdio.h>





void star_mesh_base_alg(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,/*double** vl,*/ unsigned int len/*,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n*/){
     //unsigned int egde_id=0;
     struct edge;
     typedef struct edge_lst{
          struct edge* edg;
          struct edge_lst* next;
     
     } edg_lt;
     typedef struct node{
          unsigned int num;
          edg_lt* edge_list;
          struct node* next;
          
     } node;
     typedef struct edge{
          node* n1;
          node* n2;
          double val;
     
     } edge;
     
     node* node_hd=(node*) malloc(1*sizeof(node)); node_hd->num=col[0]; node* curr_node=node_hd;
     node** node_arr=(node**) malloc(col[len-1]*sizeof(node*)); node_arr[0]=curr_node;
     for(unsigned int i=col[0]+1;i<=col[len-1];i++){//assumption here is that node numbering is sequential without skipping the numbers;
          curr_node->next=(node*)malloc(1*sizeof(node));
          curr_node=curr_node->next;
          curr_node->num=i;
          curr_node->edge_list=NULL;
          node_arr[i-1]=curr_node;
          
     }
     curr_node->next=NULL;
     edge *curr_edge=NULL; edg_lt* curr_edge_list;
     //unsigned int curr_col=0;
     for(unsigned int i=0;i<len;i++){
          if(col[i]<row[i]){
               curr_edge=(edge*) malloc(1*sizeof(edge));
               curr_edge->n1=node_arr[col[i]-1];
               curr_edge->n2=node_arr[row[i]-1];
               curr_edge->val=val[i];
               curr_edge_list=(edg_lt*)malloc(1*sizeof(edg_lt));
               curr_edge_list->edg=curr_edge;
               /*if(curr_edge->n1->edge_list==NULL){
                    curr_edge->n1->edge_list=curr_edge_list;
                    curr_edge_list->next=NULL;
                    
               }
               else{
                    curr_edge_list->next=curr_edge->n1->edge_list;
                    curr_edge->n1->edge_list=curr_edge_list;
               
               }*/
               curr_edge_list->next=curr_edge->n1->edge_list;
               curr_edge->n1->edge_list=curr_edge_list;
               curr_edge_list=(edg_lt*)malloc(1*sizeof(edg_lt));
               curr_edge_list->edg=curr_edge;
               curr_edge_list->next=curr_edge->n2->edge_list;
               curr_edge->n2->edge_list=curr_edge_list;
          
          
          }
          
          
     
     }
     
     curr_node=node_hd;
     node** nds_td_indp=(node **)malloc(col[len-1]*sizeof(node*));    
     for(unsigned int i=0;i<col[len-1];i++){
          nds_td_indp[i]=NULL;
     }
     
     
     
     
     while(curr_node!=NULL){
          nds_td_indp[curr_node->num-1]=curr_node;
          curr_node=curr_node->next;
     }
     
     
     for(unsigned int i=0;i<col[len-1];i++){
          if(nds_td_indp[i]!=NULL){
               curr_edge_list=nds_td_indp[i]->edge_list;
               curr_node=nds_td_indp[i];
               while(curr_edge_list!=NULL){
                    nds_td_indp[curr_edge_list->edg->n1->num-1]=NULL;//one of the nodes (n1 or n2) is the i^{th} node; but to avoid if else branch, set NULL to both nodes, and then restore current node using curr_node;
                    nds_td_indp[curr_edge_list->edg->n2->num-1]=NULL;
                    curr_edge_list=curr_edge_list->next;
                    
               
               }
               nds_td_indp[i]=curr_node;
          
          }
     
     }
     
     
     
     
     
     unsigned int nds_td_max=nds_td[nds_n-1],nds_n0=0;
     for(unsigned int i=col[len-1]-1;i>0;){
          while(i>0 && nds_td_indp[i]==NULL) i--;
          if(i>0 && nds_td_indp[i]->num>nds_td_max){
               nds_td_indp[i]=NULL;
               i--;
               
          
          }
          else if(i>0 && nds_td_indp[i]->num<=nds_td_max){
               for(unsigned int j=0;j<=i;j++){
                    if(nds_td_indp[j]!=NULL){
                         nds_td_indp[nds_n0]=nds_td_indp[j];
                         nds_n0++;
                    
                    }
               
               }
               break;
          
          }
     
     }
     
     
     
     
     double* nd_sum=(double*)malloc(nds_n0*sizeof(double));
     edg_lt** nd_edg_lt=(edg_lt**)malloc(2*nds_n0*sizeof(edg_lt*));
     for(unsigned int i=0;i<nds_n0;i++){
          nd_sum[i]=0;
          for(nd_edg_lt[2*i]=nds_td_indp[i]->edge_list; nd_edg_lt[2*i]->next!=NULL; nd_edg_lt[2*i]=nd_edg_lt[2*i]->next){
               nd_sum[i]+=nd_edg_lt[2*i]->edge->val;
          
          }
          
          for(nd_edg_lt[2*i]=nds_td_indp[i]->edge_list; nd_edg_lt[2*i]->next!=NULL; nd_edg_lt[2*i]=nd_edg_lt[2*i]->next){
               for(nd_edg_lt[2*i+1]=nd_edg_lt[2*i]->next; nd_edg_lt[2*i+1]!=NULL;nd_edg_lt[2*i+1]=nd_edg_lt[2*i+1]->next){
                    curr_edge=(edge*) malloc(1*sizeof(edge));
                    unsigned int n1=(nd_edg_lt[2*i+1]->edge->n1==nds_td_indp[i]->num)?(nd_edg_lt[2*i+1]->edge->n2):(nd_edg_lt[2*i+1]->edge->n1);
                    unsigned int n2=(nd_edg_lt[2*i]->edge->n1==nds_td_indp[i]->num)?(nd_edg_lt[2*i]->edge->n2):(nd_edg_lt[2*i]->edge->n1);
                    if(n1>n2){
                         unsigned int n3=n1; n1=n2; n2=n1;
                    }
                    curr_edge->n1=node_arr[n1-1];
                    curr_edge->n2=node_arr[n2-1];
                    curr_edge->val=1./((1./nd_edg_lt[2*i]->edge-val)*(1./nd_edg_lt[2*i+1]->edge->val)*nd_sum[i]);
                    curr_edge_list=curr_edge->n1->edge_lst;
                    edg_lt* curr_edge_list1, *curr_edg_lst_prev;
                    while(curr_edge_list!=NULL; curr_edge_list->edge->n1->num > n1){
                         curr_edg_lst_prev=curr_edge_list;
                         curr_edge_list=curr_edge_list->next;
                    
                    }
                    if(curr_edge_list==NULL){
                    
                    }
                    else if(curr_edge_list->edge->n1->num==n1){
                         
                         whi
                         curr_edge_list1=curr_edge_list->next;
                         
                    
                    }
                    
                    
                    if(curr_edge_list->edge->n1->num == n1){
                         
                    
                    }
                    for(unsigned
                    
                         
                    
               
               
               
               }
          
          }
          
     
     
     }
     
     
     
     
     
     rw=(unsigned int**)0; cl=rw; rw=cl;
     
         

}

