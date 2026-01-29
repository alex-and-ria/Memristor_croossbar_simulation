































                                
#include"node_schr.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


struct edge_1{
     unsigned int nn1;
     double wgth;
     struct edge_1* next;
     

};

struct edge_0{
     unsigned int nn0;
     struct edge_1* e1_head;
     unsigned int sz;
     struct edge_0* next;
     


};


void col_fnd(unsigned int *col, unsigned int *len, unsigned int nd_td_i,unsigned int* idx_l, unsigned int* idx_r){
     unsigned int base_idx, offset;
     base_idx=((*len)+1)/2; offset=base_idx;
     for(unsigned int j=0; j<log2(*len); j++){//binary search for fragment in "col" where information about node nds_td[i] stored;
          if(col[base_idx]<nd_td_i){
               offset=(offset+1)/2;
               base_idx+=offset;
          
          }
          else if (col[base_idx]>nd_td_i){
               offset=(offset+1)/2;
               base_idx-=offset;
          
          }
          else{
               break;
          
          }
     
     }
     
     
     *idx_l=base_idx; *idx_r=base_idx;
     while((*idx_l)<(*len) && col[*idx_l]>=nd_td_i)
          (*idx_l)--;
     if((*idx_l)>(*len)) (*idx_l)=0;//if overflow of unsigned int;
     else (*idx_l)++;//else point to first occurence of node location in "col";
     while((*idx_r)<(*len) && col[*idx_r]<=nd_td_i) (*idx_r)++;//points to first occurence bigger then nds_td[i], or to *len;

}

void col_fnd_full(unsigned int *col, unsigned int *len, unsigned int** idx_r, unsigned int* idx_len){
     (*idx_r)[0]=0; unsigned int curr_col=col[0], cnt_tmp=1;
     for(unsigned int i=1;i<*len;i++){
          if(col[i]>curr_col){
               curr_col=col[i];
               (*idx_r)[cnt_tmp]=i;
          
          }
     }
     (*idx_r)[cnt_tmp]=*len;
     cnt_tmp++;
     *idx_len=cnt_tmp;
     
     (*idx_r)=(unsigned int*) realloc((*idx_r),(*idx_len)*sizeof(unsigned int));
     
     

}

//col_find_fw()//keep idx_l, binary serarch from end to idx_l; TODO;

void node_analyzer(unsigned int* row, unsigned int *col, double* val, unsigned int *len, unsigned int* nds_td, unsigned int* nds_n){//adj_m starts count from node 1 (Matlab indexing); it is represented in sparse Matlab format via [row,col,val]; col are sorted; nds_td should count from one (Matlab idexing), and be sorted; it is a list of nodes to delete; check nodes and delete one of them from nds_td if two are neighbours;

     unsigned int *idx_l, *idx_r, idxl,idxr;  idx_l=&idxl; idx_r=&idxr;
     for(unsigned int i=0;i<*nds_n;i++){
          if(nds_td[i]==0) continue;
          /*base_idx=((*len)+1)/2; offset=base_idx;
          for(unsigned int j=0; j<log2(*len); j++){//binary search for fragment in "col" where information about node nds_td[i] stored;
               if(col[base_idx]<nds_td[i]){
                    offset=(offset+1)/2;
                    base_idx+=offset;
               
               }
               else if (col[base_idx]>nds_td[i]){
                    offset=(offset+1)/2;
                    base_idx-=offset;
               
               }
               else{
                    break;
               
               }
          
          }
          
          idx_l=idx_r=base_idx;
          while(idx_l<(*len) && col[idx_l]>=nds_td[i])
               idx_l--;
          if(idx_l>(*len)) idx_l=0;//if overflow of unsigned int;
          else idx_l++;//else point to first occurence of node location in "col";
          while(idx_r<(*len) && col[idx_r]<=nds_td[i]) idx_r++;//points to first occurence bigger then nds_td[i], or to *len;*/
          col_fnd(col,len, nds_td[i],idx_l,idx_r);
          
          for(unsigned int j=*idx_l,k=i+1; j<*idx_r&&k<*nds_n;){
               if(row[j]<nds_td[k]){
                    j++;
               }
               else if(row[j]>nds_td[k]){
                    k++;
               }
               else{
                    nds_td[k]=0;
                    j++; k++;
               }
          }
          
     }
     
     /*unsigned int i_idx=0,j_idx=0;
     for(unsigned int i=0;i<*nds_n;i++){
          if(nds_td[i]==0) continue;
          for(unsigned int j=i+1;j<*nds_n;j++){
               if(nds_td[j]==0) continue;
               i_idx=nds_td[i]-1; j_idx=nds_td[j]-1;
               if(adj_m[i_idx][j_idx]!=0){//if nodes are neighbours, delete them from the list of nds_td (nodes that can be deleted in parralel);
                    nds_td[j]=0;
                    
               
               }
          
          }
     
     }
     */
     
     //prune zeros;
     

}

/*
void star_mesh_prl1(unsigned int* adj_m_row, unsigned int* dim,nds_n[i],unsigned int** adj_new){



}
*/


/*void star_mesh_base(unsigned int *row, unsigned int *col, double *val, unsigned int *len, unsigned int *nds_td, unsigned int *nds_n){
     unsigned int cnds_sum=0;
     for(

}*/




void star_mesh_base(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val, double** vl,unsigned int *len, unsigned int *nds_td, unsigned int *nds_n/*unsigned int** adj_m, unsigned int* dim*/){
     double *cnds_sum; cnds_sum=(double*) malloc((*nds_n)*sizeof(double));
     unsigned int *idx_l; idx_l=(unsigned int*) malloc((*nds_n)*sizeof(unsigned int));
     unsigned int *idx_r; idx_r=(unsigned int*) malloc((*nds_n)*sizeof(unsigned int));
     //unsigned int** edge_nodes; edge_nodes=(unsigned int**) malloc((*nds_n)*sizeof(unsigned int*));
     double** edge_weights; edge_weights=(double**) malloc((*nds_n)*sizeof(double*));
     //unsigned int sz_nd_d=0,new_sz=0;
     unsigned int* indx_r=(*unsigned int) malloc((*len)*sizeof(unsigned int)); unsigned int idx_len;
     col_fnd_full(col,len,&indx_r,&idx_len);
     
     
     
     
     for(unsigned int i=1,j=0;i<idx_len && j<*nds_n;i++){
          if(col[indx_r[i]-1]==nds_td[j]){
               idx_l[j]=indx_r[i-1];
               idx_r[j]=indx_r[i];
               j++;
          }
          else{
               for(unsigned int kk=indx_r[i-1], ll=0; kk<indx_r[i] && ll<*nds_n; ll++){
                    if(row[kk]==nds_td[ll]){
                         row[kk]=0;
                         kk++;
                    
                    }
               
               }//clearing all neighbours that should be deleted;
            
          
          }
     
     }

     for(unsigned int i=0;i<*nds_n;i++){
          //col_fnd(col,len,nds_td[i],&(idx_l[i]),&(idx_r[i]));
          unsigned int mesh_sz=((idx_r[i]-idx_l[i])*(idx_r[i]-idx_l[i]-1)+0.)/2.; unsigned int cnt_tmp=0;
          //sz_nd_d+=(idx_r[i]-idx_l[i]); new_sz+=mesh_sz;
          //edge_nodes[i]=(unsigned int*) malloc(mesh_sz*sizeof(unsigned int));
          edge_weights[i]=(double*) malloc(mesh_sz*sizeof(double));
          cnds_sum[i]=0; 
          for(unsigned int j=idx_l[i];j<idx_r[i];j++){
               cnds_sum[i]+=val[j];
          
          }
          for(unsigned int j=idx_l[i];j<idx_r[i]-1;j++){
               for(unsigned int k=j+1;k<idx_r[i];k++){
                    edge_weights[i][cnt_tmp]=1./((1./val[j])*(1./val[k])*cnds_sum[i]);
                    cnt_tmp++;
                    
               
               }
          
          }
          
          
     
          /*cnds_sum=0;
          for(unsigned int j=0;j<*dim;j++) cnds_sum+=adj_m[nds_td[i]-1][j];
          ///////////////////////////////////////////////////parrallel loop;
          for(unsigned int j=0;j<*dim-1;j++){
               for(unsigned int k=j+1;k<*dim;k++){
                    if(adj_m[nds_td[i]-1][j]!=0 && adj_m[nds_td[i]-1][k]!=0){
                         adj_m[j][k]=1./((1./adj_m[nds_td[i]-1][j])*(1./adj_m[nds_td[i]-1][k])*cnds_sum);
                    
                    
                    }
                    
                    
               
               }
          
          }
          //////////////////////////////////////////*/
     
     
     
     }
     
     
     
     /*unsigned int cnt_tmp0=0, len_1=0;
     unsigned int col_prev=0;
     for(unsigned int i=0,j=0; i<*len && j<*nds_n;){
          if(col[i]==nds_td[j]){
               cnt_tmp0+=((idx_r[j]-idx_l[j])*(idx_r[j]-idx_l[j]-1)+0.)/2.;
               j=idx_r[j];
               i+=(idx_r[j]-idx_l[j]);
          
          }
          else{
               cnt_tmp0++;
               i++;
          
          }
          
     
     }*/
     
     /*new_sz=new_sz+(*len)-sz_nd_d;//2*(total new edges)+(edge that is not effected);
     (*cl)=(unsigned int*) malloc(new_sz*sizeof(unsigned int));
     (*rw)=(unsigned int*) malloc(new_sz*sizeof(unsigned int));
     (*vl)=(double*) malloc(new_sz*sizeof(double));
     unsigned int cnt_tmp_g=0,cnt_tmp_l=0;
     for(unsigned int i=0,j=0; i<*len && j<*nds_n;){
          if(col[i]==nds_td[j]){
               cnt_tmp_l=0;
               for(unsigned int k=idx_l[j];k<idx_r[j]-1;k++){
                    for(unsigned int ll=k+1;ll<idx_r[j];ll++){
                         cl[cnt_tmp_g]=row[k]; rw[cnt_tmp_g]=row[ll]; vl[cnt_tmp_g]=edge_weights[j][cnt_tmp_l];
                         cl[cnt_tmp_g+1]=row[ll]; rw[cnt_tmp_g+1]=row[k]; vl[cnt_tmp_g+1]=edge_weights[j][cnt_tmp_l];
                         cnt_tmp_g+=2; cnt_tmp_l++;
                    
                    }
               
               }
               j++;
               i=idx_r[j];
               
          
          }
          else{
               
          
          }
     
     }*/
     ///////////////////////////////////////////////////////////////////////////////////////
     /*
     struct edge_0 *e0_head=NULL, *e0_curr=NULL, *e0_prev=NULL, *e0_curr_seq=NULL;
     struct edge_1 *e1_head=NULL, *e1_curr=NULL, *e1_prev=NULL;
     unsigned int curr_cnt=0;
     unsigned int curr_col=0;
     
   
     for(unsigned int i=0,j=0; i<*len && j<*nds_n;){
          if(col[i]==nds_td[j]){//this node deleted; include corresponding mesh;
               curr_cnt=0;
               for(unsigned int kk=idx_l[j];k<idx_r[j]-1;k++){
                    for(unsigned int ll=kk+1;ll<idx_r[j];ll++){//row[kk] is left node (smaller index), row[ll] is right node (later index);
                         e1_prev=NULL; e0_prev=NULL;
                         for(e0_curr=e0_head;e0_curr!=NULL;e0_curr=e0_curr->next){
                              if(e0_curr->nn0==row[kk]){//left node exists;
                                   for(e1_curr=e0_curr->e1_head;e1_curr!=NULL;e1_curr=e1_curr->next){
                                        if(e1_curr->nn1==row[ll]){//right node already exists; means parallel resistors; conductance should be added;
                                             e1_curr->wgth+=edge_weights[j][curr_cnt];
                                             curr_cnt++;
                                             break;
                                        
                                        }
                                        else if(e1_curr->nn1>row[ll]){//need insert new node before current;
                                             if(e1_prev==NULL){
                                                  e1_prev=(struct edge_1*) malloc(1*sizeof(struct edge_1));
                                                  e1_prev->next=e0_curr->e1_head;
                                                  e0_curr->e1_head=e1_prev;
                                                  e1_curr=e0_curr->e1_head;
                                                  
                                             }
                                             else{
                                                  e1_prev->next=(struct edge_1*) malloc(1*sizeof(struct edge_1));
                                                  e1_prev->next->next=e1_curr;
                                                  e1_curr=e1_prev->next;
                                                  
                                             }
                                             e0_curr->sz++;
                                             e1_curr->nn1=row[ll];
                                             e1_curr->wght=edge_weights[j][curr_cnt];
                                             curr_cnt++;
                                             break;
                                             
                                        }
                                        e1_prev=e1_curr;
                                   
                                   }
                                   if(e1_curr==NULL){//e0_curr!=NULL, hence at least e1_curr=e0_curr->e1_head was executed, and it is some edge (not NULL); all other endpoints are smaller (add to the tail);
                                        e1_prev->next=(struct edge_1*) malloc(1*sizeof(struct edge_1));
                                        e1_curr=e1_prev->next;
                                        e0_curr->sz++;
                                        e1_curr->nn1=row[ll];
                                        e1_curr->wght=edge_weights[j][curr_cnt];
                                        curr_cnt++;
                              
                                   
                                   }//at this point edge should be established;
                                   break;
                                   
                                   
                              
                              }
                              else if(e0_curr->nn0>row[kk]){//need to set up new node before the current;
                                   if(e0_prev==NULL){//e0_curr==e0_head;
                                        e0_prev=(struct edge_0*) malloc(1*sizeof(struct edge_0));
                                        e0_prev->next=e0_head;
                                        e0_head=e0_prev;
                                        e0_curr=e0_head;
                                        
                                        
                                   
                                   }
                                   else{
                                        e0_prev->next=(struct edge_0*) malloc(1*sizeof(struct edge_0));
                                        e0_prev->next->next=e0_curr;
                                        e0_curr=e0_prev->next;
                                   
                                   }
                                   e0_curr->nn0=row[kk];
                                   e0_curr->sz=1;
                                   e0_curr->e1_head=(struct edge_1*) malloc(1*sizeof(struct edge_1));
                                   e1_curr=e0_curr->e1_head;
                                   e1_curr->nn1=row[ll];
                                   e1_curr->wght=edge_weights[j][curr_cnt];
                                   curr_cnt++;
                                   e1_curr->next=NULL;
                                   break;
                                   
                              
                              }
                              e0_prev=e0_curr;
                         
                         }
                         if(e0_curr==NULL){
                              if(e0_head=NULL){//e0_prev==NULL;
                                   e0_head=(struct edge_0*) malloc(1*sizeof(struct edge_0));
                                   e0_curr=e0_head;
                                   
                              }
                              else{
                                   e0_prev->next=(struct edge_0*) malloc(1*sizeof(struct edge_0));
                                   e0_curr=e0_prev->next;
                              
                              }
                              e0_curr->nn0=row[kk];
                              e0_curr->sz=1;
                              e0_curr->e1_head=(struct edge_1*) malloc(1*sizeof(struct edge_1));
                              e1_curr=e0_curr->e1_head;
                              e1_curr->nn1=row[ll];
                              e1_curr->wght=edge_weights[j][curr_cnt];
                              curr_cnt++;
                              e1_curr->next=NULL;
                              
                         
                         }//at this point edge should be established;
                         
                    
                    }
               
               }

               j++;
               i=idx_r[j];


          }
          else{
               unsigned int kk=0; unsigned int l_node,r_node;
               e0_curr_seq=e0_head;
               e1_prev=NULL; e0_prev=NULL;
               curr_col=col[i];
               while(curr_col=col[i]){
                    for(kk<*nds_n;){
                         if(row[i]==nds_td[kk]){//row[i] is one of nodes to delete, skip it;
                              i++; kk++;
                              break;
                              
                         
                         
                         
                         
                         }
                         else if(row[i]<nds_td[kk]){
                              if(col[i]<row[i]){//place row[i] and col[i] in structure;
                                   l_node=col[i]; r_node=row[i];
                              }
                              else{
                                   l_node=row[i]; r_node=col[i];
                              }
                              for(;e0_curr_seq!=NULL;e0_prev=e0_curr_seq,e0_curr_seq->next){
                                   if(e0_curr_seq->nn0==l_node){
                                        
                                   
                                   }
                                   else if(e0_curr_seq->nn0>l_node){
                                   
                                   }
                              
                              }
                              
                              
                         
                         }
                         
                    
                    }
                    if(kk==*nds_n){
                    
                    }
               
               }
               


          }
     }*/
     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     for(unsigned int i=0;i<*nds_n;i++){
          //free(edge_nodes[i]);
          free(edge_weights[i]);
     
     }
     
     free(cnds_sum); free(idx_l); free(idx_r); /*free(edge_nodes);*/ free(edge_weights);
     

//else if == nds_td[j] calc mesh_zs else(>) prune +sz;
}
/*
void node_dispatcher(unsigned int** adj_m, unsigned int* dim, unsigned int* nds_td, unsigned int* nds_n){

     for(unsigned int i=0;i<*nds_n;i++){//parrallel loop;
     /////////////////////////////////////////
          star_mesh(adj_m[nds_n[i]-1],dim,nds_n[i],unsigned int** adj_new);
    
          for(unsigned int j=0;j<*dim;j++){//delete old connections to the node; nodes in nds_td are not neighbours, so conncetions can be deleted parallely without interferring with other nodes;
               adj_m[nds_td[i]-1][j]=0;
               adj_m[j][nds_td[i]-1]=0;
          
          }
    /////////////////////////////////////////
    

     
     }
     

     for(unsigned int i=0;i<*nds_n;i++){//megre; parallel loop;
     
     }
     
     
     //fold rest;
     star_mesh_prl2(unsigned int** adj_m, unsigned int* dim,unsigned int* nds_td, unsigned int* nds_n);



}
*/





int tst_pnt(double** adj_m,unsigned int* dim){
     printf("q\n");
     for(unsigned int i=0;i<*dim*(*dim);i++) printf("%g ",*(adj_m[0]+i));
     //unsigned int q; scanf("%d",&q);
     /*printf("%d ",adj_m[0][0]);
     for(unsigned int i=0;i<*dim;i++){
          for(unsigned int j=0;j<*dim;j++){
               printf("%d",adj_m[i][j]);
          
          }
          printf("\nsz=%d",*dim);
     
     }*/
     printf("\nsz=%d",*dim);




     
     return 32;


}
