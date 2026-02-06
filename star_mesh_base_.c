































                                
#include"node_schr.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"



void nds_sprt(unsigned int *col, unsigned int *len, unsigned int** idx_r_td, unsigned int** idx_r_nb, unsigned int *idx_nb_len, unsigned int* nds_td, unsigned int* nds_n){
     unsigned int curr_col=col[0], cnt_nb=0;
     //unsigned int i=1;
     if(col[0]==nds_td[0]){
          (*idx_r_td)[0]=0;
     }
     else{
          (*idx_r_nb)[0]=0;
     }
     
     for(unsigned int i=1,j=0;i<*len;i++){
          if(curr_col!=col[i]){
               curr_col=col[i];
               if(j<*nds_n && col[i-1]==nds_td[j]){
                    (*idx_r_td)[2*j+1]=i;
                    j++;
               
               }
               else{
                    (*idx_r_nb)[2*cnt_nb+1]=i;
                    cnt_nb++;
               
               }
               if(j<*nds_n && col[i]==nds_td[j]){
                    (*idx_r_td)[2*j]=i;
               
               }
               else{
                    (*idx_r_nb)[2*cnt_nb]=i;
                    
               
               }
               
          
          }
     
     }
     if(col[(*len)-1]==nds_td[(*nds_n)-1]){
          (*idx_r_td)[2*(*nds_n)-1]=*len;
     
     }
     else{
     
          (*idx_r_nb)[2*cnt_nb+1]=*len;
          cnt_nb++;
          
     
     }
     *idx_nb_len=cnt_nb;
     (*idx_r_nb)=(unsigned int*) realloc((*idx_r_nb),(2*cnt_nb)*sizeof(unsigned int));//realloc to free extra allocated memory;

} 




void star_mesh_base(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,double** vl, unsigned int *len,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n){
     double *cnds_sum; cnds_sum=(double*) malloc((*nds_n)*sizeof(double));
     unsigned int *idx_r; idx_r=(unsigned int*) malloc(2*(*nds_n)*sizeof(unsigned int));// unsigned int idx_len;//indexes for nodes to delete; will index up to *nds_n; should contain first and next after last index for each node to delete;
     
     unsigned int *mesh_bdrs; mesh_bdrs=(unsigned int*) malloc(((*nds_n)+1)*sizeof(unsigned int));
     unsigned int* idx_r_nb=(unsigned int*) malloc(2*(*len)*sizeof(unsigned int)); unsigned int idx_len_nb;//indexes for neighbouring and protected nodes;
     nds_sprt(col, len, &idx_r, &idx_r_nb, &idx_len_nb, nds_td, nds_n);
 //clearing all neighbours that that should be deleted, //(*now just all neighbours) that located after the j==*nds_n-1;
     for(unsigned int i=0; i<idx_len_nb;i++){
          for(unsigned int j=0,kk=idx_r_nb[2*i]; j<*nds_n && kk<idx_r_nb[2*i+1];){
               if(row[kk]==nds_td[j]){
                    row[kk]=0;
                    kk++;
                    j++;
                    
               
               }
               else if(row[kk]<nds_td[j]){
                    kk++;
               
               }
               else j++;
          
          }
     
     }
     
     unsigned int* cnt_tmps=(unsigned int*)malloc((*nds_n)*sizeof(unsigned int));
     mesh_bdrs[0]=0;
     for(unsigned int i=0;i<*nds_n;i++){
          mesh_bdrs[i+1]=mesh_bdrs[i]+((idx_r[2*i+1]-idx_r[2*i])*(idx_r[2*i+1]-idx_r[2*i]-1)+0.);//new mesh boundary strats where ends previous one plus mesh size (number of edges); each edge is recorded twice to comply with adjacency matrix requirenments;
          
     }
     (*rw)=(unsigned int*)malloc((mesh_bdrs[*nds_n]+((*len)-2*(*nds_n)))*sizeof(unsigned int));
     (*cl)=(unsigned int*)malloc((mesh_bdrs[*nds_n]+((*len)-2*(*nds_n)))*sizeof(unsigned int));
     (*vl)=(double*)malloc((mesh_bdrs[*nds_n]+((*len)-2*(*nds_n)))*sizeof(double));
     unsigned int* rw_buff=(unsigned int*)malloc(mesh_bdrs[*nds_n]*sizeof(unsigned int));
     unsigned int* cl_buff=(unsigned int*)malloc(mesh_bdrs[*nds_n]*sizeof(unsigned int));
     double* vl_buff=(double*)malloc(mesh_bdrs[*nds_n]*sizeof(double));
     unsigned int* rw_buff1=(unsigned int*)malloc(mesh_bdrs[*nds_n]*sizeof(unsigned int));
     unsigned int* cl_buff1=(unsigned int*)malloc(mesh_bdrs[*nds_n]*sizeof(unsigned int));
     double* vl_buff1=(double*)malloc(mesh_bdrs[*nds_n]*sizeof(double));
     //#pragma omp parallel for num_threads((*nds_n)<100?(*nds_n):100)
     for(unsigned int i=0;i<*nds_n;i++){
          cnds_sum[i]=0;
          for(unsigned int j=idx_r[2*i];j<idx_r[2*i+1];j++){
               cnds_sum[i]+=val[j];
          
          }
          for(unsigned int j=idx_r[2*i];j<idx_r[2*i+1]-1;j++){
               for(unsigned int k=j+1;k<idx_r[2*i+1];k++){
                    cnt_tmps[i]=mesh_bdrs[i]+(idx_r[2*i+1]-idx_r[2*i]-1)*(j-idx_r[2*i])+(k-idx_r[2*i]-1);
                    cl_buff[cnt_tmps[i]]=row[j];//column number smaller than row number;
                    rw_buff[cnt_tmps[i]]=row[k];
                    vl_buff[cnt_tmps[i]]=1./((1./val[j])*(1./val[k])*cnds_sum[i]);
                    cnt_tmps[i]=mesh_bdrs[i]+(idx_r[2*i+1]-idx_r[2*i]-1)*(k-idx_r[2*i])+(j-idx_r[2*i]);
                    cl_buff[cnt_tmps[i]]=row[k];//column number bigger than row number;
                    rw_buff[cnt_tmps[i]]=row[j];
                    vl_buff[cnt_tmps[i]]=1./((1./val[j])*(1./val[k])*cnds_sum[i]);
                    
               
               }
          
          }
          
     }
     ////////////////////////////////////////////////////
     unsigned int curr_i, idx_m0=0, idx_res=0; cnt_tmps[0]=0;
     struct nds_lst{unsigned int i; struct nds_lst* next;}; struct nds_lst* nds_lst_head=(struct nds_lst*) malloc(1*sizeof(struct nds_lst)); struct nds_lst* nl_prev=nds_lst_head, *nl_curr, *nl_m0; nds_lst_head->i=0; nds_lst_head->next=NULL;
     for(unsigned int i=1;i<*nds_n;i++){
          nl_curr=(struct nds_lst*) malloc(1*sizeof(struct nds_lst)); nl_curr->i=i; nl_prev->next=nl_curr; nl_prev=nl_curr;
          cnt_tmps[i]=mesh_bdrs[i];
     }
     nl_curr->next=NULL;
     //if(*nds_n==1) nds_lst_head->next=NULL; else nl_curr->next=NULL;
     nl_curr=nds_lst_head;
     while(nds_lst_head->next!=NULL){
          idx_m0=nds_lst_head->i;
          nl_prev=nds_lst_head;
          nl_m0=NULL;
          for(
          nl_curr=nds_lst_head->next
          ;nl_curr!=NULL;){
               curr_i=nl_curr->i;
               if(
               cl_buff[cnt_tmps[idx_m0]] == 
               cl_buff[cnt_tmps[curr_i]] 
               && rw_buff[cnt_tmps[idx_m0]] ==
                rw_buff[cnt_tmps[curr_i]]){
                    vl_buff[cnt_tmps[idx_m0]]+=vl_buff[cnt_tmps[curr_i]];
                    cnt_tmps[curr_i]++;
                    if(cnt_tmps[curr_i]==mesh_bdrs[curr_i+1]){
                         /*if(nl_curr==nds_lst_head){
                              nl_prev=nds_lst_head; nds_lst_head=nds_lst_head->next; free(nl_prev);
                              nl_prev=NULL;
                              
                         
                         }
                         else*/{
                              nl_prev->next=nl_curr->next; 
                              free(nl_curr); 
                              nl_curr=nl_prev->next;
                         
                         }
                    
                    }
                    else{
                         nl_prev=nl_curr; nl_curr=nl_curr->next;
                    
                    }
               }
               else{
                    if((cl_buff[cnt_tmps[idx_m0]] > cl_buff[cnt_tmps[curr_i]]) || ((cl_buff[cnt_tmps[idx_m0]] == cl_buff[cnt_tmps[curr_i]]) && (rw_buff[cnt_tmps[idx_m0]] > rw_buff[cnt_tmps[curr_i]]))){
                         idx_m0=curr_i;
                         nl_m0=nl_prev;
                         
               

                    }
                    nl_prev=nl_curr; nl_curr=nl_curr->next;
               
               }
               
          
          }
          cl_buff1[idx_res]=cl_buff[cnt_tmps[idx_m0]];
          rw_buff1[idx_res]=rw_buff[cnt_tmps[idx_m0]];
          vl_buff1[idx_res]=vl_buff[cnt_tmps[idx_m0]];
          idx_res++;
          cnt_tmps[idx_m0]++;
          if(cnt_tmps[idx_m0]==mesh_bdrs[idx_m0+1]){
               if(nl_m0==NULL){
                    nl_prev=nds_lst_head; nds_lst_head=nds_lst_head->next; free(nl_prev);
                    
               
               }
               else{
                    nl_prev=nl_m0->next; nl_m0->next=nl_m0->next->next; free(nl_prev);
               
               }
          
          }
          
          
         
     
     }//at this point sorted arrays should be merged, except nds_lst_head;
     for(unsigned int i=cnt_tmps[nds_lst_head->i];i<mesh_bdrs[(nds_lst_head->i)+1];i++){
          cl_buff1[idx_res]=cl_buff[i];
          rw_buff1[idx_res]=rw_buff[i];
          vl_buff1[idx_res]=vl_buff[i];
          idx_res++;
          
     
     }
     free(nds_lst_head);
     //at this point sorted array from deleted nodes should be merged themselves;
     ///////////////////////////////////////////////
     unsigned int idx_tmp=0; unsigned int i_tm=0,j_tm=0,k_tm=0;
     for(;k_tm<idx_res && i_tm<idx_len_nb;i_tm++){
          for(j_tm=idx_r_nb[2*i_tm];j_tm<idx_r_nb[2*i_tm+1] && k_tm<idx_res;){
               if(row[j_tm]==0){
                    j_tm++; continue;
                    
               
               }
               if(col[j_tm]<cl_buff1[k_tm] || (col[j_tm]==cl_buff1[k_tm] && row[j_tm]<rw_buff1[k_tm])){
                    (*cl)[idx_tmp]=col[j_tm];
                    (*rw)[idx_tmp]=row[j_tm];
                    (*vl)[idx_tmp]=val[j_tm];
                    idx_tmp++; j_tm++;
               
               }
               else if(col[j_tm]==cl_buff1[k_tm] && row[j_tm]==rw_buff1[k_tm]){
                    (*cl)[idx_tmp]=col[j_tm];
                    (*rw)[idx_tmp]=row[j_tm];
                    (*vl)[idx_tmp]=val[j_tm]+vl_buff1[k_tm];
                    idx_tmp++; j_tm++; k_tm++;
               
               }
               else{
                    (*cl)[idx_tmp]=cl_buff1[k_tm];
                    (*rw)[idx_tmp]=rw_buff1[k_tm];
                    (*vl)[idx_tmp]=vl_buff1[k_tm];
                    idx_tmp++; k_tm++;
               
               }
               
          
          }
     
     }
     if(k_tm==idx_res){
          i_tm--;
          for(;j_tm<idx_r_nb[2*i_tm+1];j_tm++){
               (*cl)[idx_tmp]=col[j_tm];
               (*rw)[idx_tmp]=row[j_tm];
               (*vl)[idx_tmp]=val[j_tm];
               idx_tmp++;
               
               
          
          }
          i_tm++;
          for(;i_tm<idx_len_nb;i_tm++){
               for(j_tm=idx_r_nb[2*i_tm];j_tm<idx_r_nb[2*i_tm+1];j_tm++){
                    (*cl)[idx_tmp]=col[j_tm];
                    (*rw)[idx_tmp]=row[j_tm];
                    (*vl)[idx_tmp]=val[j_tm];
                    idx_tmp++; 
                    
                    
               
               }
          
          } 
          
     
     }
     for(;k_tm<idx_res;k_tm++){
          (*cl)[idx_tmp]=cl_buff1[k_tm];
          (*rw)[idx_tmp]=rw_buff1[k_tm];
          (*vl)[idx_tmp]=vl_buff1[k_tm];
          idx_tmp++; k_tm++;
     
     }//at this point arrays from deletion of the node should be merged with original points;
     (*rw)=(unsigned int*)realloc((*rw),idx_tmp*sizeof(unsigned int));
     (*cl)=(unsigned int*)realloc((*cl),idx_tmp*sizeof(unsigned int));
     (*vl)=(double*)realloc((*vl),idx_tmp*sizeof(double));
     *ln=idx_tmp;
     
     
     free(cnds_sum); free(idx_r); free(idx_r_nb);
     free(mesh_bdrs); free(cnt_tmps);
     free(rw_buff); free(cl_buff); free(vl_buff);
     free(rw_buff1); free(cl_buff1); free(vl_buff1);
     








     
     
}

#include"node_schr2.c"
int tst_pnt(unsigned int dim) {return dim;}
