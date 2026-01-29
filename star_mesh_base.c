































                                
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
     else{//TODO col can change between *nds_n and *len;
     
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
               else if(row[kk]<nds_td[j]){//TODO here if col change from *nds_n to *len;
                    kk++;
               
               }
               else j++;
          
          }
     
     }
     
     unsigned int* cnt_tmps=(unsigned int*)malloc((*nds_n)*sizeof(unsigned int));
     unsigned int tot_sz=0;
     mesh_bdrs[0]=0;
     for(unsigned int i=0;i<*nds_n;i++){
          mesh_bdrs[i+1]=mesh_bdrs[i]+((idx_r[2*i+1]-idx_r[2*i])*(idx_r[2*i+1]-idx_r[2*i]-1)+0.)/2.;//new mesh boundary strats where ends previous one plus mesh size (number of edges);
          
     }
     tot_sz=2*mesh_bdrs[*nds_n]+*len;//each edge will containe 2 nodes, but pair is recorded twice in adjacency matrix; but for now record data only once;
     (*rw)=(unsigned int*)malloc(tot_sz*sizeof(unsigned int));
     (*cl)=(unsigned int*)malloc(tot_sz*sizeof(unsigned int));
     (*vl)=(double*)malloc(tot_sz*sizeof(double));
     //#pragma omp parallel for num_threads((*nds_n)<100?(*nds_n):100)
     for(unsigned int i=0;i<*nds_n;i++){
          cnds_sum[i]=0; cnt_tmps[i]=mesh_bdrs[i];
          for(unsigned int j=idx_r[2*i];j<idx_r[2*i+1];j++){
               cnds_sum[i]+=val[j];
          
          }
          for(unsigned int j=idx_r[2*i];j<idx_r[2*i+1]-1;j++){
               for(unsigned int k=j+1;k<idx_r[2*i+1];k++){
                    (*cl)[cnt_tmps[i]]=row[j];//keep column number smaller than row number;
                    (*rw)[cnt_tmps[i]]=row[k];
                    (*vl)[cnt_tmps[i]]=1./((1./val[j])*(1./val[k])*cnds_sum[i]);
                    (cnt_tmps[i])++;
               
               }
          
          }
          
     }
     ////////////////////////////////////////////////////
     for(unsigned int i=0;i<mesh_bdrs[*nds_n]-1;i++){
          for(unsigned int j=i+1;j<mesh_bdrs[*nds_n];j++){
               if((*cl)[i]==(*cl)[j] && (*rw)[i]!=0 && (*rw)[i]==(*rw)[j]){
                    (*vl)[i]+=(*vl)[j];//parallel conductances be added;
                    (*rw)[j]=0;//mark added parallel edge;
                    
               }
          }
     
     }//at this point parrallel conductances of new meshes should summed;
     unsigned int cnt_tmp_g=mesh_bdrs[*nds_n]; unsigned char fl_dct=0;
     for(unsigned int i=0;i<idx_len_nb;i++){
          for(unsigned int j=idx_r_nb[2*i], cnt_lcl=idx_r_nb[2*i+1]-1; j<idx_r_nb[2*i+1]; j++,cnt_lcl--){//cnt_lcl inverse indexing needed to go from biggest row index to smallest;//TODO j is not used in this loop;
               if(col[cnt_lcl]<row[cnt_lcl]){//to avoid duplicates; && row[j]!=0 should be satisfyed;
                    fl_dct=0;
                    for(unsigned int k=0;k<mesh_bdrs[*nds_n];k++){
                         if((*rw)[k]==row[cnt_lcl] && (*cl)[k]==col[cnt_lcl]){
                              (*vl)[k]+=val[cnt_lcl];
                              row[k]=0; fl_dct=1;
                              break;
                         
                         }
                    
                    
                    }
                    if(fl_dct==0){
                         (*cl)[cnt_tmp_g]=col[cnt_lcl];
                         (*rw)[cnt_tmp_g]=row[cnt_lcl]; (*vl)[cnt_tmp_g]=val[cnt_lcl]; cnt_tmp_g++;
                         
                    
                    } 
               
               }
               else if(row[cnt_lcl]!=0 && col[cnt_lcl]>row[cnt_lcl]) break;//to avoid duplicates; go to new column; cnt_lcl inverse indexing alows to break here, since all begger indexes should be visited in earlier part of loop;
          
          }
     
     }//at this point edges should be merged;
  
     
     unsigned int cnt_nz=0;
     for(unsigned int i=0;i<cnt_tmp_g;){
          while((*rw)[i]==0){
               i++;
          
          }
          while((*rw)[i]!=0 && i<cnt_tmp_g){
               (*cl)[cnt_nz]=(*cl)[i]; (*rw)[cnt_nz]=(*rw)[i]; (*vl)[cnt_nz]=(*vl)[i];
               cnt_nz++; i++;
          
          }
          
     
     }//shifts;
     
     for(unsigned int i=0;i<cnt_nz;i++){
          (*cl)[cnt_nz+i]=(*rw)[i]; (*rw)[cnt_nz+i]=(*cl)[i]; (*vl)[cnt_nz+i]=(*vl)[i];
     
     }//to comply with adjacency matrix format, repeat values with switched rows and columns;
     (*rw)=(unsigned int*)realloc((*rw),2*cnt_nz*sizeof(unsigned int));
     (*cl)=(unsigned int*)realloc((*cl),2*cnt_nz*sizeof(unsigned int));
     (*vl)=(double*)realloc((*vl),2*cnt_nz*sizeof(double));
     *ln=2*cnt_nz;
     
     for(unsigned int i=0;i<*nds_n;i++){
          //free(edge_nodes[i]);
          //free(edge_weights[i]);
     
     }
     
     free(cnds_sum); /*free(idx_l);*/ free(idx_r); free(idx_r_nb); /*free(edge_nodes);*/ /*free(edge_weights);*/
     free(mesh_bdrs); free(cnt_tmps);
     //free(indx_r);
     
     
}

#include"node_schr1.c"

