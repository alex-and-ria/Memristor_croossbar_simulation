































                                
#include"node_schr.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

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

/*struct node_type{
     unsigned int n;
     struct edge_type* enges;

}

struct edge_type{
     struct node_type* node0;
     struct node_type* node1;
     double w;

}

*/

/*
#include <time.h>
struct timespec curr_time;
ret=clock_gettime(CLOCK_MONOTONIC,&curr_time);
uint64_t tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;


*/
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
               cnt_tmp++;
          
          }
     }
     (*idx_r)[cnt_tmp]=*len;
     cnt_tmp++;
     *idx_len=cnt_tmp;
     
     (*idx_r)=(unsigned int*) realloc((*idx_r),(*idx_len)*sizeof(unsigned int));
     
     

}

/*void col_fnd_dual(unsigned int *col, unsigned int *len, unsigned int** idx_r_td, unsigned int* idx_td_len, unsigned int** idx_r_nb, unsigned int *idx_nb_len, unsigned int* nds_td, unsigned int* nds_n){
     unsigned int curr_col=col[0], cnt_tmp_td=1;
     unsigned int cnt_tmp_nb=1;
     unsigned char fl_td=0,fl_nb=0;
     if(col[0]==nds_td[0]){
          (*idx_r_td)[0]=0;
          fl_td=1;
     
     }
     else{
          (*idx_r_nb)[0]=0;
          fl_nb=1;
     
     }
     for(unsigned int i=1, j=0;i<*len && j<*nds_n;i++){
          if(col[i]>curr_col){
               curr_col=col[i];
               if(col[i-1]==nds_td[j]){
                    if(fl_td==0)fl_td=1;
                    if(fl_nb==0) (*idx_r_nb)[0]=i;
                    (*idx_r_td)[cnt_tmp_td]=i;
                    cnt_tmp_td++;
                    j++;//if last node in nds_td is not last node in col, then at some point j==*nds_n, while i<*len;
                    
               }
               else{
                    if(fl_td==0) (*idx_r_td)[0]=i;
                    if(fl_nb==0) fl_nb=1;

                    (*idx_r_nb)[cnt_tmp_nb]=i;
                    cnt_tmp_nb++;
                    
               }
               
          }
     }
     
     if(col[(*len)-1]>nds_td[(*nds_n)-1]){
          for(unsigned int i=(*idx_r_td)[cnt_tmp_td-1]+1;i<(*len);i++){//record indexes of neighbourhouds that are unrecorded (have left) from previous loop;
               if(col[i]>curr_col){
                    curr_col=col[i];
                    (*idx_r_nb)[cnt_tmp_nb]=i;
                    cnt_tmp_nb++;
               
               }
               
          
          }
          (*idx_r_nb)[cnt_tmp_nb]=(*len);
          cnt_tmp_nb++;
     
     }
     else{
          (*idx_r_td)[cnt_tmp_td]=(*len);
          cnt_tmp_td++;
     
     }
     *idx_td_len=cnt_tmp_td;
     *idx_nb_len=cnt_tmp_nb;
     (*idx_r_td)=(unsigned int*) realloc((*idx_r_td),(*idx_td_len)*sizeof(unsigned int));//realloc to free extra allocated memory;
     (*idx_r_nb)=(unsigned int*) realloc((*idx_r_nb),(*idx_nb_len)*sizeof(unsigned int));
     

}*/

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


/*void node_analyzer(unsigned int* row, unsigned int *col, double* val, unsigned int *len, unsigned int* nds_td, unsigned int* nds_n, double* nb_koef){//adj_m starts count from node 1 (Matlab indexing); it is represented in sparse Matlab format via [row,col,val]; col are sorted; nds_td should count from one (Matlab idexing), and be sorted; it is a list of nodes to delete; check nodes and delete one of them from nds_td if two are neighbours;

     unsigned int *idx_l, *idx_r, idxl,idxr;  idx_l=&idxl; idx_r=&idxr;
     unsigned int cnt_zr=0;
     for(unsigned int i=0;i<*nds_n;i++){
          if(nds_td[i]==0) continue;
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
                    cnt_zr++;
               }
          }
          
     }
     (*nb_koef)=(cnt_zr+0.)/(*nds_n+0.);
}*/


void star_mesh_base(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,double** vl, unsigned int *len,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n){
     double *cnds_sum; cnds_sum=(double*) malloc((*nds_n)*sizeof(double));
     //unsigned int *idx_l; idx_l=(unsigned int*) malloc((*nds_n)*sizeof(unsigned int));
     unsigned int *idx_r; idx_r=(unsigned int*) malloc(2*(*nds_n)*sizeof(unsigned int));// unsigned int idx_len;//indexes for nodes to delete; will index up to *nds_n; should contain first and next after last index for each node to delete;
     
     unsigned int *mesh_bdrs; mesh_bdrs=(unsigned int*) malloc(((*nds_n)+1)*sizeof(unsigned int));
     //unsigned int** edge_nodes; edge_nodes=(unsigned int**) malloc((*nds_n)*sizeof(unsigned int*));
     //double** edge_weights; edge_weights=(double**) malloc((*nds_n)*sizeof(double*));
     //unsigned int sz_nd_d=0,new_sz=0;
     //unsigned int* indx_r=(unsigned int*) malloc((*len)*sizeof(unsigned int)); unsigned int idx_len;
     unsigned int* idx_r_nb=(unsigned int*) malloc(2*(*len)*sizeof(unsigned int)); unsigned int idx_len_nb;//indexes for neighbouring and protected nodes;
     //col_fnd_full(col,len,&indx_r,&idx_len);
     //col_fnd_dual(col, len, &idx_r, &idx_len, &idx_r_nb, &idx_len_nb, nds_td, nds_n);
     nds_sprt(col, len, &idx_r, &idx_r_nb, &idx_len_nb, nds_td, nds_n);
     
     /*for(unsigned int i=1,j=0;i<idx_len && j<*nds_n;i++){
          if(col[indx_r[i]-1]==nds_td[j]){
               idx_l[j]=indx_r[i-1];
               idx_r[j]=indx_r[i];
               j++;
          }
          else{
               for(unsigned int kk=indx_r[i-1], ll=0; kk<indx_r[i] && ll<*nds_n;){
                    if(row[kk]==nds_td[ll]){
                         row[kk]=0;
                         kk++;
                    
                    }
                    else if(row[kk]<nds_td[ll]) {kk++; ll--;}
                    ll++;
               
               }//clearing all neighbours that should be deleted;
            
          
          }
     
     }
     for(unsigned int i=idx_r[(*nds_n)-1]; i<*len;i++){
          for(unsigned int j=0;j<*nds_n;j++){
               if(row[i]==nds_td[j]){
                    if(i<*len-1 && row[i]>row[i+1]){//detecting approximate column boundaries to speed up computations;
                         row[i]=0;
                         break;
                    }
                    row[i]=0;
                    
               
               }
               else if(row[i]<nds_td[j]) break;
          
          }
     
     }*/ //clearing all neighbours that that should be deleted, //(*now just all neighbours) that located after the j==*nds_n-1;
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
          //col_fnd(col,len,nds_td[i],&(idx_l[i]),&(idx_r[i]));
          //unsigned int mesh_sz=((idx_r[i+1]-idx_r[i])*(idx_r[i+1]-idx_r[i]-1)+0.)/2.; unsigned int cnt_tmp=0;
          //mesh_szs[i]=mesh_sz;
          //sz_nd_d+=(idx_r[i]-idx_l[i]); new_sz+=mesh_sz;
          //edge_nodes[i]=(unsigned int*) malloc(mesh_sz*sizeof(unsigned int));
          //edge_weights[i]=(double*) malloc(mesh_sz*sizeof(double));
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
     /*unsigned int tot_sz=0;
     for(unsigned int i=0;i<*nds_n;i++){
          tot_sz+=2*mesh_szs[i];//each edge will containe 2 nodes, but pair is recorded twice in adjacency matrix; but for now record date only once;
     
     }*/
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
          for(unsigned int j=idx_r_nb[2*i], cnt_lcl=idx_r_nb[2*i+1]-1; j<idx_r_nb[2*i+1]; j++,cnt_lcl--){//cnt_lcl inverse indexing needed to go from biggest row index to smallest;
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
     
     /*unsigned int cnt_tmp_g=0,cnt_tmp=0;
     for(unsigned int i=1;i<idx_len_nb;i++){
          for(unsigned int j=idx_r_nb[i-1];j<idx_r_nb[i];j++){
               if(row[j]!=0 && col[j]<row[j]){
                    (*cl)[cnt_tmp_g]=col[j]; (*rw)[cnt_tmp_g]=row[j]; (*vl)[cnt_tmp_g]=val[j]; cnt_tmp_g++;//since this nodes have neighbours that is not affected by deletion, we should keep them; they will show up twice in adjacency matrix, for now we record that only once (only if col[j]<row[j], or, only if col[j]>row[j]);
               
               }
          
          }
          
     
     }*///recording nondeleted edges;
     
     
     /*mesh_szs=(unsigned int*) realloc(mesh_szs,((*nds_n)+1)*sizeof(unsigned int));
     unsigned int m_sz_prev=mesh_szs[0],m_sz_next=0; mesh_szs[0]=cnt_tmp_g;
     for(unsigned int i=0;i<(*nds_n);i++){
          m_sz_next=mesh_szs[i+1];
          mesh_szs[i+1]=m_sz_prev+mesh_szs[i];
          m_sz_prev=m_sz_next;
     
     }//now mesh_szs store not actual mesh sizes, but accumulated sum of mesh sizes to get relative initial indexes; we just reuse the array to reduce memory consumption;
     */
     /////////////////////////////////////////
     /*for(unsigned int i=1,j=0;i<idx_len && j<*nds_n;){
          if(col[indx_r[i]-1]==nds_td[j]){
               cnt_tmp=0;
               for(unsigned int kk=idx_l[j];kk<idx_r[j]-1;kk++){
                    for(unsigned int ll=kk+1;ll<idx_r[j];ll++){
                         (*cl)[cnt_tmp_g]=row[kk]; (*rw)[cnt_tmp_g]=row[ll]; (*vl)[cnt_tmp_g]=edge_weights[j][cnt_tmp]; cnt_tmp_g++; cnt_tmp++;
                      
                    }
               
               }
               i++; j++;
               
          }
          else{
               for(unsigned int kk=indx_r[i-1]; kk<indx_r[i]; kk++){//TODO duplicates;
                    if(row[kk]!=0 && col[kk]<row[kk]){//neither col[kk] nor row[kk] gets deleted, hence col has neightbour in row, and when we get to row, it should contain corresponding col as a neighbour, hece record it once at a time to avoid doubles;
                         (*cl)[cnt_tmp_g]=col[kk]; (*rw)[cnt_tmp_g]=row[kk]; (*vl)[cnt_tmp_g]=val[kk]; cnt_tmp_g++;
                    
                    }
               }
               i++;
          }
     
     }//get meshes merged, now need to merge doubles and add second pair to make adjacency matrix; at this point col should containe smaller node and row bigger node number in numbering scheme;
     
     unsigned int cnt_nz=0;
     for(unsigned int i=0;i<cnt_tmp_g;i++){
          for(unsigned int j=i+1;j<cnt_tmp_g;j++){
               if((*cl)[i]==(*cl)[j] && (*rw)[i]==(*rw)[j]){
                    (*cl)[j]=0; (*vl)[i]+=(*vl)[j]; //parrallel conductances are added;
                    if(cnt_nz==0) cnt_nz=j;
               
               }
          
          }
     
     }
     for(unsigned int i=cnt_nz,j=cnt_nz+1;j<cnt_tmp_g;){//shift;
          while((*cl)[j]!=0){
               (*cl)[i]=(*cl)[j]; (*rw)[i]=(*rw)[j]; (*vl)[i]=(*vl)[j];
               i++; j++;
               cnt_nz++;
          
          }
          while((*cl)[j]==0) j++;
     
     }
     for(unsigned int i=cnt_nz;i<2*cnt_nz;i++){
          (*cl)[i]=(*rw)[i-cnt_nz]; (*rw)[i]=(*cl)[i-cnt_nz]; (*vl)[i]=(*vl)[i-cnt_nz];//double edges to correspond to adjacency matrix format;
          
     
     }
     */
     
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
