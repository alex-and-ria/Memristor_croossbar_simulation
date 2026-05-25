
































#include<stdio.h>


typedef struct node{
     unsigned int num;
     unsigned int n;
     unsigned int cap;
     unsigned int* nums;
     double* vals;
     struct node* next;
} node;

typedef struct{//TODO per node timing, time neighb, fix mode 3, cap;
     unsigned int ***rw; unsigned int*** cl; double***vl; unsigned int** ln; unsigned int* nds_tgt; unsigned int tgt_n1; unsigned int max_m_sz; unsigned int* n_th;
} mode_3_param;

struct nd_data{unsigned int cap; unsigned int* nums; double* vals;};

#define min(a,b) ((a<b)?(a):(b))
#define max(a,b) ((a>b)?(a):(b))

#include"mode_2_3_alg.c"

void mode_1_alg(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,double** vl, unsigned int len,unsigned int *ln, unsigned int *nds_td, unsigned int nds_n, double thr_koef, unsigned char out_fl, mode_3_param* mode3_inp){
     unsigned int max_nds=col[len-1];
     unsigned int min_cap=64, max_cap=max_nds-1;
     struct timespec curr_time; long long unsigned int tick,dt_time, nb_time=0, transf_time=0, frr_time=0;
     clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     node** node_arr=(node**) malloc(max_nds*sizeof(node*));
     unsigned int* nds_td0=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int* nds_td_rem=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int* nds_td2=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     struct nd_data* mrg=(struct nd_data*) malloc(1*sizeof(struct nd_data));
     mrg->cap=min_cap; mrg->nums=(unsigned int*) malloc(mrg->cap*sizeof(unsigned int));
     mrg->vals=(double*) malloc(mrg->cap*sizeof(double));
     double *sums=(double*)malloc(nds_n*sizeof(double));

     node* node_mem=(node*) malloc(max_nds*sizeof(node));//for better cache locality;
     node* node_hd=&(node_mem[0]);
     node_hd->num=col[0]; node_hd->cap=min_cap;
     node_hd->nums=(unsigned int*)malloc(node_hd->cap*sizeof(unsigned int));
     node_hd->vals=(double*)malloc(node_hd->cap*sizeof(double));
     node* curr_node=node_hd;
     node_arr[0]=curr_node;
     for(unsigned int i=col[0]+1;i<=col[len-1];i++){//assumption here is that node numbering is sequential without skipping the numbers;
          node_arr[i-1]=&(node_mem[i-1]);
          node_arr[i-1]->num=i;
          node_arr[i-1]->cap=min_cap;
          node_arr[i-1]->nums=(unsigned int*)malloc(node_arr[i-1]->cap*sizeof(unsigned int));
          node_arr[i-1]->vals=(double*)malloc(node_arr[i-1]->cap*sizeof(double));
          node_arr[i-2]->next=node_arr[i-1];
          
     }
     node_arr[col[len-1]-1]->next=NULL;
     ///////////////////////////////at this point nodes should be set up;
     
     
     
     
     //edge *curr_edge=NULL;
     unsigned int curr_col=col[0];
     unsigned int i0=0;
     for(unsigned int i=1;i<len;i++){
          if(curr_col!=col[i]){
               node_arr[col[i-1]-1]->n=i-i0;
               if(node_arr[col[i-1]-1]->n>node_arr[col[i-1]-1]->cap){
                    node_arr[col[i-1]-1]->cap=min(max(node_arr[col[i-1]-1]->n,2*node_arr[col[i-1]-1]->cap),max_cap);
                    node_arr[col[i-1]-1]->nums=(unsigned int*)realloc(node_arr[col[i-1]-1]->nums,node_arr[col[i-1]-1]->cap*sizeof(unsigned int));
                    node_arr[col[i-1]-1]->vals=(double*)realloc(node_arr[col[i-1]-1]->vals,node_arr[col[i-1]-1]->cap*sizeof(double));
               }
               for(unsigned int j=i0;j<i;j++){
                    node_arr[col[i-1]-1]->nums[j-i0]=row[j];
                    node_arr[col[i-1]-1]->vals[j-i0]=val[j];
               
               }
               i0=i;
               curr_col=col[i];
          }
          
     
     }
     node_arr[max_nds-1]->n=len-i0;
     if(node_arr[max_nds-1]->n>node_arr[max_nds-1]->cap){
          node_arr[max_nds-1]->cap=min(max(node_arr[max_nds-1]->n,2*node_arr[max_nds-1]->cap),max_cap);
          node_arr[max_nds-1]->nums=(unsigned int*)realloc(node_arr[max_nds-1]->nums,node_arr[max_nds-1]->cap*sizeof(unsigned int));
          node_arr[max_nds-1]->vals=(double*)realloc(node_arr[max_nds-1]->vals,node_arr[max_nds-1]->cap*sizeof(double));
     }
     for(unsigned int j=i0;j<len;j++){
          node_arr[max_nds-1]->nums[j-i0]=row[j];
          node_arr[max_nds-1]->vals[j-i0]=val[j];
     
     }
     /////////////////////////////////////at this point edges should be set up too;
     clock_gettime(CLOCK_MONOTONIC,&curr_time);
     dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
     printf("\nmode1:\nset up time: %llu",dt_time);
     
     
     
     
     
     
     
     
     
     
     
     unsigned int nds_n0=0,nds_n_rem=0,nds_n2=nds_n;
     for(unsigned int i=0;i<nds_n;i++){
          nds_td2[i]=nds_td[i];
     
     }
     
     double curr_thr_koef=0;
     while(1/*dbg_cnt<dbg_max*/){
     clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     for(unsigned int i=0;i<max_nds;i++){
          nds_td0[i]=0; nds_td_rem[i]=0;
     }
     for(unsigned int i=0;i<nds_n2;i++){
          nds_td0[nds_td2[i]-1]=nds_td2[i]; nds_td_rem[nds_td2[i]-1]=nds_td2[i];
     }
     for(unsigned int i=0;i<max_nds;){
          while(i<max_nds && nds_td0[i]==0){
               i++;
          
          }
          if(i<max_nds){
               curr_node=node_arr[i];//nds_td0[i] should be equal to i+1 if it is not zero; curr_node->num should also be equal to i+1;
               for(unsigned int j=0;j<curr_node->n;j++){
                    if(nds_td0[curr_node->nums[j]-1]!=0 && nds_td0[i] < curr_node->nums[j]) nds_td0[curr_node->nums[j]-1]=0;
                    for(unsigned int k=0;k<node_arr[curr_node->nums[j]-1]->n;k++){
                         if(nds_td0[node_arr[curr_node->nums[j]-1]->nums[k]-1]!=0 && nds_td0[i] < node_arr[curr_node->nums[j]-1]->nums[k]){
                              nds_td0[node_arr[curr_node->nums[j]-1]->nums[k]-1]=0;
                         
                         }
                    
                    }
               
               }
               i++;
          }
          
     }
     for(unsigned int i=0;i<max_nds;i++){
          if(nds_td0[i]!=0){
               nds_td_rem[nds_td0[i]-1]=0;
               nds_td0[nds_n0]=nds_td0[i];
               nds_n0++;
          
          }
          if(nds_td_rem[i]!=0){
               nds_td_rem[nds_n_rem]=nds_td_rem[i];
               nds_n_rem++;
          
          }
     
     }
     if(nds_n0>=1) curr_thr_koef=(nds_n0+0.)/(nds_n0+nds_n_rem);
     if(nds_n0<=1 || curr_thr_koef<thr_koef){
          break;
     
     }
     /////////////////////////////////////////at this point nds_td0 (nodes to delete) and nds_td_rem (nodes left, delete on next iteration) should set up;
     clock_gettime(CLOCK_MONOTONIC,&curr_time);
     nb_time+=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
     //printf("\nneighbour search time: %llu",dt_time);
     
     
     clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     unsigned int n1;
     unsigned int bf_mg_cnt=0;
     unsigned int *ui_ptr; unsigned int ui_val; double* d_ptr;
    
     for(unsigned int i=0;i<nds_n0;i++){
          sums[i]=0;
          n1=nds_td0[i];
          for(unsigned int j=0;j<node_arr[n1-1]->n;j++){
               sums[i]+=node_arr[n1-1]->vals[j];
               
          
          }
     }
     
     for(unsigned int i=0;i<nds_n0;i++){
          for(unsigned int j=0;j<node_arr[nds_td0[i]-1]->n;j++){
               n1=node_arr[nds_td0[i]-1]->nums[j];
               double cache_scl=node_arr[nds_td0[i]-1]->vals[j]/sums[i];
               /*unsigned int idx0=((2*node_arr[nds_td0[i]-1]->n-j-1)*j)/2;//sum from n-1 to n-j (number of elements that should be calculated before calculating values for j-th index);
               for(unsigned int k=j+1;k<node_arr[nds_td0[i]-1]->n;k++){
                    curr_vals[i*nds_n0+idx0+k-j-1]=(node_arr[nds_td0[i]-1]->vals[j]*node_arr[nds_td0[i]-1]->vals[k])/sum[i];
                    int ll=(node_arr[n1-1]->n+1)/2; unsigned int diff=(ll+1)/2; unsigned char fl_fnd=0;
                    while(diff>0){
                         if(node_arr[n1-1]->nums[ll]>n2){
                              ll-=diff;
                              if(ll<0) ll=0;
                              diff=(diff+1)/2;
                              
                         
                         }
                         else if(node_arr[n1-1]->nums[ll]<n2){
                              ll+=diff;
                              if(ll>node_arr[nds_td0[i]-1]->n-1) ll=node_arr[nds_td0[i]-1]->n-1;
                              diff=(diff+1)/2;
                         
                         }
                         else{
                              fl_fnd=1;
                              break;
                              
                         
                         }
                    
                    }
                    
                    if(fl_fnd==1){
                         curr_vals[i*nds_n0+idx0]+=node_arr[n1-1]->vals
                         fl_fnd=0;
                    
                    }
               }*/
               bf_mg_cnt=0;
               if(mrg->cap<min((node_arr[n1-1]->n + node_arr[nds_td0[i]-1]->n),max_nds-1)){
                    mrg->cap=min(max(node_arr[n1-1]->n+node_arr[nds_td0[i]-1]->n,2*mrg->cap),max_cap);
                    /*mrg->nums=(unsigned int*) realloc(mrg->nums,mrg->cap*sizeof(unsigned int));
                    mrg->vals=(double*) realloc(mrg->vals,mrg->cap*sizeof(double));*/
                    free(mrg->nums); free(mrg->vals);
                    mrg->nums=(unsigned int*) malloc(mrg->cap*sizeof(unsigned int));
                    mrg->vals=(double*) malloc(mrg->cap*sizeof(double));//no need to keep (copy during realloc) old data;
                    
               
               }
               unsigned int ll=0,qq=0;
               for(; ll<node_arr[n1-1]->n && qq<node_arr[nds_td0[i]-1]->n;){
                    if(node_arr[nds_td0[i]-1]->nums[qq]==n1){//skip node from the list of neighbourds of the node to delete that has same index as the neighbout itself (form edges with all pivot's neighbours except itself);
                         qq++; continue;
                    
                    }
                    if(node_arr[n1-1]->nums[ll]==nds_td0[i]){//for node node_td0[i] and its neighbours nums=node_arr[node_td[i]-1]->nums, each of nodes in nums has node_td0[i] as their own neighbour (if a has edge with b then b has edge with a), so skip nds_td0[i] node in its neighbour list of neighbours; no other node from nds_td0 is in neigbour list because othewise this node whould be the common node (common neigbour) for two nodes from nds_td0;
                         ll++; continue;
                    
                    }
                    if(node_arr[n1-1]->nums[ll] == node_arr[nds_td0[i]-1]->nums[qq]){//edge exists, parralel conductances are added;
                         mrg->nums[bf_mg_cnt]=node_arr[n1-1]->nums[ll];
                         mrg->vals[bf_mg_cnt]=node_arr[n1-1]->vals[ll]+cache_scl*node_arr[nds_td0[i]-1]->vals[qq];//(node_arr[nds_td0[i]-1]->vals[j]*node_arr[nds_td0[i]-1]->vals[qq])/sums[i];
                         bf_mg_cnt++; ll++; qq++;
                    }
                    else if(node_arr[n1-1]->nums[ll] < node_arr[nds_td0[i]-1]->nums[qq]){//record old edge;
                         mrg->nums[bf_mg_cnt]=node_arr[n1-1]->nums[ll];
                         mrg->vals[bf_mg_cnt]=node_arr[n1-1]->vals[ll];
                         bf_mg_cnt++; ll++;
                    
                    }
                    else{//calculate and recocrd new edge;
                         mrg->nums[bf_mg_cnt]=node_arr[nds_td0[i]-1]->nums[qq];
                         mrg->vals[bf_mg_cnt]=cache_scl*node_arr[nds_td0[i]-1]->vals[qq];
                         bf_mg_cnt++; qq++;
                         
                    }
               
               }
               while(ll<node_arr[n1-1]->n){
                    if(node_arr[n1-1]->nums[ll]==nds_td0[i]){
                         ll++;
                    }
                    else{
                         mrg->nums[bf_mg_cnt]=node_arr[n1-1]->nums[ll];
                         mrg->vals[bf_mg_cnt]=node_arr[n1-1]->vals[ll];
                         bf_mg_cnt++; ll++;
                    }
               
               }
               while(qq<node_arr[nds_td0[i]-1]->n){
                    if(node_arr[nds_td0[i]-1]->nums[qq]==n1){
                         qq++; continue;
                    
                    }
                    else{
                         mrg->nums[bf_mg_cnt]=node_arr[nds_td0[i]-1]->nums[qq];
                         mrg->vals[bf_mg_cnt]=cache_scl*node_arr[nds_td0[i]-1]->vals[qq];
                         bf_mg_cnt++; qq++;
                    
                    }
               
               }
               ui_val=node_arr[n1-1]->cap; node_arr[n1-1]->cap=mrg->cap; mrg->cap=ui_val;
               ui_ptr=node_arr[n1-1]->nums; node_arr[n1-1]->nums=mrg->nums; mrg->nums=ui_ptr;
               d_ptr=node_arr[n1-1]->vals; node_arr[n1-1]->vals=mrg->vals; mrg->vals=d_ptr;
               node_arr[n1-1]->n=bf_mg_cnt;
          
          }
     
     }
     /////////////////////////////////////////////////////////////////////at this point star mesh transform should be completed for this iteration
     clock_gettime(CLOCK_MONOTONIC,&curr_time);
     transf_time+=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
     //printf("\nstar_to_mesh: %llu",dt_time);
     
     clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     node* prev_node=NULL;curr_node=node_hd;
     for(unsigned int i=0; i<nds_n0 && curr_node!=NULL;){
          if(curr_node->num==nds_td0[i]){
               free(curr_node->nums);
               free(curr_node->vals);
               node_arr[curr_node->num-1]=NULL;
               if(prev_node==NULL){
                    node_hd=curr_node->next;
                    //free(curr_node);
                    curr_node=node_hd;
               
               }
               else{
                    prev_node->next=curr_node->next;
                    //free(curr_node);
                    curr_node=prev_node->next;
               
               }
               i++;
          
          }
          else if(curr_node->num<nds_td0[i]){
               prev_node=curr_node;
               curr_node=curr_node->next;
          
          }
     }
     /////////////////////////////////at this point nums array are edited to match new graph connectivity;
     
     
     
     ui_ptr=nds_td2; nds_td2=nds_td_rem; nds_td_rem=ui_ptr;
     nds_n2=nds_n_rem; nds_n0=0; nds_n_rem=0;
     clock_gettime(CLOCK_MONOTONIC,&curr_time);
     frr_time+=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
     //printf("\nclean up time: %llu",dt_time);
     //dbg_cnt++;
     }
     
     if(out_fl==1){
          *ln=0;
          for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
               *ln+=curr_node->n;
          
          }
          (*cl)=(unsigned int*) malloc((*ln)*sizeof(unsigned int));
          (*rw)=(unsigned int*) malloc((*ln)*sizeof(unsigned int));
          (*vl)=(double*) malloc((*ln)*sizeof(double));
          unsigned int cnt_out=0;
          for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
               for(unsigned int i=0;i<curr_node->n;i++){
                    (*cl)[cnt_out]=curr_node->num;
                    (*rw)[cnt_out]=curr_node->nums[i];
                    (*vl)[cnt_out]=curr_node->vals[i];
                    cnt_out++;
                    
               
               }
          
          }
          free(nds_td0);
          free(nds_td_rem);
          free(sums);
          node* prev_node=NULL;
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               free(prev_node->vals);
               //free(prev_node);
               
            
          
          }
          free(node_arr);
          
          
          
     }
     else{
          out_fl--;
          mode_2_alg(node_arr,nds_td2,nds_n2, node_hd,max_nds,rw,cl,vl,ln,out_fl,mode3_inp,0);
          //mode_2_alg(node** node_arr,unsigned int *nds_td, unsigned int nds_n, node** node_hd_, unsigned int max_nds, unsigned int** rw, unsigned int** cl, double **vl, unsigned int *ln, unsigned char out_fl);
          free(nds_td0);
          free(nds_td_rem);
          free(sums);
          free(node_arr);
          
          
     
     
     }
     
     free(node_mem);
     free(nds_td2);
     free(mrg->nums);
     free(mrg->vals);
     free(mrg);
     printf("\nnodes processed: %u",nds_n-nds_n2);
     printf("\nneighbour search time: %llu",nb_time);
     printf("\nstar_to_mesh: %llu",transf_time);
     printf("\nclean up time: %llu",frr_time);
     

}

