






























#include <time.h>
void mode_2_alg(node** node_arr,unsigned int *nds_td, unsigned int nds_n, node* node_hd, unsigned int max_nds, unsigned int** rw, unsigned int** cl, double **vl, unsigned int *ln, unsigned char out_fl,mode_3_param* mode3_inp,unsigned int shift){
     struct timespec curr_time; long long unsigned int tick,dt_time;
     printf("\nmode2");
     unsigned int min_cap=64, max_cap=max_nds-1;
     struct nd_data* mrg=(struct nd_data*) malloc(1*sizeof(struct nd_data));
     mrg->cap=min_cap; mrg->nums=(unsigned int*) malloc(mrg->cap*sizeof(unsigned int));
     mrg->vals=(double*) malloc(mrg->cap*sizeof(double));
     node *prev_node, *curr_node;
     prev_node=node_hd;
     unsigned int *ui_ptr; unsigned int ui_val; double* d_ptr;
     unsigned int mrg_cnt=0;
     clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     for(unsigned int i=0;i<nds_n;i++){
          double sum=0;
          unsigned int n0=nds_td[i],n1;
          for(unsigned int j=0;j<node_arr[n0-1]->n;j++){
               sum+=node_arr[n0-1]->vals[j];
               
          }
          for(unsigned int j=0;j<node_arr[n0-1]->n;j++){
               n1=node_arr[n0-1]->nums[j];
               double cache_scl=node_arr[n0-1]->vals[j]/sum;
               mrg_cnt=0;
               unsigned int q=0,w=0;
               if(mrg->cap<min((node_arr[n1-1]->n + node_arr[n0-1]->n),max_nds-1)){
                    mrg->cap=min(max(node_arr[n1-1]->n+node_arr[n0-1]->n,2*mrg->cap),max_cap);
                    /*mrg->nums=(unsigned int*) realloc(mrg->nums,mrg->cap*sizeof(unsigned int));
                    mrg->vals=(double*) realloc(mrg->vals,mrg->cap*sizeof(double));*/
                    free(mrg->nums); free(mrg->vals);
                    mrg->nums=(unsigned int*) malloc(mrg->cap*sizeof(unsigned int));
                    mrg->vals=(double*) malloc(mrg->cap*sizeof(double));//no need to keep (copy during realloc) old data;
               
               }
               for(;q<node_arr[n0-1]->n && w<node_arr[n1-1]->n;){
                    if(node_arr[n0-1]->nums[q]==n1){
                         q++; continue;
                    }
                    if(node_arr[n1-1]->nums[w]==n0){
                         w++; continue;
                    }
                    if(node_arr[n0-1]->nums[q]==node_arr[n1-1]->nums[w]){
                         mrg->nums[mrg_cnt]=node_arr[n1-1]->nums[w];
                         mrg->vals[mrg_cnt]=node_arr[n1-1]->vals[w]+cache_scl*node_arr[n0-1]->vals[q];//(node_arr[n0-1]->vals[j]*node_arr[n0-1]->vals[q])/sum;
                         mrg_cnt++;
                         q++; w++;
                    }
                    else if(node_arr[n0-1]->nums[q]<node_arr[n1-1]->nums[w]){
                         mrg->nums[mrg_cnt]=node_arr[n0-1]->nums[q];
                         mrg->vals[mrg_cnt]=cache_scl*node_arr[n0-1]->vals[q];
                         mrg_cnt++;
                         q++;
                    }
                    else{
                         mrg->nums[mrg_cnt]=node_arr[n1-1]->nums[w];
                         mrg->vals[mrg_cnt]=node_arr[n1-1]->vals[w];
                         mrg_cnt++;
                         w++;
                    }
                    
               }
               while(q<node_arr[n0-1]->n){
                    if(node_arr[n0-1]->nums[q]==n1){
                         q++;
                    
                    }
                    else{
                         mrg->nums[mrg_cnt]=node_arr[n0-1]->nums[q];
                         mrg->vals[mrg_cnt]=cache_scl*node_arr[n0-1]->vals[q];
                         mrg_cnt++;
                         q++;
                    
                    }
               
               }
               while(w<node_arr[n1-1]->n){
                    if(node_arr[n1-1]->nums[w]==n0){
                         w++;
                    
                    }
                    else{
                         mrg->nums[mrg_cnt]=node_arr[n1-1]->nums[w];
                         mrg->vals[mrg_cnt]=node_arr[n1-1]->vals[w];
                         mrg_cnt++;
                         w++;
                    
                    }
               
               }
               ui_val=node_arr[n1-1]->cap; node_arr[n1-1]->cap=mrg->cap; mrg->cap=ui_val;
               ui_ptr=node_arr[n1-1]->nums; node_arr[n1-1]->nums=mrg->nums; mrg->nums=ui_ptr;
               d_ptr=node_arr[n1-1]->vals; node_arr[n1-1]->vals=mrg->vals; mrg->vals=d_ptr;
               node_arr[n1-1]->n=mrg_cnt;
               
          
          }
          curr_node=node_arr[n0-1];
          if(curr_node==node_hd){
               node_hd=curr_node->next;
               prev_node=node_hd;
               
          }
          else{
               while(prev_node->next!=curr_node){
                    prev_node=prev_node->next;

               }
               prev_node->next=curr_node->next;
               
          
          }
          free(curr_node->nums);
          free(curr_node->vals);
          //free(curr_node);
          
     }
     free(mrg->nums);
     free(mrg->vals);
     free(mrg);
     clock_gettime(CLOCK_MONOTONIC,&curr_time);
     dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
     if(out_fl==1){
          printf("\nnodes processed: %u",nds_n);
          printf("\ntotal time: %llu",dt_time);
     
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
                    (*cl)[cnt_out]=curr_node->num+shift;
                    (*rw)[cnt_out]=curr_node->nums[i]+shift;
                    (*vl)[cnt_out]=curr_node->vals[i];
                    cnt_out++;
                    
                    
               
               }
          
          }
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               free(prev_node->vals);
               //free(prev_node);
          
          }
          
          
     
     }
     else{
          
          (*(mode3_inp->n_th))=((mode3_inp->tgt_n1)%(mode3_inp->max_m_sz)==0)?(mode3_inp->tgt_n1)/(mode3_inp->max_m_sz):((mode3_inp->tgt_n1)/(mode3_inp->max_m_sz)+1);
          unsigned int ***rw=mode3_inp->rw; unsigned int ***cl=mode3_inp->cl; double ***vl=mode3_inp->vl;
          (*rw)=(unsigned int**) malloc((*(mode3_inp->n_th))*sizeof(unsigned int*));
          (*cl)=(unsigned int**) malloc((*(mode3_inp->n_th))*sizeof(unsigned int*));
          (*vl)=(double**) malloc((*(mode3_inp->n_th))*sizeof(double*));
          (*(mode3_inp->ln))=(unsigned int*) malloc((*(mode3_inp->n_th))*sizeof(unsigned int));
          max_nds=max_nds-node_hd->num+1;
          dt_time=0;
          for(unsigned int i=0;i<(*(mode3_inp->n_th));i++){
               clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
               node** node_arr0=(node**) malloc(max_nds*sizeof(node*));
               node* node_mem=(node*) malloc(max_nds*sizeof(node));
               node* nd_ptr=(node*) malloc(1*sizeof(node)); prev_node=nd_ptr;
               for(curr_node=node_hd; curr_node!=NULL; curr_node=curr_node->next){
                    node_arr0[curr_node->num-node_hd->num]=&(node_mem[curr_node->num-node_hd->num]);
                    node_arr0[curr_node->num-node_hd->num]->num=curr_node->num-node_hd->num+1;
                    node_arr0[curr_node->num-node_hd->num]->n=curr_node->n;
                    node_arr0[curr_node->num-node_hd->num]->cap=curr_node->cap;
                    node_arr0[curr_node->num-node_hd->num]->nums=(unsigned int*) malloc(curr_node->cap*sizeof(unsigned int));
                    node_arr0[curr_node->num-node_hd->num]->vals=(double*) malloc(curr_node->cap*sizeof(double));
                    for(unsigned int j=0;j<curr_node->n;j++){
                         node_arr0[curr_node->num-node_hd->num]->nums[j]=curr_node->nums[j]-node_hd->num+1;
                         node_arr0[curr_node->num-node_hd->num]->vals[j]=curr_node->vals[j];
                    }
                    prev_node->next=node_arr0[curr_node->num-node_hd->num];
                    prev_node=node_arr0[curr_node->num-node_hd->num];
               
               }
               prev_node->next=NULL;
               free(nd_ptr);
               unsigned int *nds_td00; unsigned int nds_n00;
               
               if(i<(*(mode3_inp->n_th))-1){
                    nds_n00=mode3_inp->tgt_n1-mode3_inp->max_m_sz;
                    nds_td00=(unsigned int*) malloc(nds_n00*sizeof(unsigned int));
                    for(unsigned int j=0;j<i*(mode3_inp->max_m_sz);j++){
                         nds_td00[j]=mode3_inp->nds_tgt[j]-node_hd->num+1;
                    
                    }
                    for(unsigned int j=(i+1)*(mode3_inp->max_m_sz);j<mode3_inp->tgt_n1;j++){
                         nds_td00[j-1*(mode3_inp->max_m_sz)]=mode3_inp->nds_tgt[j]-node_hd->num+1;
                    
                    }
                    
               }
               else{
                    nds_n00=i*(mode3_inp->max_m_sz);
                    nds_td00=(unsigned int*) malloc(nds_n00*sizeof(unsigned int));
                    for(unsigned int j=0;j<i*(mode3_inp->max_m_sz);j++){
                         nds_td00[j]=mode3_inp->nds_tgt[j]-node_hd->num+1;
                    
                    }
                    
               
               }
               node* node_hd0=node_arr0[0];
               out_fl=1;
               clock_gettime(CLOCK_MONOTONIC,&curr_time);
               dt_time+=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
               
               mode_2_alg(node_arr0,nds_td00,nds_n00,node_hd0,max_nds, &((*rw)[i]), &((*cl)[i]), &((*vl)[i]),&((*(mode3_inp->ln))[i]),out_fl,NULL,node_hd->num-1);
               free(node_arr0);
               free(nds_td00);
               free(node_mem);
               
               
               
          
          }
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               free(prev_node->vals);
               //free(prev_node);
          
          }
          printf("\nmode3 set up time: %llu",dt_time);
     
     
     
     }
     
     
}

