






























#include <time.h>
void mode_2_alg(node** node_arr,unsigned int *nds_td, unsigned int nds_n, node* node_hd, unsigned int max_nds, unsigned int** rw, unsigned int** cl, double **vl, unsigned int *ln, unsigned char out_fl,mode_3_param* mode3_inp,unsigned int shift,FILE* fp){
     struct timespec curr_time; long long unsigned int tick,dt_time;
     //if(fp!=NULL) fprintf(fp,"\nmode2; nds_td j iter\n");
     //unsigned int min_cap=64, max_cap=max_nds-1;
     node *prev_node, *curr_node;
     prev_node=node_hd;
     unsigned int *ui_ptr; unsigned int ui_val; double* d_ptr;
     unsigned int mrg_cnt=0; int num_thr;
     struct nd_data* mrg;
     unsigned int mrg_sz=(/*fl_pl==1 && */ fl_mem_pl==1 && omp_in_parallel()==0)?max_nds-node_hd->num:1;//after deletion in mode1 there is at most max_nds-node_hd->num+1 nodes; hence each node (node can not be a neighbour of itself) can have at most max_nds-node_hd->num neighbourhoods;
     mrg=(struct nd_data*) malloc(mrg_sz*sizeof(struct nd_data));
     for(unsigned int i=0;i<mrg_sz;i++){
          mrg[i].cap=0;
     }
     for(unsigned int i=0;i<nds_n;i++){
          double sum=0;
          unsigned int n0=nds_td[i],n1;
          node* nd0=node_arr[nds_td[i]-1];
          //#pragma omp parallel for schedule(static) reduction(+:sum)//should be ralativelly small;
          for(unsigned int j=0;j<nd0->n;j++){
               sum+=nd0->vals[j];
          }
          num_thr=nd0->n;
          //clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          #if (fl_pl==1)
               #pragma omp parallel for if(omp_in_parallel()==0 && max_nds>2*128*128+128+1) schedule(static) num_threads(num_thr<=omp_get_max_threads()?num_thr:omp_get_max_threads()) private(n1,mrg_cnt,ui_val,ui_ptr,d_ptr)
          #endif
          for(unsigned int j=0;j<nd0->n;j++){
               n1=nd0->nums[j];
               node* nd1=node_arr[n1-1];
               struct nd_data* curr_mrg=(mrg_sz>1)?&(mrg[j]):mrg;
               unsigned n0_n1_cap=min(max(nd0->n+nd1->n,2*curr_mrg->cap),max_nds-node_hd->num);
               
               if(curr_mrg->cap<min(nd0->n+nd1->n,max_nds-node_hd->num)){
                    if(curr_mrg->cap>0){
                         free(curr_mrg->nums); free(curr_mrg->vals);
                    }
                    curr_mrg->cap=n0_n1_cap;
                    curr_mrg->nums=(unsigned int*) malloc(curr_mrg->cap*sizeof(unsigned int));
                    curr_mrg->vals=(double*) malloc(curr_mrg->cap*sizeof(double));
               }
               double cache_scl=nd0->vals[j]/sum;
               mrg_cnt=0;
               unsigned int q=0,w=0;
               for(;q<nd0->n && w<nd1->n;){
                    if(nd0->nums[q]==n1){
                         q++; continue;
                    }
                    if(nd1->nums[w]==n0){
                         w++; continue;
                    }
                    if(nd0->nums[q]==nd1->nums[w]){
                         curr_mrg->nums[mrg_cnt]=nd1->nums[w];
                         curr_mrg->vals[mrg_cnt]=nd1->vals[w]+cache_scl*nd0->vals[q];//(nd0->vals[j]*nd0->vals[q])/sum;
                         mrg_cnt++;
                         q++; w++;
                    }
                    else if(nd0->nums[q]<nd1->nums[w]){
                         curr_mrg->nums[mrg_cnt]=nd0->nums[q];
                         curr_mrg->vals[mrg_cnt]=cache_scl*nd0->vals[q];
                         mrg_cnt++;
                         q++;
                    }
                    else{
                         curr_mrg->nums[mrg_cnt]=nd1->nums[w];
                         curr_mrg->vals[mrg_cnt]=nd1->vals[w];
                         mrg_cnt++;
                         w++;
                    }
                    
               }
               while(q<nd0->n){
                    if(nd0->nums[q]==n1){
                         q++;
                    
                    }
                    else{
                         curr_mrg->nums[mrg_cnt]=nd0->nums[q];
                         curr_mrg->vals[mrg_cnt]=cache_scl*nd0->vals[q];
                         mrg_cnt++;
                         q++;
                    
                    }
               
               }
               while(w<nd1->n){
                    if(nd1->nums[w]==n0){
                         w++;
                    
                    }
                    else{
                         curr_mrg->nums[mrg_cnt]=nd1->nums[w];
                         curr_mrg->vals[mrg_cnt]=nd1->vals[w];
                         mrg_cnt++;
                         w++;
                    
                    }
               
               }
               ui_val=nd1->cap; nd1->cap=curr_mrg->cap; curr_mrg->cap=ui_val;
               ui_ptr=nd1->nums; nd1->nums=curr_mrg->nums; curr_mrg->nums=ui_ptr;
               d_ptr=nd1->vals; nd1->vals=curr_mrg->vals; curr_mrg->vals=d_ptr;
               nd1->n=mrg_cnt;
               
          
          }
          //clock_gettime(CLOCK_MONOTONIC,&curr_time);
          //dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
          //if(fp!=NULL) fprintf(fp,"%llu,",dt_time);
          
     }
     node nd_ptr;
     nd_ptr.next=node_hd;
     curr_node=node_hd; prev_node=&nd_ptr;
     for(unsigned int i=0;i<nds_n;){
          if(curr_node->num==nds_td[i]){
               prev_node->next=curr_node->next;
               if(curr_node==node_hd){
                    node_hd=curr_node->next;
               
               }
               free(curr_node->nums);
               free(curr_node->vals);
               curr_node=prev_node->next;
               i++;
          
          }
          else{
               prev_node=curr_node;
               curr_node=curr_node->next;
          
          
          }
     
     }
     /*free(mrg->nums);
     free(mrg->vals);
     free(mrg);*/
     for(unsigned int i=0;i<mrg_sz;i++){
          if(mrg[i].cap>0){
               free(mrg[i].nums); free(mrg[i].vals);
          }
     
     }
     free(mrg);
     if(out_fl==1){
          ////if(fp!=NULL) fprintf(fp,"\nnodes processed: %u",nds_n);
          //printf("\ntotal time: %llu",dt_time);
     
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
          out_fl=1;
          num_thr=*(mode3_inp->n_th);
          //fprintf(fp,"\n\nmode3; num_thr=,%d\n",num_thr);
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          #if fl_pl==1
               #pragma omp parallel for schedule(static) num_threads((num_thr<=omp_get_max_threads())?num_thr:omp_get_max_threads()) private(curr_node,prev_node) if((mode3_inp->tgt_n1>=64 && mode3_inp->max_m_sz>3) || mode3_inp->tgt_n1>64)
          #endif
          for(unsigned int i=0;i<(*(mode3_inp->n_th));i++){
               node** node_arr0=(node**) malloc(max_nds*sizeof(node*));
               node* node_mem=(node*) malloc(max_nds*sizeof(node));
               node nd_ptr; prev_node=&nd_ptr;
               for(curr_node=node_hd; curr_node!=NULL; curr_node=curr_node->next){
                    node_arr0[curr_node->num-node_hd->num]=&(node_mem[curr_node->num-node_hd->num]);
                    node_arr0[curr_node->num-node_hd->num]->num=curr_node->num-node_hd->num+1;
                    node_arr0[curr_node->num-node_hd->num]->n=curr_node->n;
                    node_arr0[curr_node->num-node_hd->num]->cap=curr_node->cap;//TODO see if memory capped, try different (node->n) capacities;
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
               //out_fl=1;
               //clock_gettime(CLOCK_MONOTONIC,&curr_time);
               //dt_time+=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
               
               mode_2_alg(node_arr0,nds_td00,nds_n00,node_hd0,max_nds, &((*rw)[i]), &((*cl)[i]), &((*vl)[i]),&((*(mode3_inp->ln))[i]),out_fl,NULL,node_hd->num-1,NULL);
               free(node_arr0);
               free(nds_td00);
               free(node_mem);
               
               
               
          
          }
          clock_gettime(CLOCK_MONOTONIC,&curr_time);
          dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
          fprintf(fp,"%llu,%d,%d",dt_time,(mode3_inp->max_m_sz),(*(mode3_inp->n_th)));
          
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               free(prev_node->vals);
               //free(prev_node);
          
          }
     
     
     
     }
     
     
}

