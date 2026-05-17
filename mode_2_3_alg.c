






























#include <time.h>
void mode_2_alg(node** node_arr,double* edge_arr,unsigned int *nds_td, unsigned int nds_n, node* node_hd, unsigned int max_nds, unsigned int** rw, unsigned int** cl, double **vl, unsigned int *ln, unsigned char out_fl,mode_3_param* mode3_inp){
     struct timespec curr_time; long long unsigned int tick,dt_time;
     printf("\nmode2");
     unsigned int *buff_mrg=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int mrg_cnt;
     node *prev_node, *curr_node;
     unsigned int *uip_tmp;
     clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     for(unsigned int i=0;i<nds_n;i++){
          double sum=0;
          unsigned int n0=nds_td[i],n1,n2;
          for(unsigned int j=0;j<node_arr[n0-1]->n;j++){
               n1=node_arr[n0-1]->nums[j];
               sum+=edge_arr[(n0-1)*max_nds+n1-1];
               
          }
          for(unsigned int j=0;j<node_arr[n0-1]->n;j++){
               n1=node_arr[n0-1]->nums[j]; mrg_cnt=0;
               for(unsigned int k=0;k<node_arr[n0-1]->n;k++){
                    if(j==k){
                         continue;
                    
                    }
                    else{
                         n2=node_arr[n0-1]->nums[k];
                         if(n1<n2 && edge_arr[(n1-1)*max_nds+n2-1]>0){
                              edge_arr[(n1-1)*max_nds+n2-1]+=(edge_arr[(n0-1)*max_nds+n1-1]*edge_arr[(n0-1)*max_nds+n2-1])/sum;
                         
                         }
                         else if(n1<n2){
                              edge_arr[(n1-1)*max_nds+n2-1]=(edge_arr[(n0-1)*max_nds+n1-1]*edge_arr[(n0-1)*max_nds+n2-1])/sum;
                         
                         }
                         else{
                              edge_arr[(n1-1)*max_nds+n2-1]=edge_arr[(n2-1)*max_nds+n1-1];
                         
                         }
                         
                    }
                    
               }
               unsigned int q=0,w=0;
               for(;q<node_arr[n0-1]->n && w<node_arr[n1-1]->n;){
                    if(node_arr[n0-1]->nums[q]==n1){
                         q++; continue;
                    }
                    if(node_arr[n1-1]->nums[w]==n0){
                         w++; continue;
                    }
                    if(node_arr[n0-1]->nums[q]<node_arr[n1-1]->nums[w]){
                         buff_mrg[mrg_cnt]=node_arr[n0-1]->nums[q];
                         q++; mrg_cnt++;
                    }
                    else if(node_arr[n0-1]->nums[q]==node_arr[n1-1]->nums[w]){
                         buff_mrg[mrg_cnt]=node_arr[n0-1]->nums[q];
                         q++; mrg_cnt++;
                         w++;
                    }
                    else{
                         buff_mrg[mrg_cnt]=node_arr[n1-1]->nums[w];
                         w++; mrg_cnt++;
                    }
                    
               
               }
               while(q<node_arr[n0-1]->n){
                    if(node_arr[n0-1]->nums[q]==n1){
                         q++; continue;
                    
                    }
                    else{
                         buff_mrg[mrg_cnt]=node_arr[n0-1]->nums[q];
                         q++; mrg_cnt++;
                    
                    }
               
               }
               while(w<node_arr[n1-1]->n){
                    if(node_arr[n1-1]->nums[w]==n0){
                         w++; continue;
                    
                    }
                    else{
                         buff_mrg[mrg_cnt]=node_arr[n1-1]->nums[w];
                         w++; mrg_cnt++;
                    
                    }
               
               }
               node_arr[n1-1]->n=mrg_cnt;
               uip_tmp=node_arr[n1-1]->nums; node_arr[n1-1]->nums=buff_mrg; buff_mrg=uip_tmp;
               
          
          }
          /*for(unsigned int j=0;j<node_arr[n0-1]->n;j++){//no neighbout of n0 should be accesed, sicne n0 is deleted;
               n1=node_arr[n0-1]->nums[j];
               edge_arr[(n0-1)*max_nds+n1-1]=-1;
               edge_arr[(n1-1)*max_nds+n0-1]=-1;
          
          }*/
          curr_node=node_arr[n0-1];
          if(curr_node==node_hd){
               node_hd=curr_node->next;
               
          }
          else{
               prev_node=node_hd;
               while(prev_node->next!=curr_node){
                    prev_node=prev_node->next;
               
               }
               prev_node->next=curr_node->next;
               
          
          }
          free(curr_node->nums);
          free(curr_node);
          
     }
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
                    (*cl)[cnt_out]=curr_node->num;
                    (*rw)[cnt_out]=curr_node->nums[i];
                    (*vl)[cnt_out]=edge_arr[(curr_node->num-1)*max_nds+curr_node->nums[i]-1];//->val;
                    cnt_out++;
                    
                    
               
               }
          
          }
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               free(prev_node);
          
          }
          
          
     
     }
     else{
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          (*(mode3_inp->n_th))=((mode3_inp->tgt_n1)%(mode3_inp->max_m_sz)==0)?(mode3_inp->tgt_n1)/(mode3_inp->max_m_sz):((mode3_inp->tgt_n1)/(mode3_inp->max_m_sz)+1);
          unsigned int ***rw=mode3_inp->rw; unsigned int ***cl=mode3_inp->cl; double ***vl=mode3_inp->vl;
          (*rw)=(unsigned int**) malloc((*(mode3_inp->n_th))*sizeof(unsigned int*));
          (*cl)=(unsigned int**) malloc((*(mode3_inp->n_th))*sizeof(unsigned int*));
          (*vl)=(double**) malloc((*(mode3_inp->n_th))*sizeof(double*));
          (*(mode3_inp->ln))=(unsigned int*) malloc((*(mode3_inp->n_th))*sizeof(unsigned int));
          for(unsigned int i=0;i<(*(mode3_inp->n_th));i++){
               node** node_arr0=(node**) malloc(max_nds*sizeof(node*));
               double* edge_arr0=(double*) malloc(max_nds*max_nds*sizeof(double));
               for(unsigned int j=0;j<max_nds;j++){
                    node_arr0[j]=NULL;
                    for(unsigned int k=0;k<max_nds;k++){
                         edge_arr0[j*max_nds+k]=-1;
                    
                    }
               
               }
               curr_node=node_hd;
               node_arr0[curr_node->num-1]=(node*) malloc(1*sizeof(node));
               node_arr0[curr_node->num-1]->num=curr_node->num;
               node_arr0[curr_node->num-1]->n=curr_node->n;
               node_arr0[curr_node->num-1]->nums=(unsigned int*) malloc((max_nds-1)*sizeof(unsigned int));
               for(unsigned int j=0;j<curr_node->n;j++){
                    node_arr0[curr_node->num-1]->nums[j]=curr_node->nums[j];
                    edge_arr0[(curr_node->num-1)*max_nds+curr_node->nums[j]-1]=edge_arr[(curr_node->num-1)*max_nds+curr_node->nums[j]-1];
               }
               prev_node=node_arr0[curr_node->num-1];
               curr_node=curr_node->next;
               for(;curr_node!=NULL;curr_node=curr_node->next){
                    node_arr0[curr_node->num-1]=(node*) malloc(1*sizeof(node));
                    node_arr0[curr_node->num-1]->num=curr_node->num;
                    node_arr0[curr_node->num-1]->n=curr_node->n;
                    node_arr0[curr_node->num-1]->nums=(unsigned int*) malloc((max_nds-1)*sizeof(unsigned int));
                    for(unsigned int j=0;j<curr_node->n;j++){
                         node_arr0[curr_node->num-1]->nums[j]=curr_node->nums[j];
                         edge_arr0[(curr_node->num-1)*max_nds+curr_node->nums[j]-1]=edge_arr[(curr_node->num-1)*max_nds+curr_node->nums[j]-1];
                    }
                    prev_node->next=node_arr0[curr_node->num-1];
                    prev_node=node_arr0[curr_node->num-1];
               }
               prev_node->next=NULL;
               unsigned int *nds_td00; unsigned int nds_n00;
               
               
               
               if(i<(*(mode3_inp->n_th))-1){
                    nds_n00=mode3_inp->tgt_n1-mode3_inp->max_m_sz;
                    nds_td00=(unsigned int*) malloc(nds_n00*sizeof(unsigned int));
                    for(unsigned int j=0;j<i*(mode3_inp->max_m_sz);j++){
                         nds_td00[j]=mode3_inp->nds_tgt[j];
                    
                    }
                    for(unsigned int j=(i+1)*(mode3_inp->max_m_sz);j<mode3_inp->tgt_n1;j++){
                         nds_td00[j-1*(mode3_inp->max_m_sz)]=mode3_inp->nds_tgt[j];
                    
                    }
                    
               }
               else{
                    nds_n00=i*(mode3_inp->max_m_sz);
                    nds_td00=(unsigned int*) malloc(nds_n00*sizeof(unsigned int));
                    for(unsigned int j=0;j<i*(mode3_inp->max_m_sz);j++){
                         nds_td00[j]=mode3_inp->nds_tgt[j];
                    
                    }
                    
               
               }
               node* node_hd0=node_arr0[node_hd->num-1];
               out_fl=1;
               clock_gettime(CLOCK_MONOTONIC,&curr_time);
               dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
               printf("\nmode3 set up time: %llu",dt_time);
               
               mode_2_alg(node_arr0,edge_arr0,nds_td00,nds_n00,node_hd0,max_nds, &((*rw)[i]), &((*cl)[i]), &((*vl)[i]),&((*(mode3_inp->ln))[i]),out_fl,NULL);
               free(node_arr0); free(edge_arr0);
               free(nds_td00);
               
               
               
          
          }
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               free(prev_node);
          
          }
     
     
     
     }
     free(buff_mrg);
     
     
}

