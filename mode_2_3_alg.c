































void mode_2_alg(node** node_arr,double* edge_arr,unsigned int *nds_td, unsigned int nds_n, node** node_hd_, unsigned int max_nds, unsigned int** rw, unsigned int** cl, double **vl, unsigned int *ln, unsigned char out_fl){
     node* node_hd=(*node_hd_);
     unsigned int *buff_mrg=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int mrg_cnt;
     node *prev_node, *curr_node;
     unsigned int *uip_tmp;
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
     //TODO mode 3 (recursive);
          
     
     
     
     *node_hd_=node_hd;
     
     }
     free(buff_mrg);
     
     
}
