































void mode_2_alg_w(node** node_arr,unsigned int *nds_td, unsigned int nds_n, node* node_hd, unsigned int max_nds, unsigned int** rw, unsigned int** cl, double **vl, unsigned int *ln, unsigned char out_fl,mode_3_param* mode3_inp,unsigned int shift,int fd, unsigned int* max_offst_elm,unsigned int m,unsigned int n){//nds_n0=0 as terminat, set max size at the beggining
     node *prev_node, *curr_node;
     prev_node=node_hd;
     unsigned int *ui_ptr; unsigned int ui_val;
     unsigned int mrg_cnt=0;
     struct nd_data* mrg=(struct nd_data*) malloc(1*sizeof(struct nd_data)); mrg->cap=0;
     
     (void)rw; (void)cl; (void)vl; (void)mode3_inp; (void)ln;
     
     for(unsigned int i=0;i<nds_n;i++){
          unsigned int n0=nds_td[i],n1;
          unsigned int ui_curr=0;
          node* nd0=node_arr[nds_td[i]-1];
          write(fd,&(nd0->n),sizeof(unsigned int));
          ui_curr=nd0->num+shift;
          write(fd,&ui_curr,sizeof(unsigned int));
          if(shift==0){
               write(fd,nd0->nums,(nd0->n)*sizeof(unsigned int));
          }
          else{
               for(unsigned int j=0;j<nd0->n;j++){
                    ui_curr=nd0->nums[j]+shift;
                    write(fd,&ui_curr,sizeof(unsigned int));
               }
          }
          if(max_offst_elm!=NULL && (*max_offst_elm)<(nd0->n+2)) (*max_offst_elm)=(nd0->n+2);
          
          for(unsigned int j=0;j<nd0->n;j++){
               n1=nd0->nums[j];
               node* nd1=node_arr[n1-1];
               unsigned n0_n1_cap=min(max(nd0->n+nd1->n,2*mrg->cap),max_nds-node_hd->num);
               if(mrg->cap<min(nd0->n+nd1->n,max_nds-node_hd->num)){
                    if(mrg->cap>0){
                         free(mrg->nums);
                    }
                    mrg->cap=n0_n1_cap;
                    mrg->nums=(unsigned int*) malloc(mrg->cap*sizeof(unsigned int));
               }
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
                         mrg->nums[mrg_cnt]=nd1->nums[w];
                         mrg_cnt++;
                         q++; w++;
                    }
                    else if(nd0->nums[q]<nd1->nums[w]){
                         mrg->nums[mrg_cnt]=nd0->nums[q];
                         mrg_cnt++;
                         q++;
                    }
                    else{
                         mrg->nums[mrg_cnt]=nd1->nums[w];
                         mrg_cnt++;
                         w++;
                    }
                    
               }
               while(q<nd0->n){
                    if(nd0->nums[q]==n1){
                         q++;
                    
                    }
                    else{
                         mrg->nums[mrg_cnt]=nd0->nums[q];
                         mrg_cnt++;
                         q++;
                    
                    }
               
               }
               while(w<nd1->n){
                    if(nd1->nums[w]==n0){
                         w++;
                    
                    }
                    else{
                         mrg->nums[mrg_cnt]=nd1->nums[w];
                         mrg_cnt++;
                         w++;
                    
                    }
               
               }
               ui_val=nd1->cap; nd1->cap=mrg->cap; mrg->cap=ui_val;
               ui_ptr=nd1->nums; nd1->nums=mrg->nums; mrg->nums=ui_ptr;
               nd1->n=mrg_cnt;
               
          
          }
          
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
               curr_node=prev_node->next;
               i++;
          
          }
          else{
               prev_node=curr_node;
               curr_node=curr_node->next;
          
          
          }
     
     }
     if(mrg->cap>0){
          free(mrg->nums);
     }
     free(mrg);
     
     
     
     
     unsigned int ui_curr=0;
     write(fd,&(ui_curr),sizeof(unsigned int));//write zero to signal thansition to next mode (mode3 or output);
     if(out_fl==1){
          ui_curr=0;//len;
          for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
               ui_curr+=curr_node->n;
          
          }
          write(fd,&ui_curr,sizeof(unsigned int));//write output length;
          for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
               for(unsigned int i=0;i<curr_node->n;i++){
                    ui_curr=curr_node->num+shift;
                    write(fd,&(ui_curr),sizeof(unsigned int));
                    ui_curr=curr_node->nums[i]+shift;
                    write(fd,&(ui_curr),sizeof(unsigned int));
                    
               }
          
          }
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
          
          }
     
     }
     else{//TODO seprate files for m,n,mx_m_sz; separate max_szes (probably grow) per file; do offsets so that read can be done in parallel (prerecord, calculate, overwrite), write not in parrallel;
          unsigned int buff_sz=64; char ch_buff[buff_sz];
          snprintf(ch_buff,buff_sz,"more_fold_%ux%u_%u",m,n,mode3_inp->max_m_sz);
          int fd0=open(ch_buff, O_RDWR | O_CREAT | O_TRUNC,0644);
          (*(mode3_inp->n_th))=((mode3_inp->tgt_n1)%(mode3_inp->max_m_sz)==0)?(mode3_inp->tgt_n1)/(mode3_inp->max_m_sz):((mode3_inp->tgt_n1)/(mode3_inp->max_m_sz)+1);
          //unsigned int ***rw=mode3_inp->rw; unsigned int ***cl=mode3_inp->cl; double ***vl=mode3_inp->vl;
          //(*rw)=(unsigned int**) malloc((*(mode3_inp->n_th))*sizeof(unsigned int*));
          //(*cl)=(unsigned int**) malloc((*(mode3_inp->n_th))*sizeof(unsigned int*));
          //(*vl)=(double**) malloc((*(mode3_inp->n_th))*sizeof(double*));
          //(*(mode3_inp->ln))=(unsigned int*) malloc((*(mode3_inp->n_th))*sizeof(unsigned int));
          off_t* offsets=(off_t*) malloc(((*(mode3_inp->n_th))+1)*sizeof(off_t));
          //off_t curr_fl_pos=lseek(fd, 0, SEEK_CUR);//we write at the begining of the file, so current position is 0;
          write(fd0,offsets,((*(mode3_inp->n_th))+1)*sizeof(off_t));
          offsets[0]=lseek(fd0, 0, SEEK_CUR);//TODO check if it is consitent number in bytes;
          max_nds=max_nds-node_hd->num+1;
          out_fl=1;
          //unsigned int num_thr=*(mode3_inp->n_th);
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
                    //node_arr0[curr_node->num-node_hd->num]->vals=(double*) malloc(curr_node->cap*sizeof(double));
                    for(unsigned int j=0;j<curr_node->n;j++){
                         node_arr0[curr_node->num-node_hd->num]->nums[j]=curr_node->nums[j]-node_hd->num+1;
                         //node_arr0[curr_node->num-node_hd->num]->vals[j]=curr_node->vals[j];
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
               
               //mode_2_alg(node_arr0,nds_td00,nds_n00,node_hd0,max_nds, &((*rw)[i]), &((*cl)[i]), &((*vl)[i]),&((*(mode3_inp->ln))[i]),out_fl,NULL,node_hd->num-1,NULL);
               mode_2_alg_w(node_arr0,nds_td00,nds_n00,node_hd0,max_nds,NULL,NULL,NULL,NULL,out_fl,NULL,node_hd->num-1,fd0,NULL,m,n);
               offsets[i+1]=lseek(fd0, 0, SEEK_CUR);
               free(node_arr0);
               free(nds_td00);
               free(node_mem);
               
               
               
          
          }
          
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               //free(prev_node->vals);
          
          }
          lseek(fd0, 0, SEEK_SET);
          write(fd0,offsets,((*(mode3_inp->n_th))+1)*sizeof(off_t));
          
          free(offsets);
          close(fd0);
     
     
     
     }
     
     
}
