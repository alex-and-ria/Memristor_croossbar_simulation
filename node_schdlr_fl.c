
































#define flnm_str "%ux%u_%.5f.txt"

void node_analyzer_fl(unsigned int* row, unsigned int *col, unsigned int len, unsigned int* nds_td, unsigned int nds_n, unsigned int* nds_td1,unsigned int* nds_n1, double* nb_koef, unsigned int* n_indp){//adj_m starts count from node 1 (Matlab indexing); it is represented in sparse Matlab format via [row,col,val]; col are sorted; nds_td should count from one (Matlab idexing), and be sorted; it is a list of nodes to delete; check nodes and delete one of them from nds_td if two are neighbours;
     unsigned int cnt_zr=0;
     unsigned int j=0;
     for(unsigned int i=0;i<nds_n;i++){
          if(nds_td[i]==0) continue;
          while(j<len && col[j]<nds_td[i]) j++;
          
          for(unsigned int k=i+1; j<len && k<nds_n && col[j]==nds_td[i];){
               if(row[j]<nds_td[k]){
                    j++;
               }
               else if(row[j]>nds_td[k]){
                    k++;
               }
               else{
                    nds_td1[cnt_zr]=nds_td[k];
                    nds_td[k]=0;
                    j++; k++;
                    cnt_zr++;
               }
          }
          
     }
     (*nb_koef)=(cnt_zr+0.)/(nds_n+0.);
     *n_indp=(nds_n-cnt_zr);
     *nds_n1=cnt_zr;
}

void schr_rec(unsigned int *rw,unsigned int *cl, unsigned int len, unsigned int *nds_td, unsigned int nds_n,double th_nb_koef,unsigned int m, unsigned int n, unsigned char** fl_nm, unsigned int*** tot_mem_buff, unsigned int* mem_buff_cnt){
     (*tot_mem_buff)=(unsigned int**) malloc(16*sizeof(unsigned int*)); *mem_buff_cnt=0;
     (*fl_nm)=(unsigned char*) malloc(64*sizeof(unsigned char));// (*tot_mem_buff)[(*mem_buff_cnt)]=(*fl_nm); (*mem_buff_cnt)++;
     sprintf((char*)(*fl_nm),flnm_str,m,n,th_nb_koef);
     FILE* fp=fopen((char*)(*fl_nm),"wb");
     unsigned int max_edges=(cl[len-1]*(cl[len-1]-1));//maximun node number is in the end position of cl array; in fully connected network with n nodes there are n*(n-1)/2 edges; in adjacency matrix each edge is recorded twice (col_i row_j and row_i col_j);
     unsigned int len0=len, nds_n0=nds_n;
     fprintf(fp,"mode1: %u\n",nds_n0);
     unsigned int* rw0=(unsigned int*) malloc(max_edges*sizeof(unsigned int));(*tot_mem_buff)[(*mem_buff_cnt)]=rw0; (*mem_buff_cnt)++;
     unsigned int* cl0=(unsigned int*) malloc(max_edges*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=cl0; (*mem_buff_cnt)++;
     unsigned int* nds_td0=(unsigned int*) malloc(nds_n0*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=nds_td0; unsigned int n0=(*mem_buff_cnt); (*mem_buff_cnt)++;
     for(unsigned int i=0;i<len0;i++){
          rw0[i]=rw[i]; cl0[i]=cl[i];
     
     }
     for(unsigned int i=0;i<nds_n0;i++){
          nds_td0[i]=nds_td[i];
          fprintf(fp,"%u ",nds_td0[i]);
     
     }
     //allocating maximun necesssary memory whenever possible, realloc later if needed or to save some memory;
     unsigned int* nds_td1=(unsigned int*)malloc(nds_n0*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=nds_td1; (*mem_buff_cnt)++; //nodes that are not deleted now, but should be deleted on next iterations;
     double nb_koef; unsigned int n_indp; unsigned int nds_n1;
     node_analyzer_fl(rw0,cl0,len0,nds_td0, nds_n0, nds_td1,&nds_n1, &nb_koef, &n_indp);
     
     unsigned int *mesh_bdrs=(unsigned int*) malloc((nds_n0+2)*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=mesh_bdrs;  unsigned int n1=(*mem_buff_cnt); (*mem_buff_cnt)++;//mesh boundaries for meshes arrays stored contiguously;
     unsigned int* idx_r_td=(unsigned int*) malloc(2*nds_n0*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=idx_r_td;  unsigned int n2=(*mem_buff_cnt); (*mem_buff_cnt)++;//each node of nds_n0 has opening and closing index; for the nodes to delete in current iterarion;
     unsigned int* idx_r_nb=(unsigned int*) malloc(2*len0*sizeof(unsigned int)); unsigned int idx_rnb_len;// (*tot_mem_buff)[(*mem_buff_cnt)]=idx_r_nb; (*mem_buff_cnt)++;//for the nodes to delete in folloving iterations (so keep them in current iteration);
     unsigned int* cnt_tmps=(unsigned int*) malloc((nds_n0+1)*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=cnt_tmps;  unsigned int n3=(*mem_buff_cnt); (*mem_buff_cnt)++;//array for indexing for meshes;
     //unsigned int* nb_idx_arr=(unsigned int*) malloc(max_edges*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=nb_idx_arr; (*mem_buff_cnt)++;
     unsigned int* rw_buff=(unsigned int*)malloc(max_edges*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=rw_buff; (*mem_buff_cnt)++;
     unsigned int* cl_buff=(unsigned int*)malloc(max_edges*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=cl_buff; (*mem_buff_cnt)++;
     //unsigned int* rw_cl_idx=(unsigned int*)malloc(3*max_edges*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=rw_cl_idx; (*mem_buff_cnt)++;//idexes to save to file to form rw_buff1 and cl_buff1;
     unsigned int* rw_buff1=(unsigned int*)malloc(max_edges*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=rw_buff1; (*mem_buff_cnt)++;
     unsigned int* cl_buff1=(unsigned int*)malloc(max_edges*sizeof(unsigned int)); (*tot_mem_buff)[(*mem_buff_cnt)]=cl_buff1; (*mem_buff_cnt)++;
     unsigned int* idx_mrg=(unsigned int*)malloc(2*max_edges*sizeof(unsigned int));
     unsigned int* idx_el=(unsigned int*)malloc(max_edges*sizeof(unsigned int));
     unsigned int* idx_remd=(unsigned int*)malloc(2*max_edges*sizeof(unsigned int));
     
     
     
     unsigned int n_iter_dbg=1, n_iter=0;
     while (nb_koef<th_nb_koef && n_indp>1 && n_iter<n_iter_dbg){//continue mode 1 until threshold or until number of indipendent nodes more then 1;
          n_iter++;
          unsigned int tmp_cnt=0;
          for(unsigned int i=0;i<nds_n0;i++){//grouping nds_t0 to avoid '0';
               if(nds_td0[i]!=0){
                    nds_td0[tmp_cnt]=nds_td0[i];
                    tmp_cnt++;
               
               }
          
          }
          nds_n0=tmp_cnt;
          nds_td0=(unsigned int*) realloc(nds_td0,nds_n0*sizeof(unsigned int)); (*tot_mem_buff)[n0]=nds_td0;
          mesh_bdrs=(unsigned int*) realloc(mesh_bdrs,(nds_n0+2)*sizeof(unsigned int)); (*tot_mem_buff)[n1]=mesh_bdrs;
          idx_r_td=(unsigned int*) realloc(idx_r_td,2*nds_n0*sizeof(unsigned int)); (*tot_mem_buff)[n2]=idx_r_td;
          cnt_tmps=(unsigned int*) realloc(cnt_tmps,(nds_n0+1)*sizeof(unsigned int)); (*tot_mem_buff)[n3]=cnt_tmps;
          //fprintf nds_n0, nds_td0 //can be ommited if resulting row and col is recorded, then recorded and correspondingly used for transition;
          nds_sprt(cl0, &len0, &idx_r_td,&idx_r_nb,&idx_rnb_len, nds_td0, &nds_n0);
          mesh_bdrs[0]=0; unsigned int cnt_remd=0;
          for(unsigned int i=0;i<nds_n0;i++){
               mesh_bdrs[i+1]=mesh_bdrs[i]+((idx_r_td[2*i+1]-idx_r_td[2*i])*(idx_r_td[2*i+1]-idx_r_td[2*i]-1)+0.);//new mesh boundary strats where ends previous one plus mesh size (number of edges); each edge is recorded twice to comply with adjacency matrix requirenments;
               
          }
          
          for(unsigned int i=0; i<idx_rnb_len;i++){//marking with zeroes values in row that stay for nodes that are not deleted in this iteration;
               for(unsigned int j=0,kk=idx_r_nb[2*i]; j<nds_n0 && kk<idx_r_nb[2*i+1];){
                    while(kk<idx_r_nb[2*i+1] && rw0[kk]<nds_td0[j]){
                         cl_buff[mesh_bdrs[nds_n0]+cnt_remd]=cl0[kk];
                         rw_buff[mesh_bdrs[nds_n0]+cnt_remd]=rw0[kk];
                         idx_remd[cnt_remd]=kk;
                         cnt_remd++;
                         kk++;
                    
                    }
                    if(kk<idx_r_nb[2*i+1] && rw0[kk]==nds_td0[j]){
                         //rw0[kk]=0;
                         kk++;
                         j++;
                    }
                    else if(kk<idx_r_nb[2*i+1] && rw0[kk]>nds_td0[j]){
                         j++;
                    }
                    if(j==nds_n0){
                         for(;kk<idx_r_nb[2*i+1];){
                              cl_buff[mesh_bdrs[nds_n0]+cnt_remd]=cl0[kk];
                              rw_buff[mesh_bdrs[nds_n0]+cnt_remd]=rw0[kk];
                              idx_remd[cnt_remd]=kk;
                              cnt_remd++;
                              kk++;

                         }

                    }
                    
                    
               
               }
          
          }
          mesh_bdrs[nds_n0+1]=mesh_bdrs[nds_n0]+cnt_remd;
          
          fprintf(fp,"\ncl_buff, rw_buff (remdr): %u to %u\n",mesh_bdrs[nds_n0],mesh_bdrs[nds_n0]+cnt_remd);
          for(unsigned int i=0;i<cnt_remd;i++){
               fprintf(fp,"%u ",idx_remd[i]);
          
          }
          fprintf(fp,"\ncurr_nds_n (idx_r_td,mesh_bdrs): %u\n",nds_n0);
          for(unsigned int i=0;i<nds_n0;i++){
               fprintf(fp,"%u %u %u ",idx_r_td[2*i],idx_r_td[2*i+1], mesh_bdrs[i]);
          
          }
          #ifdef _OPENMP
               unsigned int num_threads_val=(nds_n0<(unsigned int) omp_get_max_active_levels())?nds_n0:((unsigned int)omp_get_max_active_levels());
          #endif
          #pragma omp parallel for num_threads(num_threads_val)
          for(unsigned int i=0;i<nds_n0;i++){
               for(unsigned int j=idx_r_td[2*i];j<idx_r_td[2*i+1]-1;j++){
                    for(unsigned int k=j+1;k<idx_r_td[2*i+1];k++){
                         cnt_tmps[i]=mesh_bdrs[i]+(idx_r_td[2*i+1]-idx_r_td[2*i]-1)*(j-idx_r_td[2*i])+(k-idx_r_td[2*i]-1);
                         cl_buff[cnt_tmps[i]]=rw0[j];//column number smaller than row number;
                         rw_buff[cnt_tmps[i]]=rw0[k];
                         cnt_tmps[i]=mesh_bdrs[i]+(idx_r_td[2*i+1]-idx_r_td[2*i]-1)*(k-idx_r_td[2*i])+(j-idx_r_td[2*i]);
                         cl_buff[cnt_tmps[i]]=rw0[k];//column number bigger than row number;
                         rw_buff[cnt_tmps[i]]=rw0[j];
                         
                    
                    }
               
               }
               
          }
          /////////////////////////////////////
          unsigned int curr_i, idx_m0=0, idx_res=0,idx_mrg_cnt=0; cnt_tmps[0]=0; 
          struct nds_lst{unsigned int i; struct nds_lst* next;}; struct nds_lst* nds_lst_head=(struct nds_lst*) malloc(1*sizeof(struct nds_lst)); struct nds_lst* nl_prev=nds_lst_head, *nl_curr, *nl_m0; nds_lst_head->i=0; nds_lst_head->next=NULL;
          for(unsigned int i=1;i<nds_n0+1;i++){
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
               for(nl_curr=nds_lst_head->next;nl_curr!=NULL;){
                    curr_i=nl_curr->i;
                    if(cl_buff[cnt_tmps[idx_m0]] == cl_buff[cnt_tmps[curr_i]] 
                    && rw_buff[cnt_tmps[idx_m0]] == rw_buff[cnt_tmps[curr_i]]){
                         //vl_buff[cnt_tmps[idx_m0]]+=vl_buff[cnt_tmps[curr_i]];
                         //frintf fl_merge if(nds_n0_m0) nb_idx_arr[cnt_tmps[idx_m0]-mesh_bdrs[nds_n0]] cnt_tmps[curr_i];
                         idx_mrg[2*idx_mrg_cnt]=cnt_tmps[idx_m0]; idx_mrg[2*idx_mrg_cnt+1]=cnt_tmps[curr_i]; idx_mrg_cnt++;
                     
                         cnt_tmps[curr_i]++;
                         if(cnt_tmps[curr_i]==mesh_bdrs[curr_i+1]){
                              nl_prev->next=nl_curr->next; 
                              free(nl_curr); 
                              nl_curr=nl_prev->next;
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
               idx_el[idx_res]=cnt_tmps[idx_m0];
               //vl_buff1[idx_res]=vl_buff[cnt_tmps[idx_m0]];
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
               idx_el[idx_res]=i;
               //vl_buff1[idx_res]=vl_buff[i];
               idx_res++;
               
          }
          free(nds_lst_head);
          //fprinf idx_res, rw_cl_idx;
          fprintf(fp,"\nmerges: %u\n",idx_mrg_cnt);
          for(unsigned int i=0;i<idx_mrg_cnt;i++){
               fprintf(fp,"%u %u",idx_mrg[2*i],idx_mrg[2*i+1]);
          
          }
          fprintf(fp,"\nindexes: %u\n",idx_res);
          for(unsigned int i=0;i<idx_res;i++){
               fprintf(fp,"%u ",idx_el[i]);
          
          }
          //at this point sorted array from deleted nodes should be merged themselves;
          ///////////////////////////////////////////////
          unsigned int* ui_ptr_tmp=cl0; cl0=cl_buff1; cl_buff1=ui_ptr_tmp;//free(ui_ptr_tmp);
          ui_ptr_tmp=rw0; rw0=rw_buff1; rw_buff1=ui_ptr_tmp;
          len0=idx_res;
          ui_ptr_tmp=nds_td0; nds_td0=nds_td1; nds_td1=ui_ptr_tmp;
          nds_n0=nds_n1;
          node_analyzer_fl(rw0,cl0,len0,nds_td0, nds_n0, nds_td1,&nds_n1, &nb_koef, &n_indp);
           
          
          
          //free(idx_r); free(idx_r_nb);
          //free(mesh_bdrs); free(cnt_tmps);
          //free(rw_buff); free(cl_buff); free(vl_buff);
          //free(rw_buff1); free(cl_buff1); free(vl_buff1);
     
          
          
     }
     
     /*fprintf(fp,"\nnds_n1: %u\n",nds_n0);
     for(unsigned int i=0;i<nds_n1;i++){
          fprintf(fp,"%u ",nds_td1[i]);
     
     }*/
     fclose(fp);
     
     
     free(idx_r_nb);
     free(idx_mrg); free(idx_el);
     free(idx_remd);
     

}

void schr_read(unsigned int *nds_td, unsigned int nds_n,double th_nb_koef,unsigned int m, unsigned int n, unsigned int *row, unsigned int *col, double* val,unsigned int len){
	unsigned char fl_nm[64];
	sprintf((char*)fl_nm,flnm_str,m,n,th_nb_koef);
	FILE* fp=fopen((char*)fl_nm,"rb");
	if(fp==NULL){
		printf("\nfp==NULL");
		return;
	
	}
	unsigned int max_edges=(col[len-1]*(col[len-1]-1));//maximun node number is in the end position of cl array; in fully connected network;
	unsigned int fb_sz=128; unsigned char fl_buff[fb_sz];
	unsigned int* nds_td0=(unsigned int*) malloc(nds_n*sizeof(unsigned int));
	unsigned int nds_n0;
	if(fscanf(fp,"mode1: %u",&nds_n0)==1){
	     if(nds_n==nds_n0){
	          for(unsigned int i=0;i<nds_n0;i++){
	               fscanf(fp,"%u",&(nds_td0[i]));
	               if(nds_td0[i]!=nds_td[i]){
	                    printf("nds_td0[i]!=nds_td[i]");
	                    goto cl_ret;
	                    //return;
	               }
	          }
	          unsigned int len0=len,cnt_remd_0,cnt_remd_1,idx_tmp0,idx_res;
	          unsigned int* rw0=(unsigned int*) malloc(max_edges*sizeof(unsigned int));//TODO (xn);
	          unsigned int* cl0=(unsigned int*) malloc(max_edges*sizeof(unsigned int));
	          double* vl0=(double*) malloc(max_edges*sizeof(double));
	          unsigned int* idx_mrg=(unsigned int*) malloc(2*sizeof(unsigned int)); unsigned int idx_mrg_n=0;
	          //unsigned int* idx_el=(unsigned int*) malloc(len0*(len0-1)*sizeof(unsigned int)); unsigned int idx_el_n=0;
	          unsigned int* cl_buff=(unsigned int*) malloc(max_edges*sizeof(unsigned int));
	          unsigned int* rw_buff=(unsigned int*) malloc(max_edges*sizeof(unsigned int));
	          double* vl_buff=(double*) malloc(max_edges*sizeof(double));
	          unsigned int* idx_r_td=(unsigned int*) malloc(2*nds_n0*sizeof(unsigned int));
	          unsigned int* mesh_bdrs=(unsigned int*) malloc(nds_n0*sizeof(unsigned int));
	          unsigned int* cnt_tmps=(unsigned int*) malloc(nds_n0*sizeof(unsigned int));
	          double* cnds_sum=(double*) malloc(nds_n0*sizeof(double));
	          
	        
	          for(unsigned int i=0;i<len0;i++){
	               cl0[i]=col[i]; rw0[i]=row[i]; vl0[i]=val[i];
	          }
	          //fgets((char*)fl_buff,fb_sz,fp);//read the '\n' symbol that is left in the line;
	          fgets((char*)fl_buff,fb_sz,fp);
	          while(strncmp("\ncl_buff, rw_buff (remdr):", (char*)fl_buff, strlen("cl_buff, rw_buff (remdr):"))==0){
	               sscanf((char*)fl_buff,"cl_buff, rw_buff (remdr): %u to %u",&cnt_remd_0,&cnt_remd_1);
	               for(unsigned int i=cnt_remd_0;i<cnt_remd_1;i++){
	                    fscanf(fp,"%u",&idx_tmp0);
	                    cl_buff[i]=cl0[idx_tmp0]; rw_buff[i]=rw0[idx_tmp0];
	                    vl_buff[i]=val[idx_tmp0];
	                    
	               
	               }
	               unsigned int qq=fscanf(fp,"\ncurr_nds_n (idx_r_td,mesh_bdrs): %u",&nds_n0); qq++;
	               for(unsigned int i=0;i<nds_n0;i++){
	                    fscanf(fp,"%u %u %u",&(idx_r_td[2*i]),&(idx_r_td[2*i+1]),&(mesh_bdrs[i]));
	                    
	                    
	               }
	                #ifdef _OPENMP
                         unsigned int num_threads_val=(nds_n0<(unsigned int) omp_get_max_active_levels())?(nds_n0):((unsigned int)omp_get_max_active_levels());
                    #endif
                    #pragma omp parallel for num_threads(num_threads_val)
                    for(unsigned int i=0;i<nds_n0;i++){
                         cnds_sum[i]=0;
                         for(unsigned int j=idx_r_td[2*i];j<idx_r_td[2*i+1];j++){
                              cnds_sum[i]+=vl0[j];
                         
                         }
                         for(unsigned int j=idx_r_td[2*i];j<idx_r_td[2*i+1]-1;j++){
                              for(unsigned int k=j+1;k<idx_r_td[2*i+1];k++){
                                   cnt_tmps[i]=mesh_bdrs[i]+(idx_r_td[2*i+1]-idx_r_td[2*i]-1)*(j-idx_r_td[2*i])+(k-idx_r_td[2*i]-1);
                                   cl_buff[cnt_tmps[i]]=rw0[j];//column number smaller than row number;
                                   rw_buff[cnt_tmps[i]]=rw0[k];
                                   vl_buff[cnt_tmps[i]]=1./((1./vl0[j])*(1./vl0[k])*cnds_sum[i]);
                                   cnt_tmps[i]=mesh_bdrs[i]+(idx_r_td[2*i+1]-idx_r_td[2*i]-1)*(k-idx_r_td[2*i])+(j-idx_r_td[2*i]);
                                   cl_buff[cnt_tmps[i]]=rw0[k];//column number bigger than row number;
                                   rw_buff[cnt_tmps[i]]=rw0[j];
                                   vl_buff[cnt_tmps[i]]=1./((1./vl0[j])*(1./vl0[k])*cnds_sum[i]);
                                   
                              
                              }
                         
                         }
                         
                    }
	               /////////////////////////////////////////////////
	               fscanf(fp,"\nmerges: %u",&idx_mrg_n);
	               for(unsigned int i=0;i<idx_mrg_n;i++){
	                    fscanf(fp,"%u %u",&(idx_mrg[0]),&(idx_mrg[1]));
	                    vl0[idx_mrg[0]]+=vl0[idx_mrg[1]];
	               }
	               fscanf(fp,"\nindexes: %u",&len0);
	               for(unsigned int i=0;i<len0;i++){
	                    fscanf(fp,"%u",&idx_res);
	                    cl0[i]=cl_buff[idx_res]; rw0[i]=rw_buff[idx_res];
	                    vl0[i]=vl_buff[idx_res];
	                    
	               
	               }
	               fgets((char*)fl_buff,fb_sz,fp);
	          
	          }
	          
	          
	          /*
	          
     for(unsigned int i=0;i<*nds_n;i++){
          cnds_sum[i]=0;
          for(unsigned int j=idx_r[2*i];j<idx_r[2*i+1];j++){
               cnds_sum[i]+=val[j];
               
               add from idx_r[2*i] to idx_r[2*i+1];
               
          
          }
          for(unsigned int j=idx_r[2*i];j<idx_r[2*i+1]-1;j++){
               for(unsigned int k=j+1;k<idx_r[2*i+1];k++){
                    cnt_tmps[i]=mesh_bdrs[i]+(idx_r[2*i+1]-idx_r[2*i]-1)*(j-idx_r[2*i])+(k-idx_r[2*i]-1);
                    cl_buff[cnt_tmps[i]]=row[j];//column number smaller than row number;
                    rw_buff[cnt_tmps[i]]=row[k];
                    vl_buff[cnt_tmps[i]]=1./((1./val[j])*(1./val[k])*cnds_sum[i]);
                    
                    rec: cnt_tmps[i], j,k
                    
                    cnt_tmps[i]=mesh_bdrs[i]+(idx_r[2*i+1]-idx_r[2*i]-1)*(k-idx_r[2*i])+(j-idx_r[2*i]);
                    cl_buff[cnt_tmps[i]]=row[k];//column number bigger than row number;
                    rw_buff[cnt_tmps[i]]=row[j];
                    vl_buff[cnt_tmps[i]]=1./((1./val[j])*(1./val[k])*cnds_sum[i]);
                    
                    rec: cnt_tmps[i],
                    
               
               }
          
          }
          
     }
     
     */
	          
	          
	          
	          
	          
	          
	          
	          /*
	               fgets((char*)fl_buff,fb_sz,fp);
	               if(feof(fp)){
	                    break;
	               }
	               else{
	                    if(strncmp("merges:", (char*)fl_buff, strlen("merges:"))==0){
	                          
	                         
	                         
	                    
	                    }
	                    else if(strncmp("nds_n1:", (char*)fl_buff, strlen("nds_n1:"))==0){
	                    
	                    }
	                    else{
	                         printf("\nno_match_found");
	                         break;
	                    
	                    
	                    }
	               
	               }
	               
	               if(fscanf(fp,"mode1: %u",&nds_n0)==0){
	               }
	               else if(fscanf(fp,"mode1: %u",&nds_n0)==0){
	               }
	               else if(feof(fp)){
	               
	               }
	               
	               
	               break;
	               //return;
	          
	          
	          }*/
	          
	          
	          
	          
	          
	          
	          
	          
	          
	     }
	     else{
	          printf("\nnds_n");
	          //return;
	          
	     }
	     
	}
	else{
	     printf("\nmode1 not detected");
	     //return;
	     
	}
	
	
	
	
	
	
	
	/*unsigned char* bf_r_addr=NULL, *token=NULL;//fgets reads at maximum fb_sz-1 characters from string, and 2*300*300 as maximum node number, shoud fit in padd_sz characters;
	if(strcmp("mode1:\n",fgets((char*)fl_buff,fb_sz,fp))==0){
		while(1){
			unsigned int i=0;
			unsigned int fl_str_len=strlen(fgets((char*)fl_buff,fb_sz,fp));
			while(fl_buff[fl_str_len-1]!='\n'){
				for(int j=fl_str_len-1;j>=0;j--){
					if(fl_buff[j]==' '){
						bf_r_addr=fl_buff+j;//record address where last complete token ends;
						fl_str_len=fl_str_len-j-1;//how many symbols left after;
						token=(unsigned char*)strtok((char*)fl_buff," ");
						break;
						
					}
					
				}
				while((token+strlen((char*)token))<bf_r_addr){//read tokens that are in string;
					if(nds_td0[i]!=(unsigned int)(atoi((char*)token))){
						printf("\n(nds_td0[i]!=(unsigned int)(atoi(token)))");
						return;
					
					}
					i++;
					token=(unsigned char*)strtok(NULL," ");
				
				}
				if((token+strlen((char*)token))==bf_r_addr){
					if(nds_td0[i]!=(unsigned int)(atoi((char*)token))){//read last complete token;
						printf("\n(nds_td0[i]!=(unsigned int)(atoi(token)))");
						return;
					
					}
					i++;
					for(unsigned int j=0;j<fl_str_len;j++){
						fl_buff[j]=*(bf_r_addr+j+1);//copy remaining characters to the beginning of the string;
					
					}
					fl_str_len=strlen(fgets((char*)(fl_buff+fl_str_len),fb_sz,fp));
					if(feof(fp)==0){
						break;
					
					}
				
				}
			
			
			}
		
		}
	
	}
	else{
		printf("mode1 (strcmp)\n");
	
	}*/

cl_ret:
	fclose(fp);
	
	free(nds_td0);
	

}

void schr_rec_free(unsigned int** tot_mem_buff, unsigned int mem_buff_cnt,unsigned char* fl_nm){
     for(unsigned int i=0;i<mem_buff_cnt;i++){
          free(tot_mem_buff[i]);
     
     }
     free(tot_mem_buff);
     free(fl_nm);

}
