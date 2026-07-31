































#include<stdio.h>               
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>





void mode_1_alg_offline_r(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,double** vl, unsigned int len,unsigned int *ln, unsigned char out_fl, mode_3_param* mode3_inp,unsigned int m, unsigned int n){
     //TODO OpenMP, nmap, run in Matlab;  
     //struct timespec curr_time; long long unsigned int tick,dt_time;
     //unsigned char fl_omp=0;
     unsigned int buff_sz=64; char ch_buff[buff_sz]; char buff_rd[buff_sz];
     snprintf(ch_buff,buff_sz,"mode1_%ux%u",m,n);
     int fd=open(ch_buff, O_RDONLY);
     snprintf(ch_buff,buff_sz,"mode1:\n");
     read(fd,buff_rd,strlen(ch_buff));
     if(strncmp(buff_rd,ch_buff,strlen(ch_buff))!=0){//'\0' was not writen in writer function;
          printf("%s!=%s\n",buff_rd,ch_buff);
          close(fd);
          return;
     }
     unsigned int max_offst_elm;
     read(fd,&max_offst_elm,sizeof(unsigned int));
     unsigned int max_nds=col[len-1];
     double* nd_arr=(double*) calloc(((size_t)max_nds*(max_nds-1))/2,sizeof(double));//calloc should initialize pointers to NULL;
     size_t* arr_i=(size_t*) malloc((max_nds)*sizeof(size_t));//these are precalculated indexes for node n1 (smaller one) in nd_arr; to get full index for adge between nodes n1 and n2 (n1<n2), arr_i[n1] should be added apropriate shift (depending on value of n2); we need calculate them for max_nds-1 nodes, but to simplify interface (node numbers start from 1 in Matlab indexing) we store precalculated partial indexes in indexes of this array corresponding to node number (for node n1, its precalculated partial index should be stored in arr_i[n1]), arr_i[0] is unused;
     for(unsigned int i=1;i<max_nds;i++){
          arr_i[i]=((i-1)*(2*(size_t)max_nds-i))/2;
          
     
     }
     #define arr_ij(i,j) (arr_i[i]+(size_t)(j)-(i)-1)
     //if maximum number of nodes is max_nds, then first node can have at maximun max_nds-1 neighbours, second -- max_nds-2 (since graph is undirected, and, hence adjacency matrix is simmetric, if there is an edge between node 1 and node 2, it should be recorded as edge of node 1 and its neighbour node 2, no need to repeate same value for node 2 and its neighbour node 1), n-th -- max_nds-n, (max_nds-1)-th -- 1 (only node number max_nds itself). Total number of doubles that is needed is sum_{i=1}^{max_nds-1} (max_nds-i)=(max_nds*(max_nds-1))/2;
     //to get index where edge of node i and its neighbour node j should be recoreded (assuming i<j) we need to sum number of edges that precede node i-th neighbours and add position of jth neighbur; before i-th node there are sum_{k=1}^{i-1} (max_nds-k)=max_nds*(i-1)-(i*(i-1))/2=(i-1)*(2*max_nds-i)/2=idx_base edges recorded, edge with (i+1)-th neighbour should be recorded at this index, (i+2)-th at idx_base+1 and so on, hence edge with j_th node should be recorded at idx_base+(j-i-1);
     for(unsigned int i=0;i<len;i++){
          if(col[i]<row[i]){
               nd_arr[arr_ij(col[i],row[i])]=val[i];
          }
     }
     unsigned int nds_n0;
     unsigned int* offsets=NULL; unsigned int ofst_n=0;
     unsigned int rd_buff_sz=((m<128 && n<128)?1:omp_get_max_threads());//run OpenMP for mode 1 if(m>=128 || n>=128);
     unsigned int* rd_dt=(unsigned int*)calloc(rd_buff_sz*max_offst_elm,sizeof(unsigned int));//allocate with calloc, so that free() on unused buffer is successful (free(null_ptr));
     
     
     off_t curr_fl_pos;
     while(1){
          if(read(fd,&nds_n0,sizeof(unsigned int))==0 || nds_n0==0){//probably reached end of file, or transit to next mode, or to the ouput;
               break;
          }
          if(ofst_n<(nds_n0+1)){//should be triggered one (first iteration) when the indipendent set is the biggest; but occasionally indidipendent set might be several elements bigger (but as graph gets dense, indipendent set naturally decreases) on the next iteration, so it is here for safety;
               if(ofst_n>0){
                    free(offsets);//reallocate the buffer, but no need to copy old data;
               }
               offsets=(unsigned int*) malloc((nds_n0+1)*sizeof(unsigned int));//(nds_n0+0) because offsets in files need to hold offests for each block plus the end offset (to determine size); first offset is 0;
               ofst_n=nds_n0+1;
               
          }
          read(fd,offsets,(nds_n0+1)*sizeof(unsigned int));
          curr_fl_pos=lseek(fd, 0, SEEK_CUR);
          
          #pragma omp parallel for schedule(static) num_threads((nds_n0<=(unsigned int)omp_get_max_threads())?nds_n0:(unsigned int)omp_get_max_threads()) if(m>=128 || n>=128)
          for(unsigned int i=0;i<nds_n0;i++){
//clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;     */          
               pread(fd, (rd_dt+omp_get_thread_num()*max_offst_elm),(offsets[i+1]-offsets[i]),curr_fl_pos+offsets[i]);
/*clock_gettime(CLOCK_MONOTONIC,&curr_time);
dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
printf(",rd_dt[thr_id] read=%llu",dt_time);*/
               unsigned int* thr_data=rd_dt+omp_get_thread_num()*max_offst_elm;
               unsigned int neighb_sz=thr_data[0];
               unsigned int n0=thr_data[1];//n0 is a number of pivot node;
               unsigned int n1,n2,n1_min,n1_max,n2_min,n2_max;
               double sum_inv=0; double scl_mult=0;
//clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
               for(unsigned int j=0;j<neighb_sz;j++){
                    n1=min(n0,thr_data[j+2]);
                    n2=max(n0,thr_data[j+2]);
                    sum_inv+=nd_arr[arr_ij(n1,n2)];
               
               }
               sum_inv=1/sum_inv;//division is more expensive, so store inverse and multiply by it in the loop;
/*clock_gettime(CLOCK_MONOTONIC,&curr_time);
dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
printf(",sum calc=%llu",dt_time);
clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;*/
               //unsigned int curr_nb_ofst=2*neighb_sz+1;
               for(unsigned int j=0;j<neighb_sz-1;j++){
                    n1=thr_data[j+2]; n1_min=min(n0,n1); n1_max=max(n0,n1);
                    scl_mult=nd_arr[arr_ij(n1_min,n1_max)]*sum_inv;
                    for(unsigned int k=j+1;k<neighb_sz;k++){
                         n2=thr_data[k+2];
                         n2_min=min(n0,n2); n2_max=max(n0,n2);
                         
                         nd_arr[arr_ij(n1,n2)] += scl_mult * nd_arr[arr_ij(n2_min,n2_max)];
                         //curr_nb_ofst+=6;
                    }
                    
               
               }
//clock_gettime(CLOCK_MONOTONIC,&curr_time);
//dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;

//printf(",edge calc=%llu",dt_time);
               
          
          }
          lseek(fd, offsets[nds_n0], SEEK_CUR);
     
     }
     
     
     if(out_fl==1){//ouput after mode1;
          //to improve speed, and minimize output file size, transitions are signified with unisigned int val that is equal to 0 (nds_n0, or node->n), not char* buffers; it is now responsibility of a reader to read correct file (with proper setting of out_fl);
          /*snprintf(ch_buff,buff_sz,"mode1_out:\n");
          read(fd,buff_rd,strlen(ch_buff));
          if(strncmp(buff_rd,ch_buff,strlen(ch_buff))!=0){//'\0' was not writen in writer function;
               printf("%s!=%s\n",buff_rd,ch_buff);
               return;
          }*/
          read(fd,ln,sizeof(unsigned int));
          (*cl)=(unsigned int*) malloc((*ln)*sizeof(unsigned int));
          (*rw)=(unsigned int*) malloc((*ln)*sizeof(unsigned int));
          (*vl)=(double*) malloc((*ln)*sizeof(double));
          unsigned int n1,n2;
          unsigned int *out_buff=(unsigned int*) malloc((2*(*ln))*sizeof(unsigned int));
          read(fd,out_buff,(2*(*ln))*sizeof(unsigned int));
          for(unsigned int i=0;i<*ln;i++){
               (*cl)[i]=out_buff[2*i];
               (*rw)[i]=out_buff[2*i+1];
               n1=min((*cl)[i],(*rw)[i]); n2=max((*cl)[i],(*rw)[i]);
               (*vl)[i]=nd_arr[arr_ij(n1,n2)];
          
          }
          free(out_buff);
          
          
          
     }
     else{//continue to mode2
          out_fl--;
          unsigned int curr_n=0;
          while(read(fd,&curr_n,sizeof(unsigned int))!=0 && curr_n!=0){
               read(fd,rd_dt,(curr_n+1)*sizeof(unsigned int));
               unsigned int n0=rd_dt[0];
               unsigned int n1,n2,n1_min,n1_max,n2_min,n2_max;
               double sum_inv=0; double scl_mult=0;
               for(unsigned int j=0;j<curr_n;j++){
                    n1=min(n0,rd_dt[j+1]);
                    n2=max(n0,rd_dt[j+1]);
                    sum_inv+=nd_arr[arr_ij(n1,n2)];
               
               }
               sum_inv=1/sum_inv;
               for(unsigned int j=0;j<curr_n-1;j++){
                    n1=rd_dt[j+1];
                    n1_min=min(n0,n1); n1_max=max(n0,n1);
                    scl_mult=nd_arr[arr_ij(n1_min,n1_max)]*sum_inv;
                    for(unsigned int k=j+1;k<curr_n;k++){
                         n2=rd_dt[k+1];//since neighbour listst are sorted, smaller index means smaller node number, so n1<n2 by choice of n1 and n2;
                         n2_min=min(n0,n2); n2_max=max(n0,n2);
                         
                         nd_arr[arr_ij(n1,n2)] += scl_mult*nd_arr[arr_ij(n2_min,n2_max)];
                    }
                    
               
               }
               
          
          }
          if(out_fl==1){
               read(fd,ln,sizeof(unsigned int));
               (*cl)=(unsigned int*) malloc((*ln)*sizeof(unsigned int));
               (*rw)=(unsigned int*) malloc((*ln)*sizeof(unsigned int));
               (*vl)=(double*) malloc((*ln)*sizeof(double));
               unsigned int n1,n2;
               unsigned int *out_buff=(unsigned int*) malloc((2*(*ln))*sizeof(unsigned int));
               read(fd,out_buff,(2*(*ln))*sizeof(unsigned int));
               for(unsigned int i=0;i<*ln;i++){
                    (*cl)[i]=out_buff[2*i];
                    (*rw)[i]=out_buff[2*i+1];
                    n1=min((*cl)[i],(*rw)[i]); n2=max((*cl)[i],(*rw)[i]);
                    (*vl)[i]=nd_arr[arr_ij(n1,n2)];
               
               }
               free(out_buff);
          
          }
          else{//mode 3
               unsigned int buff_sz=64; char ch_buff[buff_sz];
               snprintf(ch_buff,buff_sz,"more_fold_%ux%u_%u",m,n,mode3_inp->max_m_sz);
               int fd0=open(ch_buff, O_RDONLY);
               (*(mode3_inp->n_th))=((mode3_inp->tgt_n1)%(mode3_inp->max_m_sz)==0)?(mode3_inp->tgt_n1)/(mode3_inp->max_m_sz):((mode3_inp->tgt_n1)/(mode3_inp->max_m_sz)+1);
               unsigned int ui_n_th=(*(mode3_inp->n_th));
               unsigned int ***rw=mode3_inp->rw; unsigned int ***cl=mode3_inp->cl; double ***vl=mode3_inp->vl;
               (*rw)=(unsigned int**) malloc(ui_n_th*sizeof(unsigned int*));
               (*cl)=(unsigned int**) malloc(ui_n_th*sizeof(unsigned int*));
               (*vl)=(double**) malloc(ui_n_th*sizeof(double*));
               (*(mode3_inp->ln))=(unsigned int*) malloc(ui_n_th*sizeof(unsigned int));
               unsigned int shift;
               read(fd0,&shift,sizeof(unsigned int));//node_hd->num-1;
               unsigned int max_nds0=max_nds-shift;
               size_t* arr_i0=(size_t*) malloc(max_nds0*sizeof(size_t));//precalculated indexes for smaller subgraphs (after mode 1 and mode 2); input is just node number (start with 1 in Matlab indexing), so arr_i0[0] is not used (to simplify call, even if need to allocate extra memory);
               for(unsigned int i=1;i<max_nds0;i++){
                    arr_i0[i]=((i-1)*(2*(size_t)max_nds0-i))/2;
               }
               #define arr0_ij(i,j) (arr_i0[i]+(size_t)(j)-(i)-1)
               off_t* offsets0=(off_t*) malloc((ui_n_th+1)*sizeof(off_t));
               read(fd0,offsets0,(ui_n_th+1)*sizeof(off_t));
               double** nd_arr00=(double**) malloc(ui_n_th*sizeof(double*));
               #pragma omp parallel for schedule(static) num_threads((ui_n_th<=(unsigned int)omp_get_max_threads())?ui_n_th:(unsigned int)omp_get_max_threads()) if((m>128 || n>128)||(mode3_inp->max_m_sz>=32 && ui_n_th>=4))
               for(unsigned int i=0;i<ui_n_th;i++){
                    //unsigned int buff_n=(offsets0[i+1]-offsets0[i])/sizeof(unsigned int);
                    unsigned int* m3_buff=(unsigned int*) malloc(offsets0[i+1]-offsets0[i]);//TODO check if not off by 1;
                    pread(fd0,m3_buff,offsets0[i+1]-offsets0[i],offsets0[i]);
                    unsigned int cnt_j=0;
                    unsigned int neighb_sz=m3_buff[cnt_j];
                    unsigned int n0,n1,n2,n1_min,n1_max,n2_min,n2_max;
                    nd_arr00[i]=(double*) malloc((((size_t)max_nds0*(max_nds0-1))/2)*sizeof(double));
                    for(unsigned int j=1;j<max_nds0;j++){
                         for(unsigned int k=j+1;k<(max_nds0+1);k++){
                              nd_arr00[i][arr0_ij(j,k)]=nd_arr[arr_ij((j+shift),(k+shift))];
                         
                         }
                    
                    }
                    double sum_inv=0; double scl_mult=0;
                    while(neighb_sz!=0){
                         sum_inv=0;
                         n0=m3_buff[cnt_j+1];
                         for(unsigned int k=0;k<neighb_sz;k++){
                              n1=min(n0,m3_buff[cnt_j+2+k]);
                              n2=max(n0,m3_buff[cnt_j+2+k]);
                              sum_inv+=nd_arr00[i][arr0_ij(n1,n2)];
                         
                         }
                         sum_inv=1/sum_inv;
                         for(unsigned int k=0;k<neighb_sz-1;k++){
                              n1=m3_buff[cnt_j+2+k]; n1_min=min(n0,n1); n1_max=max(n0,n1);
                              scl_mult=nd_arr00[i][arr0_ij(n1_min,n1_max)]*sum_inv;
                              for(unsigned int ll=k+1;ll<neighb_sz;ll++){
                                   n2=m3_buff[cnt_j+2+ll];
                                   n2_min=min(n0,n2); n2_max=max(n0,n2);
                                   nd_arr00[i][arr0_ij(n1,n2)] += scl_mult * nd_arr00[i][arr0_ij(n2_min,n2_max)];
                              }
                              
                              
                         
                         }
                         cnt_j+=neighb_sz+2;
                         neighb_sz=m3_buff[cnt_j];
                    
                    }
                    cnt_j++;
                    unsigned int curr_len=m3_buff[cnt_j];
                    cnt_j++;
                    (*(mode3_inp->ln))[i]=curr_len;
                    (*cl)[i]=(unsigned int*) malloc(curr_len*sizeof(unsigned int));
                    (*rw)[i]=(unsigned int*) malloc(curr_len*sizeof(unsigned int));
                    (*vl)[i]=(double*) malloc(curr_len*sizeof(double));
                    for(unsigned int j=0;j<curr_len;j++){//record output;
                         n1=m3_buff[cnt_j+2*j];
                         n2=m3_buff[cnt_j+2*j+1];
                         (*cl)[i][j]=n1+shift;
                         (*rw)[i][j]=n2+shift;
                         n1_min=min(n1,n2); n2_max=max(n1,n2);
                         (*vl)[i][j]=nd_arr00[i][arr0_ij(n1_min,n2_max)];
                    
                    }
                    free(m3_buff);
                    free(nd_arr00[i]);
                    
          
               }
               free(nd_arr00);
               free(arr_i0);
               
               close(fd0);
               free(offsets0);
               //free(nds_td0);
               //free(nds_td_rem);
               //free(node_arr);
          }
          
          
     
     
     }
     
     
     free(nd_arr);
     free(rd_dt);
     free(offsets);
     free(arr_i);
     
     close(fd);




}


