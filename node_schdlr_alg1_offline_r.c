































#include<stdio.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>





void mode_1_alg_offline_r(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,double** vl, unsigned int len,unsigned int *ln, unsigned int *nds_td, unsigned int nds_n, double thr_koef, unsigned char out_fl, mode_3_param* mode3_inp,unsigned int m, unsigned int n){
     (void)rw; (void)cl; (void)vl; (void)ln; (void)nds_td; (void)nds_n; (void)thr_koef; (void)out_fl; (void)mode3_inp;
     struct timespec curr_time; long long unsigned int tick,dt_time;
     unsigned char fl_omp=0;
     unsigned int buff_sz=64; char ch_buff[buff_sz]; char buff_rd[buff_sz];
     snprintf(ch_buff,buff_sz,"mode1_%ux%u",m,n);
     int fd=open(ch_buff, O_RDONLY);
     snprintf(ch_buff,buff_sz,"mode1:\n");
     read(fd,buff_rd,strlen(ch_buff));
     if(strncmp(buff_rd,ch_buff,strlen(ch_buff))!=0){//'\0' was not writen in writer function;
          printf("%s!=%s\n",buff_rd,ch_buff);
          return;
     }
     
     unsigned int max_nds=col[len-1];
     double* nd_arr=(double*) calloc(((size_t)max_nds*(max_nds-1))/2,sizeof(double));//calloc should initialize pointers to NULL;
     #define arr_ij(i,j) (((((size_t)(i)-1)*(2*(size_t)(max_nds)-(size_t)(i)))/2)+((size_t)(j)-(size_t)(i)-1))
     //if maximum number of nodes is max_nds, then first node can have at maximun max_nds-1 neighbours, second -- max_nds-2 (since graph is undirected, and, hence adjacency matrix is simmetric, if there is an edge between node 1 and node 2, it should be recorded as edge of node 1 and its neighbour node 2, no need to repeate same value for node 2 and its neighbour node 1), n-th -- max_nds-n, (max_nds-1)-th -- 1 (only node number max_nds itself). Total number of doubles that is needed is sum_{i=1}^{max_nds-1} (max_nds-i)=(max_nds*(max_nds-1))/2;
     //to get index where edge of node i and its neighbour node j should be recoreded (assuming i<j) we need to sum number of edges that precede node i-th neighbours and add position of jth neighbur; before i-th node there are sum_{k=1}^{i-1} (max_nds-k)=max_nds*(i-1)-(i*(i-1))/2=(i-1)*(2*max_nds-i)/2=idx_base edges recorded, edge with (i+1)-th neighbour should be recorded at this index, (i+2)-th at idx_base+1 and so on, hence edge with j_th node should be recorded at idx_base+(j-i-1);
     
     for(unsigned int i=0;i<len;i++){
          if(col[i]<row[i]){
               nd_arr[arr_ij(col[i],row[i])]=val[i];
          }
     }
     unsigned int nds_n0;
     unsigned int* offsets=NULL; unsigned int ofst_n=0;
     unsigned int rd_buff_sz=((fl_omp==0)?1:omp_get_max_threads());
     unsigned int** rd_dt=(unsigned int**)calloc(rd_buff_sz,sizeof(unsigned int*));//allocate with calloc, so that free() on unused buffer is successful (free(null_ptr));
     unsigned int* rd_dt_cap=(unsigned int*)calloc(rd_buff_sz,sizeof(unsigned int));
     off_t curr_fl_pos;
     while(1){
          if(read(fd,&nds_n0,sizeof(unsigned int))==0 || nds_n0<=1){//probably reached end of file;
               break;
          }
          
          if(ofst_n<(nds_n0+1)){
               if(ofst_n>0){
                    free(offsets);//reallocate the buffer, but no need to copy old data;
               }
               offsets=(unsigned int*) malloc((nds_n0+1)*sizeof(unsigned int));//(nds_n0+0) because offsets in files need to hold offests for each block plus the end offset (to determine size); first offset is 0;
               ofst_n=nds_n0+1;
               
          }
          read(fd,offsets,(nds_n0+1)*sizeof(unsigned int));
          
          curr_fl_pos=lseek(fd, 0, SEEK_CUR);
          
          for(unsigned int i=0;i<nds_n0;i++){
clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
               unsigned int thr_id=omp_get_thread_num();
               if(rd_dt_cap[thr_id]<(offsets[i+1]-offsets[i])){
                    if(rd_dt_cap[thr_id]>0){
                         free(rd_dt[thr_id]);
                    }
                    rd_dt_cap[thr_id]=2*(offsets[i+1]-offsets[i]);//icreasing capacity to double of what is needed to avoid frequent allocations;
                    rd_dt[thr_id]=(unsigned int*) malloc(rd_dt_cap[thr_id]);
               }
clock_gettime(CLOCK_MONOTONIC,&curr_time);
dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
printf("\nrd_dt[thr_id] realloc=%llu",dt_time);
clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;               
               pread(fd, rd_dt[thr_id],(offsets[i+1]-offsets[i]),curr_fl_pos+offsets[i]);
clock_gettime(CLOCK_MONOTONIC,&curr_time);
dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
printf(",rd_dt[thr_id] read=%llu",dt_time);
               unsigned int* thr_data=rd_dt[thr_id];
               unsigned int neighb_sz=thr_data[0];
               double sum=0;
clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
               for(unsigned int j=0;j<neighb_sz;j++){
                    sum+=nd_arr[arr_ij(thr_data[2*j+1],thr_data[2*j+2])];
               
               }
               clock_gettime(CLOCK_MONOTONIC,&curr_time);
dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
printf(",sum calc=%llu",dt_time);
clock_gettime(CLOCK_MONOTONIC,&curr_time); tick=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
               unsigned int curr_nb_ofst=2*neighb_sz+1;
               for(unsigned int j=0;j<neighb_sz-1;j++){
                    for(unsigned int k=j+1;k<neighb_sz;k++){
                         nd_arr[arr_ij(thr_data[curr_nb_ofst+4],thr_data[curr_nb_ofst+5])] += (nd_arr[arr_ij(thr_data[curr_nb_ofst],thr_data[curr_nb_ofst+1])] * nd_arr[arr_ij(thr_data[curr_nb_ofst+2],thr_data[curr_nb_ofst+3])])/sum;
                         curr_nb_ofst+=6;
                    }
                    
               
               }
clock_gettime(CLOCK_MONOTONIC,&curr_time);
dt_time=curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec-tick;
printf(",edge calc=%llu",dt_time);
               
          
          }
          lseek(fd, offsets[nds_n0], SEEK_CUR);
     
     }
     
     
     /*if(out_fl==1 && nds_n0<=1){//ouput after mode1
          snprintf(ch_buff,buff_sz,"mode1_out:\n");
          read(fd,buff_rd,strlen(ch_buff));
          if(strncmp(buff_rd,ch_buff,strlen(ch_buff))!=0){//'\0' was not writen in writer function;
               printf("%s!=%s\n",buff_rd,ch_buff);
               return;
          }
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
          //mode_2_alg(node_arr,nds_td2,nds_n2, node_hd,max_nds,rw,cl,vl,ln,out_fl,mode3_inp,0,NULL);
          //free(nds_td0);
          //free(nds_td_rem);
          //free(node_arr);
          
          
     
     
     }*/
     
     
     
     
     free(nd_arr);
     for(unsigned int i=0;i<rd_buff_sz;i++){
          free(rd_dt[i]);
     
     }
     free(rd_dt);
     free(rd_dt_cap);
     free(offsets);
     
     close(fd);




}

