
































#include<stdio.h>
#include<string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include"mode_2_alg_w.c"
void mode_1_alg_offline_w(unsigned int *row, unsigned int *col, unsigned int len, unsigned int *nds_td, unsigned int nds_n, double thr_koef, unsigned char out_fl, mode_3_param* mode3_inp,unsigned int m, unsigned int n){
     unsigned int max_nds=col[len-1];
     unsigned int min_cap=64, max_cap=max_nds-1;
     unsigned int buff_sz=64; char ch_buff[buff_sz];
     snprintf(ch_buff,buff_sz,"mode1_%ux%u",m,n);
     int fd=open(ch_buff, O_RDWR | O_CREAT | O_TRUNC,0644);
     snprintf(ch_buff,buff_sz,"mode1:\n");
     write(fd,ch_buff,strlen(ch_buff));
     node** node_arr=(node**) malloc(max_nds*sizeof(node*));
     unsigned int* nds_td0=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int* nds_td_rem=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int* nds_td2=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     struct nd_data* curr_mrg=(struct nd_data*) malloc(1*sizeof(struct nd_data));
     unsigned int* offsets; unsigned int ofst_n=0;
     unsigned int max_offst_elm=0;//find the maximum buffer need to store biggest neighbour list and additionally its size and pivot number;
     off_t curr_fl_pos=lseek(fd, 0, SEEK_CUR);
     write(fd,&max_offst_elm,sizeof(unsigned int));//write as placeholder; when calculate maximum offset size in bytes, write it at same ofset;
     //(void)val;//not really needed for write part, but used to get rid of "parameter ‘val’ set but not used" warning;
 

     
     node* node_mem=(node*) malloc(max_nds*sizeof(node));//for better cache locality;
     node* node_hd=&(node_mem[0]);
     node_hd->num=col[0]; node_hd->cap=min_cap;
     node_hd->nums=(unsigned int*)malloc(node_hd->cap*sizeof(unsigned int));
     node* curr_node=node_hd;
     node_arr[0]=curr_node;
     for(unsigned int i=col[0]+1;i<=col[len-1];i++){//assumption here is that node numbering is sequential without skipping the numbers;
          node_arr[i-1]=&(node_mem[i-1]);
          node_arr[i-1]->num=i;
          node_arr[i-1]->cap=min_cap;
          node_arr[i-1]->nums=(unsigned int*)malloc(node_arr[i-1]->cap*sizeof(unsigned int));
          node_arr[i-2]->next=node_arr[i-1];
     }
     node_arr[col[len-1]-1]->next=NULL;
     ///////////////////////////////at this point nodes should be set up;
     
     
     unsigned int curr_col=col[0];
     unsigned int i0=0;
     for(unsigned int i=1;i<len;i++){
          if(curr_col!=col[i]){
               node_arr[col[i-1]-1]->n=i-i0;
               if(node_arr[col[i-1]-1]->n>node_arr[col[i-1]-1]->cap){
                    node_arr[col[i-1]-1]->cap=min(max(node_arr[col[i-1]-1]->n,2*node_arr[col[i-1]-1]->cap),max_cap);
                    node_arr[col[i-1]-1]->nums=(unsigned int*)realloc(node_arr[col[i-1]-1]->nums,node_arr[col[i-1]-1]->cap*sizeof(unsigned int));
               }
               for(unsigned int j=i0;j<i;j++){
                    node_arr[col[i-1]-1]->nums[j-i0]=row[j];
               
               }
               i0=i;
               curr_col=col[i];
          }
          
     
     }
     node_arr[max_nds-1]->n=len-i0;
     if(node_arr[max_nds-1]->n>node_arr[max_nds-1]->cap){
          node_arr[max_nds-1]->cap=min(max(node_arr[max_nds-1]->n,2*node_arr[max_nds-1]->cap),max_cap);
          node_arr[max_nds-1]->nums=(unsigned int*)realloc(node_arr[max_nds-1]->nums,node_arr[max_nds-1]->cap*sizeof(unsigned int));
     }
     for(unsigned int j=i0;j<len;j++){
          node_arr[max_nds-1]->nums[j-i0]=row[j];
     
     }
     /////////////////////////////////////at this point edges should be set up too;
     curr_mrg->cap=node_arr[0]->cap;
     curr_mrg->nums=(unsigned int*) malloc(curr_mrg->cap*sizeof(unsigned int));
     
     
     
     
     unsigned int nds_n0=0,nds_n_rem=0,nds_n2=nds_n;
     for(unsigned int i=0;i<nds_n;i++){
          nds_td2[i]=nds_td[i];
     
     }
     double curr_thr_koef=0;
     while(1/*dbg_cnt<dbg_max*/){
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
               nds_td_rem[i]=0;
               curr_node=node_arr[i];//nds_td0[i] should be equal to i+1 if it is not zero; curr_node->num should also be equal to i+1;
               for(unsigned int j=0;j<curr_node->n;j++){
                    if(curr_node->num<curr_node->nums[j]) nds_td0[curr_node->nums[j]-1]=0;//two threads can write to same location in nds_t0 at the same time (formally, data race), but since they all write 0, it is still 0 no matter who was last to wrote it;
                    node* node_nb=node_arr[curr_node->nums[j]-1];
                    for(unsigned int k=0;k<node_nb->n;k++){
                         if(curr_node->num<node_nb->nums[k]) nds_td0[node_nb->nums[k]-1]=0;
                    
                    }
               
               }
               i++;
          }
          
     }//at this point nds_td0 and nds_td_rem should be "symmetric" in a sense that all nodes form nds_td2 should now be devided in indipendent set (nds_td0, not having common neighbour) and whatever left from nds_td2 (nds_td_rem);
     for(unsigned int i=0;i<max_nds;i++){
          if(nds_td0[i]!=0){
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
          nds_n0=0;
          write(fd,&nds_n0,sizeof(unsigned int));//use zero as delimiter to signify transition to the next mode (mode2 or to output);
          
          break;
     
     }
     /////////////////////////////////////////at this point nds_td0 (nodes to delete) and nds_td_rem (nodes left, delete on next iteration) should set up;
    
     
     unsigned int n1;
     unsigned int bf_mg_cnt=0;
     unsigned int *ui_ptr; unsigned int ui_val;
     if(ofst_n<(nds_n0+1)){
          if(ofst_n>0){
               free(offsets);//reallocate the buffer, but no need to copy old data;
          }
          offsets=(unsigned int*) malloc((nds_n0+1)*sizeof(unsigned int));//(nds_n0+0) because offsets in files need to hold offests for each block plus the end offset (to determine size); first offset is 0;
          ofst_n=nds_n0+1;
     
     }
     write(fd,&nds_n0,sizeof(unsigned int));
     offsets[0]=0;
     for(unsigned int i=1;i<=nds_n0;i++){
          unsigned int neighb_sz=node_arr[nds_td0[i-1]-1]->n;
          offsets[i]=offsets[i-1]+(neighb_sz+2)*sizeof(unsigned int);//the format: first number is number of neighbours for given node; second number is a pivot, after that stored numbers of neighbours of pivot;
          //offsets[i]=offsets[i-1]+(2*neighb_sz+1)*sizeof(unsigned int);//the format: first number is number of neighbours for given node; then 2*neib_sz are pairs of node numbers to collect sum (edges from n0 to its neibhbours);
          //offsets[i]+=(3*(neighb_sz*(neighb_sz-1)))*sizeof(unsigned int);//the format: after that there are blocks of 3 pairs of numbers. Each block is used to create or update an edge. For expample, if pivot node is n0, and we need to create an edge between its neighbours n1 and n2, (assuming n0<n1<n2), first and second pair are (n0,n1) and (n0,n2). They are used to show which edges to take. Third pair is (n1,n2) is used to show which edge is to create (update). The pairs should be ordered (smaller, bigger) for convenience for reading (since the graph is symmetrical, only half the graph is needed to be stored).
          //Total nuber of such triples is the nuber of all new edges created (neighb_sz*(neighb_sz-1)) divided by 2 (only half of graph need to be stored):  (3*2)*(neighb_sz*(neighb_sz-1))/2=(3*(neighb_sz*(neighb_sz-1)));
          if(max_offst_elm<(neighb_sz+2)) max_offst_elm=(neighb_sz+2);
          
          
     
     }
     write(fd,offsets,(nds_n0+1)*sizeof(unsigned int));
     
     for(unsigned int i=0;i<nds_n0;i++){
          node* nd_pivot=node_arr[nds_td0[i]-1];
          write(fd,&(nd_pivot->n),sizeof(unsigned int));
          write(fd,&(nd_pivot->num),sizeof(unsigned int));
          write(fd,nd_pivot->nums,(nd_pivot->n)*sizeof(unsigned int));
          /*for(unsigned int j=0;j<nd_pivot->n;j++){
               if(nd_pivot->num<nd_pivot->nums[j]){
                    write(fd,&(nd_pivot->num),sizeof(unsigned int));
                    write(fd,&(nd_pivot->nums[j]),sizeof(unsigned int));
               }
               else{
                    write(fd,&(nd_pivot->nums[j]),sizeof(unsigned int));
                    write(fd,&(nd_pivot->num),sizeof(unsigned int));
               }
          
          }//edges to calculate sum;
          */
          for(unsigned int j=0;j<nd_pivot->n;j++){//for writing we need to iterate untill j<nd_pivot->n-1,but for graph update, we still need to process (nd_pivot->n-1)-th node too;
               /*for(unsigned int k=j+1;k<nd_pivot->n;k++){
                    if(nd_pivot->num < nd_pivot->nums[j]){
                         write(fd,&(nd_pivot->num),sizeof(unsigned int));
                         write(fd,&(nd_pivot->nums[j]),sizeof(unsigned int));
                         write(fd,&(nd_pivot->num),sizeof(unsigned int));
                         write(fd,&(nd_pivot->nums[k]),sizeof(unsigned int));
                    }
                    else if(nd_pivot->num > nd_pivot->nums[j] && nd_pivot->num < nd_pivot->nums[k]){
                         write(fd,&(nd_pivot->nums[j]),sizeof(unsigned int));
                         write(fd,&(nd_pivot->num),sizeof(unsigned int));
                         write(fd,&(nd_pivot->num),sizeof(unsigned int));
                         write(fd,&(nd_pivot->nums[k]),sizeof(unsigned int));
                    }
                    else{
                         write(fd,&(nd_pivot->nums[j]),sizeof(unsigned int));
                         write(fd,&(nd_pivot->num),sizeof(unsigned int));
                         write(fd,&(nd_pivot->nums[k]),sizeof(unsigned int));
                         write(fd,&(nd_pivot->num),sizeof(unsigned int));
                    }
                    write(fd,&(nd_pivot->nums[j]),sizeof(unsigned int));
                    write(fd,&(nd_pivot->nums[k]),sizeof(unsigned int));
               }*/
               n1=nd_pivot->nums[j];
               bf_mg_cnt=0;
               if(curr_mrg->cap<min((node_arr[n1-1]->n + nd_pivot->n),max_nds-1)){
                    curr_mrg->cap=min(max(node_arr[n1-1]->n+nd_pivot->n,2*curr_mrg->cap),max_nds-1);
                    free(curr_mrg->nums);
                    curr_mrg->nums=(unsigned int*) malloc(curr_mrg->cap*sizeof(unsigned int));
               
               }
               unsigned int ll=0,qq=0;
               for(; ll<node_arr[n1-1]->n && qq<nd_pivot->n;){
                    if(nd_pivot->nums[qq]==n1){//skip node from the list of neighbourds of the node to delete that has same index as the neighbout itself (form edges with all pivot's neighbours except itself);
                         qq++; continue;
                    
                    }
                    if(node_arr[n1-1]->nums[ll]==nds_td0[i]){//for node node_td0[i] and its neighbours nums=node_arr[node_td[i]-1]->nums, each of nodes in nums has node_td0[i] as their own neighbour (if a has edge with b then b has edge with a), so skip nds_td0[i] node in its neighbour list of neighbours; no other node from nds_td0 is in neigbour list because othewise this node whould be the common node (common neigbour) for two nodes from nds_td0;
                         ll++; continue;
                    
                    }
                    if(node_arr[n1-1]->nums[ll] == nd_pivot->nums[qq]){//edge exists, parralel conductances are added;
                         curr_mrg->nums[bf_mg_cnt]=node_arr[n1-1]->nums[ll];
                         bf_mg_cnt++; ll++; qq++;
                    }
                    else if(node_arr[n1-1]->nums[ll] < nd_pivot->nums[qq]){//record old edge;
                         curr_mrg->nums[bf_mg_cnt]=node_arr[n1-1]->nums[ll];
                         bf_mg_cnt++; ll++;
                    
                    }
                    else{//calculate and recocrd new edge;
                         curr_mrg->nums[bf_mg_cnt]=nd_pivot->nums[qq];
                         bf_mg_cnt++; qq++;
                         
                    }
               
               }
               while(ll<node_arr[n1-1]->n){
                    if(node_arr[n1-1]->nums[ll]==nds_td0[i]){
                         ll++;
                    }
                    else{
                         curr_mrg->nums[bf_mg_cnt]=node_arr[n1-1]->nums[ll];
                         bf_mg_cnt++; ll++;
                    }
               
               }
               while(qq<nd_pivot->n){
                    if(nd_pivot->nums[qq]==n1){
                         qq++; continue;
                    
                    }
                    else{
                         curr_mrg->nums[bf_mg_cnt]=nd_pivot->nums[qq];
                         bf_mg_cnt++; qq++;
                    
                    }
               
               }
               ui_val=node_arr[n1-1]->cap; node_arr[n1-1]->cap=curr_mrg->cap; curr_mrg->cap=ui_val;
               ui_ptr=node_arr[n1-1]->nums; node_arr[n1-1]->nums=curr_mrg->nums; curr_mrg->nums=ui_ptr;
               node_arr[n1-1]->n=bf_mg_cnt;
               
               
          
          }
     
     }
     /////////////////////////////////////////////////////////////////////at this point star mesh transform should be completed for this iteration
     
     
     
     
     node* prev_node=NULL;curr_node=node_hd;
     for(unsigned int i=0; i<nds_n0 && curr_node!=NULL;){
          if(curr_node->num==nds_td0[i]){
               free(curr_node->nums);
               node_arr[curr_node->num-1]=NULL;
               if(prev_node==NULL){
                    node_hd=curr_node->next;
                    curr_node=node_hd;
               
               }
               else{
                    prev_node->next=curr_node->next;
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
     }
     
     if(out_fl==1){//ouput after mode1
          //to improve speed, and minimize output file size, do not signify transition with char* buffers; it is now responsibility of a reader to read correct file (with proper setting of out_fl);
          //snprintf(ch_buff,buff_sz,"mode1_out:\n");
          //write(fd,ch_buff,strlen(ch_buff));
          unsigned int ln_curr=0;
          for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
               ln_curr+=curr_node->n;
          
          }
          write(fd,&ln_curr,sizeof(unsigned int));
          for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
               for(unsigned int i=0;i<curr_node->n;i++){
                    write(fd,&(curr_node->num),sizeof(unsigned int));
                    write(fd,&(curr_node->nums[i]),sizeof(unsigned int));
                    
               
               }
          
          }
          free(nds_td0);
          free(nds_td_rem);
          node* prev_node=NULL;
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               
            
          
          }
          free(node_arr);
          
          
          
     }
     else{//continue to mode2 and mode3
          out_fl--;
          mode_2_alg_w(node_arr,nds_td2,nds_n2, node_hd,max_nds,out_fl,mode3_inp,fd,&max_offst_elm,m,n);
          free(nds_td0);
          free(nds_td_rem);
          free(node_arr);
          
          
     
     
     }
     
     
     free(node_mem);
     free(nds_td2);
     free(curr_mrg->nums);
     free(curr_mrg);
     free(offsets);
     
     pwrite(fd, &max_offst_elm,sizeof(unsigned int),curr_fl_pos);//record size fo the biggest neighbour list;
     close(fd);
     

}




